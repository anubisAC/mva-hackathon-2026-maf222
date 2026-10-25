"""Local broad drug/proxy signature comparison; never patient expression or efficacy.

Uses paired public top-gene signatures, preserving signed direction and context.
The descriptive score is the negative Pearson correlation of signed gene-list
membership (up=1, down=-1, unlisted=0) on a stated shared observed universe.
No significance or effect-size inference is made from this screening score.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import statistics
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


DATASETS = {
    "drug": "lincs-l1000-cp.gmt.gz",
    "knockout": "lincs-l1000-xpr.gmt.gz",
    "knockdown": "cm4ai.gmt.gz",
}
SOURCE = "https://perturbseqr.maayanlab.cloud/download"
CAUTION = (
    "Public BUB1B perturbations are disease-model proxies, not this patient's "
    "expression profile or the specific alleles. Directional gene-list reversal "
    "does not establish functional rescue, tolerability, or therapeutic benefit."
)


def normalized_name(value: str) -> str:
    return " ".join(str(value).casefold().split())


def read_pairs(path: Path, gene_index: dict, select_gene: str | None = None):
    """Read complete GMT universe; retain all drugs or only a specified gene.

    Integer bitsets keep the complete drug catalog tractable without retaining
    millions of duplicate Python gene strings. No aliases are inferred.
    """
    pairs = {}
    universe = 0
    counts = {"rows": 0, "duplicate_identical_rows": 0, "unpaired_contexts": 0,
              "ambiguous_direction_genes_removed": 0}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Invalid GMT row in {path}")
            base, direction = fields[0].rsplit(" ", 1)
            if direction not in {"up", "down"}:
                raise ValueError(f"Missing signed direction: {fields[0]}")
            mask = 0
            for symbol in fields[2:]:
                if symbol:
                    index = gene_index.setdefault(symbol, len(gene_index))
                    mask |= 1 << index
            universe |= mask
            counts["rows"] += 1
            term_fields = base.split("::")
            if select_gene and term_fields[0].upper() != select_gene.upper():
                continue
            pair = pairs.setdefault(base, {"term": base, "fields": term_fields,
                                           "source_identifier": fields[1]})
            if direction in pair:
                if pair[direction] != mask:
                    raise ValueError(f"Conflicting duplicate signed context: {base}")
                counts["duplicate_identical_rows"] += 1
            pair[direction] = mask
    result = []
    for pair in pairs.values():
        if "up" not in pair or "down" not in pair:
            counts["unpaired_contexts"] += 1
            continue
        overlap = pair["up"] & pair["down"]
        counts["ambiguous_direction_genes_removed"] += overlap.bit_count()
        pair["up"] &= ~overlap
        pair["down"] &= ~overlap
        result.append(pair)
    counts["paired_contexts"] = len(result)
    counts["observed_universe_genes"] = universe.bit_count()
    return result, universe, counts


def compare_signed(query: dict, drug: dict, universe: int) -> dict:
    qu, qd = query["up"] & universe, query["down"] & universe
    du, dd = drug["up"] & universe, drug["down"] & universe
    n = universe.bit_count()
    qn, dn = qu.bit_count() + qd.bit_count(), du.bit_count() + dd.bit_count()
    qs, ds = qu.bit_count() - qd.bit_count(), du.bit_count() - dd.bit_count()
    reversed_count = (qu & dd).bit_count() + (qd & du).bit_count()
    mimicked_count = (qu & du).bit_count() + (qd & dd).bit_count()
    variance = (qn - qs * qs / n) * (dn - ds * ds / n) if n else 0
    score = ((reversed_count - mimicked_count + qs * ds / n) /
             math.sqrt(variance)) if variance > 0 else None
    return {"reversal_score": score, "mimic_score": -score if score is not None else None,
            "opposite_direction_overlap": reversed_count,
            "same_direction_overlap": mimicked_count, "shared_universe_genes": n,
            "query_up_in_universe": qu.bit_count(), "query_down_in_universe": qd.bit_count(),
            "drug_up_in_universe": du.bit_count(), "drug_down_in_universe": dd.bit_count()}


def catalog_name_index(catalog: dict) -> dict:
    index = defaultdict(set)
    for drug in catalog.get("drugs", []):
        names = [drug.get("id", ""), drug.get("name", "")]
        for synonym in drug.get("synonyms", []):
            if isinstance(synonym, str):
                names.append(synonym)
            elif isinstance(synonym, dict):
                names.extend([synonym.get("molecule_synonym", ""), synonym.get("synonyms", "")])
        for name in names:
            if name:
                index[normalized_name(name)].add(drug["id"])
    return dict(index)


def match_catalog(name: str, index: dict) -> dict:
    ids = sorted(index.get(normalized_name(name), set()))
    return {"phase4_ids": ids, "mapping_status": "unique_exact" if len(ids) == 1 else
            "ambiguous_exact" if ids else "unmatched", "current_regulatory_status": "unverified"}


def aggregate_contexts(values: list[dict]) -> dict:
    """Equal weighting: drug/proxy comparison -> dose-time -> cell -> drug.

    The dose-time context key includes cell type, time and concentration. This
    avoids choosing a best dose or letting a repeatedly profiled cell dominate.
    """
    conditions = defaultdict(list)
    for row in values:
        if row["score"] is not None:
            conditions[(row["cell"], row["condition"])].append(row["score"])
    cells = defaultdict(list)
    for (cell, _), scores in conditions.items():
        cells[cell].append(statistics.median(scores))
    per_cell = {cell: statistics.median(scores) for cell, scores in sorted(cells.items())}
    scores = list(per_cell.values())
    raw_scores = [row["score"] for row in values if row["score"] is not None]
    return {"score": statistics.median(scores) if scores else None, "cell_count": len(scores),
            "positive_cell_fraction": sum(x > 0 for x in scores) / len(scores) if scores else None,
            "context_count": len(conditions), "comparison_count": len(values),
            "valid_comparison_count": len(raw_scores),
            "positive_comparison_fraction": sum(x > 0 for x in raw_scores) / len(raw_scores) if raw_scores else None,
            "min_comparison_score": min(raw_scores) if raw_scores else None,
            "max_comparison_score": max(raw_scores) if raw_scores else None,
            "per_cell_scores": per_cell, "min_cell_score": min(scores) if scores else None,
            "max_cell_score": max(scores) if scores else None}


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def download_sources(directory: Path) -> dict:
    """Acquire whole public catalogs only; never send a gene or case query."""
    directory.mkdir(parents=True, exist_ok=True)
    receipts = {}
    for key, filename in DATASETS.items():
        path = directory / filename
        url = "https://perturbseqr.s3.amazonaws.com/" + filename
        receipt_path = directory / (filename + ".receipt.json")
        if path.exists() and receipt_path.exists():
            previous = json.loads(receipt_path.read_text(encoding="utf-8"))
            if file_hash(path) != previous["sha256"]:
                raise ValueError(f"Cached file differs from its receipt: {path}")
            receipts[key] = previous
            continue
        if not path.exists():
            partial = path.with_suffix(path.suffix + ".part")
            with urllib.request.urlopen(url, timeout=120) as response, partial.open("wb") as stream:
                while block := response.read(4 * 1024 * 1024):
                    stream.write(block)
            partial.replace(path)
        receipt = {"url": url, "bytes": path.stat().st_size, "sha256": file_hash(path),
                   "receipt_recorded_at": datetime.now(timezone.utc).isoformat(),
                   "public_bulk_catalog": True, "private_query_transmitted": False}
        write_json(receipt_path, receipt)
        receipts[key] = receipt
    write_json(directory / "acquisition.json", receipts)
    return receipts


def csv_rows(path: Path, rows: list[dict]):
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, separators=(",", ":")) if isinstance(v, (list, dict))
                             else v for k, v in row.items()})


def run(signature_dir: Path, catalog_path: Path, output_dir: Path, gene: str = "BUB1B") -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    if (output_dir / "summary.json").exists():
        raise FileExistsError("Use a new output directory to preserve the previous screen")
    started = datetime.now(timezone.utc).isoformat()
    code_hash = file_hash(Path(__file__))
    input_paths = {key: signature_dir / filename for key, filename in DATASETS.items()}
    input_paths["catalog"] = catalog_path
    inputs = {key: {"path": str(path.resolve()), "sha256": file_hash(path)}
              for key, path in input_paths.items()}
    for key in DATASETS:
        inputs[key]["source_url"] = "https://perturbseqr.s3.amazonaws.com/" + DATASETS[key]
    index = {}
    drugs, drug_universe, drug_audit = read_pairs(input_paths["drug"], index)
    knockout, ko_universe, ko_audit = read_pairs(input_paths["knockout"], index, gene)
    knockdown, kd_universe, kd_audit = read_pairs(input_paths["knockdown"], index, gene)
    if not knockout or not knockdown:
        raise ValueError("Both registered public proxy sources must contain a paired signature")
    ko_by_cell = defaultdict(list)
    for proxy in knockout:
        ko_by_cell[proxy["fields"][1]].append(proxy)
    knockdown = [p for p in knockdown if p["fields"][1] == "KOLF2.1J" and
                 "KD" in p["fields"] and "CRISPRi" in p["fields"]]
    if not knockdown:
        raise ValueError("The registered hiPSC CRISPRi proxy is missing")
    ko_shared = drug_universe & ko_universe
    kd_shared = drug_universe & kd_universe
    mapping = catalog_name_index(json.loads(catalog_path.read_text(encoding="utf-8")))
    values = defaultdict(lambda: {"same_cell": [], "hipsc_cross_cell": []})
    source_ids = defaultdict(set)
    comparison_count = 0
    fields = ["drug", "cell", "cell_type", "time", "dose", "drug_context", "source_identifier",
              "route", "proxy_context", "proxy_cell", "reversal_score", "mimic_score",
              "opposite_direction_overlap", "same_direction_overlap", "shared_universe_genes",
              "query_up_in_universe", "query_down_in_universe", "drug_up_in_universe", "drug_down_in_universe"]
    with gzip.open(output_dir / "all_signature_comparisons.csv.gz", "wt", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for drug in drugs:
            term = drug["fields"] + [""] * 5
            name, cell, cell_type, time, dose = term[:5]
            source_ids[name].add(drug["source_identifier"])
            condition = "::".join(term[:5])
            routes = [("same_cell", ko_by_cell.get(cell, []), ko_shared),
                      ("hipsc_cross_cell", knockdown, kd_shared)]
            for route, proxies, universe in routes:
                for proxy in proxies:
                    result = compare_signed(proxy, drug, universe)
                    values[name][route].append({"cell": cell, "condition": condition,
                                               "score": result["reversal_score"]})
                    writer.writerow({"drug": name, "cell": cell, "cell_type": cell_type,
                                     "time": time, "dose": dose, "drug_context": drug["term"],
                                     "source_identifier": drug["source_identifier"], "route": route,
                                     "proxy_context": proxy["term"], "proxy_cell": proxy["fields"][1], **result})
                    comparison_count += 1
    rows = []
    for name, groups in values.items():
        same = aggregate_contexts(groups["same_cell"])
        hipsc = aggregate_contexts(groups["hipsc_cross_cell"])
        row = {"drug": name, **match_catalog(name, mapping)}
        row["source_identifiers"] = sorted(source_ids[name])
        row["source_identifier_count"] = len(source_ids[name])
        row["identity_status"] = ("single_source_identifier" if len(source_ids[name]) == 1 and
                                   "" not in source_ids[name] else "multiple_or_missing_source_identifiers_unresolved")
        row["identity_resolved_for_candidate_review"] = (row["mapping_status"] == "unique_exact" and
                                                         row["identity_status"] == "single_source_identifier")
        row.update({"same_cell_" + k: v for k, v in same.items()})
        row.update({"hipsc_cross_cell_" + k: v for k, v in hipsc.items()})
        row["multi_cell_coverage"] = same["cell_count"] >= 2
        row["positive_both_proxy_routes"] = bool(same["score"] is not None and same["score"] > 0 and
                                                   hipsc["score"] is not None and hipsc["score"] > 0)
        row["evidence_status"] = "public_proxy_directional_similarity_only"
        row["mechanism_gate"] = ("known_adverse_spindle_comparator" if normalized_name(name) == "colchicine"
                                  else "requires_mechanism_and_confounding_review")
        rows.append(row)
    rows.sort(key=lambda r: (not r["multi_cell_coverage"],
                            -(r["same_cell_score"] if r["same_cell_score"] is not None else -2), r["drug"].casefold()))
    for rank, row in enumerate(rows, 1):
        row["screen_order"] = rank
    exact = [row for row in rows if row["mapping_status"] == "unique_exact"]
    write_json(output_dir / "drug_rankings.json", rows)
    csv_rows(output_dir / "drug_rankings.csv", rows)
    csv_rows(output_dir / "phase4_exact_mapped.csv", exact)
    csv_rows(output_dir / "name_mapping_audit.csv", [{"drug": r["drug"], "mapping_status": r["mapping_status"],
                                                     "phase4_ids": r["phase4_ids"]} for r in rows])
    proxy_contexts = [{"source": source, "term": p["term"], "cell": p["fields"][1],
                       "up_genes": sorted(g for g, i in index.items() if p["up"] & 1 << i),
                       "down_genes": sorted(g for g, i in index.items() if p["down"] & 1 << i)}
                      for source, proxies in [("LINCS_XPR", knockout), ("CM4AI", knockdown)] for p in proxies]
    write_json(output_dir / "public_proxy_signatures.json", proxy_contexts)
    write_json(output_dir / "observed_background_genes.json", {
        "definition": "Intersection of genes appearing anywhere in the complete paired up/down GMT source libraries; not a complete measured-gene or genome background.",
        "same_cell": sorted(g for g, i in index.items() if ko_shared & 1 << i),
        "hipsc_cross_cell": sorted(g for g, i in index.items() if kd_shared & 1 << i)})
    summary = {"status": "complete", "gene": gene, "started_at": started,
               "completed_at": datetime.now(timezone.utc).isoformat(), "inputs": inputs,
               "code_sha256": code_hash, "source": SOURCE, "caution": CAUTION,
               "drug_audit": drug_audit, "knockout_audit": ko_audit, "knockdown_audit": kd_audit,
               "unique_drugs": len(rows), "paired_drug_signatures": len(drugs),
               "knockout_proxy_signatures": len(knockout), "knockout_proxy_cells": sorted(ko_by_cell),
               "hipsc_proxy_signatures": len(knockdown), "signature_comparisons": comparison_count,
               "drugs_with_same_cell_coverage": sum(r["same_cell_cell_count"] > 0 for r in rows),
               "drugs_with_multiple_matched_cells": sum(r["multi_cell_coverage"] for r in rows),
               "unique_exact_phase4_name_matches": len(exact),
               "ambiguous_phase4_name_matches": sum(r["mapping_status"] == "ambiguous_exact" for r in rows),
               "drug_names_with_multiple_source_ids": sum(r["source_identifier_count"] > 1 for r in rows),
               "exact_phase4_matches_with_single_source_id": sum(r["identity_resolved_for_candidate_review"] for r in rows),
               "adverse_mechanistic_comparators": [{"drug": r["drug"], "same_cell_score": r["same_cell_score"],
                                                   "hipsc_cross_cell_score": r["hipsc_cross_cell_score"]}
                                                  for r in rows if r["mechanism_gate"] == "known_adverse_spindle_comparator"],
               "same_cell_observed_universe_genes": ko_shared.bit_count(),
               "hipsc_observed_universe_genes": kd_shared.bit_count(),
               "statistical_status": "Descriptive screening scores; no p-values or significance claims. Multiple testing is not represented by these uncalibrated ranks.",
               "top_phase4_by_same_cell_order": [r["drug"] for r in exact[:10]]}
    write_json(output_dir / "summary.json", summary)
    report = ["# Broad public-proxy drug signature screen", "", CAUTION, "",
              f"Compared **{len(drugs):,} paired drug contexts** spanning **{len(rows):,} named compounds**. "
              f"{summary['drugs_with_same_cell_coverage']:,} compounds have at least one same-cell comparison; "
              f"{summary['drugs_with_multiple_matched_cells']:,} cover at least two matched cells. "
              f"{len(exact):,} names match one phase-4 ChEMBL record exactly; phase 4 does not verify current approval.", "",
              "## What the score measures", "",
              "A positive score means opposite signed gene-list membership; a negative score means mimicry. "
              "Each gene is +1 for a reported up list, -1 for down, and 0 otherwise. The score is the negative "
              "Pearson correlation of these memberships over the source libraries' shared observed-gene universe. "
              "Unlisted=0 means absent from the thresholded list, not a measured lack of expression change. "
              "Backgrounds are library-wide, not the detected-gene background of each individual assay. "
              "It is not a correlation of expression magnitudes. Full background gene lists are saved. "
              "No p-values, statistical significance, or probabilities of rescue are assigned.", "",
              "The downloaded LINCS lists are pruned to top genes from Characteristic Direction coefficients. "
              "The source method compares replicate perturbation profiles with other profiles in the same batch; "
              "these are not necessarily untreated controls. L1000 directly measures 978 genes and infers many "
              "others. These processing choices limit biological interpretation of signed overlaps.", "",
              "The primary comparison uses BUB1B CRISPR knockouts and drugs in the same cancer cell line. "
              "The independent secondary comparison uses BUB1B CRISPRi knockdown in KOLF2.1J hiPSCs against "
              "drug signatures from other cells. The gene and drug perturbation durations differ and are recorded. "
              "An acute knockout or knockdown does not model the patient's two alleles exactly.", "",
              "All tested drug doses and times contribute. Scores are collapsed by median within a cell/dose/time "
              "condition, then median across conditions within each cell, then median across cells. Cancer-line "
              "knockout signatures are pooled within condition, without selecting the best signature. Ordering places "
              "compounds with at least two same-cell contexts first, followed by descending same-cell score. "
              "This ordering is a review aid, not a validated therapeutic ranking. Raw comparison counts, "
              "positive fractions and score ranges expose within-cell and dose heterogeneity; these comparisons "
              "are correlated and their counts do not denote independent biological replication.", "",
              "## First exact phase-4 name matches in that order", "",
              "| Compound | Same-cell score | Matched cells | Positive raw comparisons | hiPSC cross-cell score | Source IDs |",
              "|---|---:|---:|---:|---:|---:|"]
    for r in exact[:20]:
        fmt = lambda x: "n/a" if x is None else f"{x:.4f}"
        report.append(f"| {r['drug']} | {fmt(r['same_cell_score'])} | {r['same_cell_cell_count']} | "
                      f"{fmt(r['same_cell_positive_comparison_fraction'])} | {fmt(r['hipsc_cross_cell_score'])} | "
                      f"{r['source_identifier_count']} |")
    report.extend(["", "## Interpretation gates", "",
                   "A compound can reverse a proliferative signature by arresting or killing cells. It can also "
                   "suppress a compensatory protective response. Every candidate therefore requires target-action "
                   "direction review, concentration/exposure review, viability controls, and chromosome-segregation "
                   "and checkpoint assays. Same-cell agreement and an independent hiPSC direction are evidence "
                   "checks, not substitutes for these tests. A surviving candidate should improve segregation "
                   "without simply stopping mitosis or removing abnormal cells.", "",
                   "Colchicine is a known adverse mechanistic comparator for spindle restoration because its "
                   "catalogued mechanism inhibits tubulin/microtubule assembly. It nevertheless returns a positive "
                   "reversal score in this screen. This is a failure of reversal alone as a rescue-selection rule, "
                   "not support for colchicine. Empirical-null calibration, cell-cycle and cytotoxicity controls, "
                   "and signed mechanism gates are required before selecting candidate interventions.", "",
                   "Drug names are matched only by case and whitespace normalized exact names or source synonyms. "
                   "Ambiguous matches and unmapped compounds remain visible. No salt stripping, inferred approval, "
                   "or structure identity confirmation is performed. Source BRD identifiers are preserved, and "
                   "names with multiple identifiers are explicitly flagged. Their named-compound scores pool "
                   "unresolved source records and are exploratory until forms/structures are reconciled. "
                   "Some research compounds and "
                   "discontinued products are present in the drug source. Drug identifiers, dose, time, cell, and "
                   "all pair comparisons are preserved in the compressed comparison table.", "",
                   f"Sources: [Perturb-Seqr full library downloads]({SOURCE}), "
                   "[LINCS data methodology](https://maayanlab.cloud/sigcom-lincs/static/tutorial/ldp3/LDP3Tutorial.html), "
                   "[CM4AI hiPSC perturbation atlas](https://www.nature.com/articles/s41587-026-03199-w). "
                   "The exact public files and hashes are recorded in summary.json.", ""])
    (output_dir / "SIGNATURE_SCREEN.md").write_text("\n".join(report), encoding="utf-8")
    if file_hash(Path(__file__)) != code_hash:
        raise RuntimeError("Source code changed during execution; preserve this run as incomplete")
    manifest = {"status": "complete", "inputs": inputs, "code_path": str(Path(__file__).resolve()),
                "code_sha256": code_hash,
                "outputs": [{"path": str(path.resolve()), "bytes": path.stat().st_size,
                             "sha256": file_hash(path)} for path in sorted(output_dir.iterdir()) if path.is_file()]}
    write_json(output_dir / "run_manifest.json", manifest)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--signatures", type=Path, default=Path("references/drug_screen/signatures"))
    parser.add_argument("--catalog", type=Path, default=Path("references/drug_screen/catalog.json"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--download", action="store_true", help="Download/cache entire public bulk libraries")
    parser.add_argument("--gene", default="BUB1B")
    args = parser.parse_args()
    if args.download:
        print(json.dumps(download_sources(args.signatures), indent=2))
    if args.output:
        print(json.dumps(run(args.signatures, args.catalog, args.output, args.gene), indent=2))
    elif not args.download:
        parser.error("Provide --output for a screen or --download for acquisition")


if __name__ == "__main__":
    main()
