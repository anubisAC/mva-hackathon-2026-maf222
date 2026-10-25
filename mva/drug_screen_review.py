"""Join broad screening evidence without inventing a composite efficacy score."""
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

from .intake import sha256
from .report import ACKNOWLEDGEMENT
from .strategy_b import baseline_fingerprints


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def build(root: Path, catalog_path: Path) -> dict:
    catalog = read(catalog_path)
    molecules = {d["id"]: d for d in catalog["drugs"]}
    targets = read(root / "targets/ranked_drugs.json")
    families = {row["id"]: row for row in targets}
    drug_rows = read(root / "signatures/drug_rankings.json")
    ss = read(root / "signatures/summary.json")
    ts = read(root / "targets/summary.json")
    matrix = []
    for row in drug_rows:
        if row["mapping_status"] != "unique_exact":
            continue
        molecule = molecules[row["phase4_ids"][0]]
        parent_id = molecule.get("parent_id") or molecule["id"]
        family = families[parent_id]
        matrix.append({"compound_name": row["drug"], "chembl_id": molecule["id"], "parent_id": parent_id,
                       "source_identifiers": row.get("source_identifiers", []),
                       "source_identifier_count": row.get("source_identifier_count"),
                       "identity_status": row.get("identity_status", "name_mapping_only"),
                       "source_screen_order": row["screen_order"], "same_cell_score": row["same_cell_score"],
                       "same_cell_count": row["same_cell_cell_count"],
                       "same_cell_positive_fraction": row["same_cell_positive_cell_fraction"],
                       "same_cell_min_score": row["same_cell_min_cell_score"],
                       "same_cell_max_score": row["same_cell_max_cell_score"],
                       "same_cell_raw_comparisons": row.get("same_cell_valid_comparison_count"),
                       "same_cell_raw_positive_fraction": row.get("same_cell_positive_comparison_fraction"),
                       "same_cell_raw_min_score": row.get("same_cell_min_comparison_score"),
                       "same_cell_raw_max_score": row.get("same_cell_max_comparison_score"),
                       "hipsc_cross_cell_score": row["hipsc_cross_cell_score"],
                       "both_proxy_scores_positive": row["positive_both_proxy_routes"],
                       "direct_target_module_overlap": family["module_overlap"],
                       "contextual_module_overlap": family["contextual_module_overlap"],
                       "target_overlap_q": family["degree_adjusted_q"],
                       "direction_flags": family["direction_flags"],
                       "signature_mechanism_gate": row.get("mechanism_gate", "requires_review"),
                       "withdrawal_flag_present": family["withdrawal_flag_present"],
                       "withdrawal_warning_present": molecule.get("withdrawal_warning_present"),
                       "black_box_flag_present": family["black_box_flag_present"],
                       "veterinary_flag_present": family["veterinary_flag_present"],
                       "current_approval": "unverified", "functional_rescue": "unmeasured",
                       "interpretation": "screening association only; no composite benefit score"})
    write(root / "evidence_matrix.json", matrix)
    with (root / "evidence_matrix.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(matrix[0]) if matrix else ["compound_name"])
        writer.writeheader()
        for row in matrix:
            writer.writerow({k: json.dumps(v) if isinstance(v, (dict, list)) else v for k, v in row.items()})
    prior = {name: next((r for r in drug_rows if r["drug"].casefold() == name), None)
             for name in ("azithromycin", "gentamicin", "amikacin", "colchicine")}
    write(root / "previous_probe_comparison.json", prior)
    selection = read(root / "targets/module_selection.json")
    summary = {"status": "completed", "catalog_records": len(catalog["drugs"]),
               "parent_drug_families": ts["drug_families"], "phase4_exact_signature_names": len(matrix),
               "distinct_mapped_parent_families": len({r["parent_id"] for r in matrix}),
               "signature_drug_names": ss["unique_drugs"], "drug_contexts": ss["paired_drug_signatures"],
               "signature_comparisons": ss["signature_comparisons"],
               "same_cell_covered_names": ss["drugs_with_same_cell_coverage"],
               "multiple_same_cell_names": ss["drugs_with_multiple_matched_cells"],
               "target_pathway_q_below_0_05": ts["degree_adjusted_q_below_0_05"],
               "validated_rescue_candidates": 0,
               "decision": "Use broad screening for hypothesis generation; calibrate nulls and cytostatic confounding before choosing new leads."}
    write(root / "summary.json", summary)
    fmt = lambda value: "unavailable" if value is None else f"{value:.5f}"
    lines = ["# Broad computational drug screening for the BUB1B hypothesis", "",
             "A broad drug screen is a better discovery design than selecting a few antibiotics from literature. This first implementation combines drug-target annotations with public gene-perturbation and drug-expression experiments. It can uncover candidate mechanisms and expose contradictory evidence. It does not simulate treatment in this patient, establish drug binding to either allele, or validate functional rescue.", "",
             "## What has been screened", "",
             "| Layer | Completed scope | Interpretation |", "| --- | --- | --- |",
             f"| Drug catalog | {len(catalog['drugs']):,} ChEMBL phase-four molecule records, {ts['drug_families']:,} parent families | Historical development-stage annotation; current human approval needs separate review |",
             f"| Target and pathway | {len(selection['selected'])} Reactome modules, {ts['module_genes_total']} union genes; {ts['assessable_families']:,} families with eligible direct human single-protein annotations | Biological relevance, with explicit missing coverage |",
             f"| Drug signatures | {ss['paired_drug_signatures']:,} paired up/down contexts, {ss['unique_drugs']:,} compound names | Public responses across cells, doses and times |",
             f"| Same-cell proxy | {ss['knockout_proxy_signatures']} BUB1B knockout signatures in {len(ss['knockout_proxy_cells'])} cancer cell lines; {ss['drugs_with_same_cell_coverage']:,} names covered | Loss-of-function proxy, with different drug and knockout durations |",
             f"| Secondary proxy | BUB1B knockdown in KOLF2.1J human induced pluripotent stem cells | Compared with drug experiments in other cells; weaker context match |",
             f"| Joined review | {len(matrix):,} exact drug-name or source-synonym matches, {summary['distinct_mapped_parent_families']:,} distinct parents | Other ambiguous and unmatched names remain in the full screen |", "",
             "Sources: [ChEMBL data services](https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services), [Reactome downloads](https://reactome.org/download-data), [Perturb-Seqr bulk signatures](https://perturbseqr.maayanlab.cloud/download). Reference releases and exact bytes are recorded in the manifests. Public bulk downloads are matched locally; patient files and mutation coordinates were not sent to these services.", "",
             "## How this answers the mechanism question", "",
             "1. Start from the conditional BUB1B loss hypothesis and retain the uncertainty of the missense allele and phase. Public knockout/knockdown signatures approximate some consequences of reduced function; they are not a measured patient signature.",
             "2. Screen the entire downloaded drug universe against the selected biological modules and against opposing changes in public gene-expression lists. Retain inactive, missing, ambiguous and conflicting results.",
             "3. Review the direction of target action, cell-model consistency, current formulation/approval, achievable exposure and normal-versus-tumor effects. A drug touching the pathway can aggravate its failure.",
             "4. Reserve structural docking for candidates with a defensible protein target, binding pocket and structure. Docking estimates a possible binding pose, not restoration of function. The stop allele may lose RNA through NMD; docking a full-length protein cannot represent that loss. A missense structural model also needs independent structural and functional checks.", "",
             "The human [2WVI structure](https://www.rcsb.org/structure/2WVI) covers the BUBR1 N-terminal region, while [6JKM](https://www.rcsb.org/structure/6JKM) is a Drosophila kinase-domain structure. Neither should be presented as a complete human mutant structure. The [AlphaFold guidance](https://www.ebi.ac.uk/training/online/courses/alphafold/an-introductory-guide-to-its-strengths-and-limitations/strengths-and-limitations-of-alphafold/) describes limitations for mutation effects and ligand context. No docking run is included in these results.", "",
             "## Results and the limits they reveal", "",
             f"The strict target layer finds {ts['families_with_primary_overlap']} drug families with direct single-protein overlap and {ts['families_with_contextual_overlap']} with broader contextual overlap. Only {ts['module_genes_in_target_universe']} module genes are present in its eligible target universe: {', '.join(ts['tested_module_genes'])}. None passes the corrected threshold. This is sparse coverage, not evidence that the remaining drugs cannot help. Complex, protein-family and fusion components are retained as context without treating every component as a proven binding target.", "",
             "The target test compares overlap to a hypergeometric null and a separate exact null stratified by target annotation popularity. Both include all assessable parent families in Benjamini–Hochberg correction. The strata do not remove correlated drug-class bias, establish causal pathway direction, or measure protein-interaction network distance. Selected modules are co-membership sets, not a signed downstream network.", "",
             "The signature score compares signed membership in reported up/down gene lists. Positive values indicate opposition and negative values similarity. Unlisted genes are represented as zero for calculation; that does not mean they were measured and unchanged. The score is not a correlation of raw expression magnitudes. It has no calibrated p-value, success probability or therapeutic threshold. All doses and times contribute through condition and cell medians, but this can conceal heterogeneity. Public cancer cells and acute gene depletion differ from congenital disease.", "",
             "The LINCS source measures 978 landmark genes and infers many others; these files preserve pruned top-gene lists. Its Characteristic Direction analysis uses other profiles within a batch as controls, rather than necessarily matched untreated samples. Consequently the comparison describes relative list opposition. See the [source processing methods](https://maayanlab.cloud/sigcom-lincs/static/tutorial/ldp3/LDP3Tutorial.html). Several names pool multiple BRD source identifiers; an exact ChEMBL name match does not establish identical chemical structures across those records. Source identifiers must be resolved before choosing a molecule.", "",
             "### First exact mapped names in the descriptive screen order", "",
             "These are examples for auditing, not a selected treatment shortlist. Tiny positive secondary scores are not strong replication. Ordering favors at least two same-cell contexts before the primary score; current approval, exposure and benefit have not been established. Even the leading miconazole name has just three same-cell comparisons, including one negative score; positive cell medians do not mean every condition agrees.", "",
             "| Compound | Same-cell score | Matched cells | Stem-cell proxy across other cells | Review flag |", "| --- | ---: | ---: | ---: | --- |"]
    for row in matrix[:12]:
        flags = []
        if row["withdrawal_flag_present"] or row["withdrawal_warning_present"]:
            flags.append("source withdrawal flag")
        if row["direction_flags"]:
            flags.append("mitotic mechanism concern")
        if row["hipsc_cross_cell_score"] is not None and row["hipsc_cross_cell_score"] <= 0:
            flags.append("proxy directions disagree")
        lines.append(f"| {row['compound_name']} | {fmt(row['same_cell_score'])} | {row['same_cell_count']} | {fmt(row['hipsc_cross_cell_score'])} | {', '.join(flags) or 'benefit and exposure unresolved'} |")
    colchicine = prior["colchicine"]
    if colchicine:
        lines += ["", f"**The mechanistic warning is colchicine:** its same-cell reversal score is {fmt(colchicine['same_cell_score'])}, despite a microtubule mechanism that conflicts with the restoration rationale. This known concerning comparator is a false positive if reversal is interpreted as benefit. It demonstrates that the current score cannot select therapies. [The spindle-checkpoint biology](https://reactome.org/content/detail/R-HSA-69618) and the source drug-target annotations must constrain interpretation."]
    lines += ["", "Absence of a direction or withdrawal flag is not evidence of a safe or favorable mechanism. The limited direction rules cover selected checkpoint components and microtubules; every retained lead needs an individual review.", "",
              "### What happened to the earlier probes", "",
              "| Prior probe | Same-cell score | Matched cells | Secondary score |", "| --- | ---: | ---: | ---: |"]
    for name in ("azithromycin", "gentamicin", "amikacin"):
        row = prior[name]
        lines.append(f"| {name} | {fmt(row['same_cell_score']) if row else 'no exact named context'} | {row['same_cell_cell_count'] if row else 'unavailable'} | {fmt(row['hipsc_cross_cell_score']) if row else 'unavailable'} |")
    lines += ["", "Missing gentamicin under an exact name is a coverage gap. Signature similarity does not specifically assay readthrough, so these scores neither confirm nor disprove the original mechanistic hypotheses.", "",
              "## The next discriminating tests", "",
              "- Calibrate against randomly selected or matched control gene perturbations and drug classes, then correct the full set of candidate tests. The present signature ranks are descriptive.",
              "- Test sensitivity to leaving out a cell, dose, time and pathway module. Separate proliferation arrest, stress, cell death and compensatory responses from chromosome-fidelity recovery.",
              "- Inspect raw expression effect sizes, replicate agreement and quality measures for the surviving contexts. Thresholded gene lists lose information and can contain inferred genes.",
              "- For any structural hypothesis, establish that the expressed protein and relevant domain exist, compare reference and altered models, and validate docking with known ligands/decoys before interpreting binding. Binding alone does not establish stabilization or rescue.",
              "- Prioritize experiments that restore chromosome segregation per completed division, with protein identity, viability, exposure and separate malignant-cell controls.", "",
              "## Review files", "",
              "- [Combined evidence matrix](evidence_matrix.csv)",
              "- [All signature rankings](signatures/drug_rankings.csv)",
              "- [Signature methods and source details](signatures/SIGNATURE_SCREEN.md)",
              "- [Drug target and pathway rankings](targets/ranked_drugs.csv)",
              "- [Selected and excluded pathway modules](targets/module_selection.json)",
              "- [Exact run fingerprints](run_manifest.json)", "",
              "The original strategies remain available and unchanged. This screen adds a broad, reproducible evidence layer; it has not established a new effective drug or a higher contest score.", "",
              "## Acknowledgement", "", ACKNOWLEDGEMENT, ""]
    (root / "REVIEW.md").write_text("\n".join(lines), encoding="utf-8")
    original = baseline_fingerprints(Path("outputs/review_ready"))
    old = read(Path("outputs/strategy_b/run_manifest.json"))
    if any(sha256(Path("outputs/strategy_b") / name) != expected for name, expected in old["outputs_sha256"].items()):
        raise ValueError("Prior Strategy B output changed")
    if original != old["baseline_sha256"]:
        raise ValueError("Prior baseline differs from recorded hashes")
    inputs = [catalog_path, root / "targets/run_manifest.json", root / "signatures/summary.json"]
    files = [p for p in sorted(root.rglob("*")) if p.is_file() and p != root / "run_manifest.json"]
    write(root / "run_manifest.json", {"status": "completed", "completed_at_utc": datetime.now(timezone.utc).isoformat(),
          "original_baseline_unchanged": True, "strategy_b_outputs_unchanged": True,
          "input_sha256": {str(p.resolve()): sha256(p) for p in inputs}, "review_code_sha256": sha256(Path(__file__)),
          "outputs_sha256": {p.relative_to(root).as_posix(): sha256(p) for p in files}, "summary": summary})
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("outputs/drug_screen_c"))
    parser.add_argument("--catalog", type=Path, default=Path("references/drug_screen/catalog.json"))
    args = parser.parse_args()
    print(json.dumps(build(args.root, args.catalog), indent=2))
