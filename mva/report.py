"""Connect outputs from both tracks into a local, reviewable research record."""

from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path

from .intake import sha256
from .submission import FIELDS


ACKNOWLEDGEMENT = (
    "This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible."
)


def read(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def build(output_dir: Path, config_path: Path) -> dict:
    qc = read(output_dir / "track1/qc.json")
    hypotheses = read(output_dir / "track1/hypotheses.json")
    candidates = read(output_dir / "track1/candidates.json")
    drugs = read(output_dir / "track2/track2_candidates.json")
    snapshot = output_dir / "track2/evidence_snapshot.json"
    evidence_path = snapshot if snapshot.exists() else config_path
    if sha256(evidence_path) != drugs["config_sha256"]:
        raise ValueError("Track 2 evidence differs from this run. Restore its original ledger or create a new run; historical evidence cannot be silently replaced.")
    config = read(evidence_path)
    if not snapshot.exists():
        snapshot.write_bytes(evidence_path.read_bytes())
    by_id = {c["variant_id"]: c for c in candidates}
    links = []
    for hypothesis in hypotheses:
        mechanisms = [m for m in config["mechanisms"] if m["gene"] == hypothesis["gene"]]
        links.append({"gene": hypothesis["gene"], "variant_hypothesis": hypothesis, "mechanism_branches": mechanisms, "status": "conditional_unvalidated", "interpretation": "A literature mechanism can motivate experiments but does not establish this person's variant function or a drug response."})
    (output_dir / "mechanism_bridge.json").write_text(json.dumps(links, indent=2) + "\n", encoding="utf-8")
    staging = output_dir / "submission_review"
    staging.mkdir(exist_ok=True)
    # EPCR is intentionally empty until the human scientific review: a rank score
    # must never be presented as a calibrated causal probability.
    draft = staging / "track1_UNREVIEWED.csv"
    with draft.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for hypothesis in hypotheses[:10] if qc["scan_complete"] else []:
            row = dict.fromkeys(FIELDS, "")
            row.update(proband_id="PROBAND01", finding_type="primary", notes=f"UNREVIEWED {hypothesis['gene']}; {hypothesis['type']}; phase={hypothesis['phase']}; population frequencies={hypothesis['population_frequency_status']}; not a diagnostic conclusion")
            for i, variant_id in enumerate(hypothesis["variant_ids"], 1):
                variant = by_id[variant_id]
                for key in ("chrom", "pos", "ref", "alt"):
                    row[f"{key}_{i}"] = "chr" + variant[key] if key == "chrom" else variant[key]
            writer.writerow(row)
    enrichment = []
    for path in sorted((output_dir / "enrichment").glob("*/enrichment_results.json")):
        result = read(path)
        enrichment.append({"library": path.parent.name, "sets": len(result["results"]), "significant": result["significant_set_count"], "signature": result["signature"], "top": result["results"][:3]})
    counts = drugs["counts"]
    lines = ["# MVA first analysis review", "",
             "Both research tracks have run locally. Track 1 prioritizes molecular hypotheses from the supplied VCF; Track 2 evaluates public drug-mechanism evidence. These outputs support the next scientific review and are not a submitted entry or a validated diagnosis.", "",
             "## Track 1 results", "",
             f"The scan processed **{qc['records_seen']:,} VCF records**, with **{qc['candidate_count']:,} called alternate alleles** overlapping the {qc['panel_gene_count']}-gene starter panel and **{qc['hypothesis_count']} recessive hypotheses**. Scan complete: {qc['scan_complete']}. Assembly: {qc['assembly']}.", ""]
    for i, hypothesis in enumerate(hypotheses[:5], 1):
        variants = [by_id[v] for v in hypothesis["variant_ids"]]
        detail = "; ".join(f"{v['variant_id']} ({v['protein_change'] or v['consequence']}, {v['quality_status']} QC)" for v in variants)
        lines.append(f"{i}. **{hypothesis['gene']}** — {hypothesis['type'].replace('_', ' ')}: {detail}. Phase: {hypothesis['phase']}. Population frequency: {hypothesis['population_frequency_status']}.")
    lines += ["", "This ranks a literature-selected panel, not all genes. Missing population frequency is unknown, not evidence of rarity. Neither a high score nor two heterozygous alleles proves pathogenicity or trans phase. The HPO profile is preserved, including separate family history, but no phenotype similarity model has been fitted.", "",
              "[Detailed variant review](track1/review.md) · [Allele table](track1/candidates.csv) · [Variant hypotheses](track1/hypotheses.json)", "",
              "## Track 2 results", "",
              f"The evidence screen produced **{counts['accepted']} accepted experimental priorities**, **{counts['deferred']} deferred hypotheses or comparators**, and **{counts['rejected']} rejected candidate(s)**. Accepted would mean ready for further experimental prioritization, never established treatment efficacy.", "",
              *[f"- **{c['name']}** — {c['state']}. {' '.join(c['decision_reasons'])}" for c in drugs["candidates"]], "",
              "The molecular bridge links each Track 1 hypothesis to its literature mechanism without asserting that the mechanism has been established in this case.", "",
              *([f"The leading gene hypothesis is **{links[0]['gene']}**. Its conditional mechanism branch: " + " ".join(m["description"] for m in links[0]["mechanism_branches"])] if links and links[0]["mechanism_branches"] else ["A curated mechanism branch for the leading hypothesis has not yet been connected."]), "",
              "[Drug evidence and experiments](track2/track2_research_draft.md) · [Drug table](track2/track2_candidates.csv) · [Mechanism bridge](mechanism_bridge.json)", "",
              "## What carried over from DJ1", "",
              "The workflow follows the DJ1 paper's evidence collection, gene-set enrichment, regulator interpretation and experimental validation sequence. It retains the distinction between specific experimental findings and broad signatures. Enrichment uses explicit backgrounds, exact documented alias mapping and Benjamini–Hochberg correction over every tested set.", ""]
    for entry in enrichment:
        lines.append(f"- **{entry['library']}**: {entry['sets']} sets tested; {entry['significant']} at FDR < 0.05.")
        lines.append(f"  Signature: {entry['signature']['name']}; evidence type: {entry['signature']['evidence_type']}; direction: {entry['signature']['direction']}.")
        if entry["top"]:
            first = entry["top"][0]
            lines.append(f"  Strongest exploratory overlap: {first['set_name']}; raw p={first['p_value']:.4g}, FDR={first['fdr']:.4g}.")
    lines += ["", "Literature signatures do not measure patient expression. ChEA is restricted to human experiments; the motif library contains predicted targets across differing contexts. The library's covered-gene union is an exploratory background. Cross-library results are not independent validation. Gene-set overlap cannot establish up/down regulation or nominate a drug by expression reversal.", "",
              "## Next changes with the highest value", "",
              "1. Add local population-frequency and ClinVar annotations and independently confirm transcript consequences. Verify phase and inspect supporting reads before assigning EPCR.",
              "2. Extend beyond the starter panel with genome-wide functional annotation, phenotype/inheritance scoring and structural-variant/CNV assessment. Compare ranks with and without each evidence channel.",
              "3. Find relevant measured perturbation signatures or generate allele-specific expression/protein evidence. Compare checkpoint restoration, micronuclei, chromosome segregation and viability together to avoid mistaking cell death for rescue.",
              "4. Expand the approved-drug search around the supported mechanism, with matched controls, negative evidence, exposure feasibility and separate tumor/normal-cell effects. Do not force a positive candidate from this initial set.",
              "5. Turn the reviewed evidence into the official methods report and three-minute pitch; review the repository before publishing.", "",
              "## Submission status", "",
              "No upload or competition submission has occurred. The CSV in submission_review is deliberately marked UNREVIEWED and leaves EPCR blank. The local validator will reject it until those estimates are reviewed. GitHub publication, the final report, account-specific AI plan/data-setting disclosure and the pitch video remain outstanding.", "",
              "## Research provenance", "",
              "The analysis is hypothesis-led and not blinded. During requirements research, a public search surfaced other teams' reported results. No answer coordinates were seeded into the code or panel configuration. Performance on this case cannot be claimed as independent blinded recovery. Reference downloads are cached with URLs, retrieval times and SHA-256 checksums; run_manifest.json records code/configuration and input fingerprints.", "",
              "[Official requirements](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026)", "",
              "## Acknowledgement", "", ACKNOWLEDGEMENT, ""]
    review = output_dir / "REVIEW.md"
    review.write_text("\n".join(lines), encoding="utf-8")
    result = {"review": str(review.resolve()), "draft_csv": str(draft.resolve()), "scan_complete": qc["scan_complete"], "records": qc["records_seen"], "hypotheses": len(hypotheses), "drug_states": counts, "enrichment": [{k: v for k, v in e.items() if k not in {"top", "signature"}} for e in enrichment]}
    manifest_path = output_dir / "run_manifest.json"
    if manifest_path.exists():
        manifest = read(manifest_path)
        if manifest.get("status") == "completed":
            manifest.setdefault("report_revisions", []).append({"created_at_utc": datetime.now(timezone.utc).isoformat(), "code_sha256": sha256(Path(__file__)), "evidence_sha256": sha256(evidence_path)})
            manifest["summary"] = result
            manifest["outputs_sha256"] = {p.relative_to(output_dir).as_posix(): sha256(p) for p in sorted(output_dir.rglob("*")) if p.is_file() and p != manifest_path}
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return result
