"""Public bulk reference acquisition and conversion for strategy B.

No patient data or patient-derived queries are sent to these services.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from .intake import sha256


PUBLIC_FILES = {
    "annotation.gtf.gz": "https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz",
    "cds.fa.gz": "https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz",
    "phenotype_to_genes.txt": "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/phenotype_to_genes.txt",
    "phenotype.hpoa": "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/phenotype.hpoa",
    "clinvar.vcf.gz": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
}


def download(url: str, destination: Path) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    receipt = destination.with_name(destination.name + ".source.json")
    if destination.exists() and receipt.exists():
        previous = json.loads(receipt.read_text(encoding="utf-8"))
        if previous.get("url") == url and previous.get("sha256") == sha256(destination):
            print(f"Using verified cache: {destination.name}", flush=True)
            return previous
    part = destination.with_name(destination.name + ".part")
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "MVA-local-research/0.2"})
            digest, total = hashlib.sha256(), 0
            print(f"Downloading public reference: {destination.name}", flush=True)
            with urllib.request.urlopen(request, timeout=120) as response, part.open("wb") as handle:
                resolved_url = response.geturl().split("?")[0]
                last_modified = response.headers.get("Last-Modified")
                expected_bytes = response.headers.get("Content-Length")
                for block in iter(lambda: response.read(4 * 1024 * 1024), b""):
                    handle.write(block)
                    digest.update(block)
                    total += len(block)
            if expected_bytes and total != int(expected_bytes):
                raise ValueError(f"Truncated download: {destination.name}")
            part.replace(destination)
            result = {"url": url, "resolved_url": resolved_url, "retrieved_at_utc": datetime.now(timezone.utc).isoformat(), "last_modified": last_modified, "bytes": total, "sha256": digest.hexdigest()}
            receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
            print(f"Downloaded {destination.name}: {total:,} bytes", flush=True)
            return result
        except (urllib.error.URLError, TimeoutError):
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def acquire(directory: Path) -> dict:
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {name: pool.submit(download, url, directory / name) for name, url in PUBLIC_FILES.items()}
        receipts, errors = {}, {}
        for name, future in futures.items():
            try:
                receipts[name] = future.result()
            except Exception as error:
                errors[name] = str(error)
    result = {"files": receipts, "errors": errors}
    (directory / "acquisition.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def attributes(raw: str) -> dict:
    result = {}
    for key, value in re.findall(r'(\w+) "([^"]*)"', raw):
        if key == "tag":
            result.setdefault("tag", []).append(value)
        else:
            result[key] = value
    return result


def build_gene_models(directory: Path) -> dict:
    """Select one canonical protein-coding transcript per gene by public tags.

    GTF CDS excludes the terminal stop, while Ensembl CDS fasta includes it.
    Include stop_codon intervals explicitly before validating lengths.
    """
    gtf_path = directory / "annotation.gtf.gz"
    transcripts = {}
    skipped = {}
    with gzip.open(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"transcript", "exon", "CDS", "stop_codon"}:
                continue
            info = attributes(fields[8])
            if info.get("gene_biotype") != "protein_coding" or "Ensembl_canonical" not in info.get("tag", []):
                continue
            if info.get("transcript_biotype") != "protein_coding":
                continue
            chrom, _, feature, start, end, _, strand, _, _ = fields
            if chrom not in {str(i) for i in range(1, 23)} | {"X", "Y"}:
                continue
            transcript_id = info.get("transcript_id")
            if not transcript_id:
                continue
            model = transcripts.setdefault(transcript_id, {"symbol": info.get("gene_name", info["gene_id"]), "gene_id": info["gene_id"], "chrom": chrom, "start": int(start), "end": int(end), "strand": 1 if strand == "+" else -1, "transcript_id": transcript_id, "transcript_version": info.get("transcript_version"), "exons": [], "coding_intervals": [], "transcript_selection": "Ensembl canonical protein coding transcript, release 116"})
            model["start"], model["end"] = min(model["start"], int(start)), max(model["end"], int(end))
            if feature == "exon":
                model["exons"].append({"start": int(start), "end": int(end)})
            elif feature in {"CDS", "stop_codon"}:
                model["coding_intervals"].append((int(start), int(end)))
    print(f"Canonical protein-coding transcript models found: {len(transcripts):,}", flush=True)
    current, sequence = None, []

    def save_sequence():
        if current in transcripts:
            transcripts[current]["cds_sequence"] = "".join(sequence).upper()

    with gzip.open(directory / "cds.fa.gz", "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                save_sequence()
                current, sequence = line[1:].split()[0].split(".")[0], []
            elif current in transcripts:
                sequence.append(line.strip())
        save_sequence()
    valid = []
    for model in transcripts.values():
        reason = None
        coding = model.pop("coding_intervals")
        sequence = model.get("cds_sequence", "")
        if not coding or not model["exons"] or not sequence:
            reason = "missing_sequence_or_coding_intervals"
        else:
            model["translation_start"] = min(start for start, _ in coding)
            model["translation_end"] = max(end for _, end in coding)
            expected = sum(max(0, min(e["end"], model["translation_end"]) - max(e["start"], model["translation_start"]) + 1) for e in model["exons"])
            if len(sequence) != expected:
                reason = "cds_length_mismatch"
            elif len(sequence) % 3 or sequence[:3] != "ATG" or sequence[-3:] not in {"TAA", "TAG", "TGA"}:
                reason = "incomplete_or_nonstandard_cds"
            elif not set(sequence) <= set("ACGT"):
                reason = "ambiguous_cds"
            elif any(sequence[i:i + 3] in {"TAA", "TAG", "TGA"} for i in range(0, len(sequence) - 3, 3)):
                reason = "internal_stop_or_recoded_transcript"
        if reason:
            skipped.setdefault(reason, []).append({"gene": model["symbol"], "transcript": model["transcript_id"]})
        else:
            valid.append(model)
    symbol_counts = {}
    for model in valid:
        symbol_counts[model["symbol"]] = symbol_counts.get(model["symbol"], 0) + 1
    unique = []
    for model in valid:
        if symbol_counts[model["symbol"]] != 1:
            skipped.setdefault("ambiguous_duplicate_symbol", []).append({"gene": model["symbol"], "transcript": model["transcript_id"]})
        else:
            unique.append(model)
    valid = sorted(unique, key=lambda m: (m["chrom"], m["start"], m["symbol"]))
    summary = {"assembly": "GRCh38", "release": 116, "candidate_models": len(transcripts), "usable_models": len(valid), "skipped_counts": {key: len(value) for key, value in skipped.items()}, "source_hashes": {name: sha256(directory / name) for name in ("annotation.gtf.gz", "cds.fa.gz")}, "limitations": ["Canonical nuclear coding transcripts only; mitochondrial, noncoding, incomplete/nonstandard CDS and alternative transcript effects are not covered.", "Selection rules apply to every gene without disease-specific selection."]}
    with gzip.open(directory / "genes_grch38.jsonl.gz", "wt", encoding="utf-8") as handle:
        handle.write(json.dumps(summary) + "\n")
        for model in valid:
            handle.write(json.dumps(model) + "\n")
    (directory / "gene_models_summary.json").write_text(json.dumps({**summary, "skipped_models": skipped}, indent=2) + "\n", encoding="utf-8")
    return summary


def build_phenotype_reference(hpo_path: Path, phenotype_path: Path, output_path: Path) -> dict:
    """Combine public HPO gene associations with locally read patient features."""
    genes, gene_diseases = {}, {}
    rows = 0
    with hpo_path.open(encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            # HPO release headers use ncbi_gene_id/gene_symbol/hpo_id.
            symbol = row.get("gene_symbol") or row.get("entrez-gene-symbol")
            hpo_id = row.get("hpo_id") or row.get("#hpo_id")
            if symbol and hpo_id and hpo_id.startswith("HP:"):
                genes.setdefault(symbol, set()).add(hpo_id)
                if row.get("disease_id"):
                    gene_diseases.setdefault(symbol, set()).add(row["disease_id"])
                rows += 1
    if not genes:
        raise ValueError("No gene/HPO associations parsed; inspect current HPO headers")
    hpoa_path = hpo_path.with_name("phenotype.hpoa")
    inheritance_modes = {"HP:0000006": "autosomal_dominant", "HP:0000007": "autosomal_recessive", "HP:0001417": "x_linked", "HP:0001419": "x_linked_recessive", "HP:0001423": "x_linked_dominant", "HP:0001427": "mitochondrial", "HP:0001450": "y_linked"}
    disease_modes = {}
    if hpoa_path.exists():
        with hpoa_path.open(encoding="utf-8-sig") as handle:
            reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
            for row in reader:
                identifier = row.get("database_id") or row.get("databaseid")
                mode = inheritance_modes.get(row.get("hpo_id", ""))
                if identifier and mode and row.get("qualifier", "") != "NOT":
                    disease_modes.setdefault(identifier, set()).add(mode)
    gene_inheritance, inheritance_sources = {}, {}
    for gene, diseases in gene_diseases.items():
        evidence = {disease: sorted(disease_modes[disease]) for disease in diseases if disease in disease_modes}
        if evidence:
            gene_inheritance[gene] = sorted({mode for modes in evidence.values() for mode in modes})
            inheritance_sources[gene] = evidence
    phenotype = json.loads(phenotype_path.read_text(encoding="utf-8"))
    patient, family = [], []
    for feature in phenotype["features"]:
        item = {"hpo_id": feature["hpo_id"], "label": feature["label"], "source": feature["subject"]}
        (family if feature["subject"] == "family_history" else patient).append(item)
    result = {"patient_features": patient, "family_features": family, "genes": {key: sorted(value) for key, value in genes.items()}, "gene_inheritance": gene_inheritance, "gene_inheritance_sources": inheritance_sources, "source": {"hpo_path": str(hpo_path.resolve()), "hpo_sha256": sha256(hpo_path), "hpo_release": "v2026-09-01", "inheritance_path": str(hpoa_path.resolve()) if hpoa_path.exists() else None, "inheritance_sha256": sha256(hpoa_path) if hpoa_path.exists() else None, "clinical_sha256": phenotype.get("sha256"), "association_rows": rows, "note": "HPO phenotype_to_genes release includes its own annotation propagation; ranking uses exact supplied term membership, not a bespoke semantic similarity algorithm. Inheritance is a gene-level union of associated disease annotations, not confirmation of the causal condition or mode."}}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result) + "\n", encoding="utf-8")
    return {"gene_count": len(genes), "association_rows": rows, "inheritance_genes": len(gene_inheritance), "patient_features": len(patient), "family_features": len(family), "output": str(output_path)}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["download", "models", "phenotypes"])
    parser.add_argument("--directory", type=Path, default=Path("references/alternative"))
    args = parser.parse_args()
    if args.action == "download":
        print(json.dumps(acquire(args.directory), indent=2))
    elif args.action == "models":
        print(json.dumps(build_gene_models(args.directory), indent=2))
    else:
        print(json.dumps(build_phenotype_reference(args.directory / "phenotype_to_genes.txt", Path("outputs/review_ready/intake/phenotype.json"), Path("outputs/strategy_b_inputs/phenotype_reference.json")), indent=2))
