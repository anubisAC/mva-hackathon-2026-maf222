"""Local input inventory and clinical-feature extraction."""

from __future__ import annotations

import gzip
import hashlib
import json
import re
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_phenotype(path: Path) -> dict:
    with ZipFile(path) as document:
        root = ElementTree.fromstring(document.read("word/document.xml"))
    namespace = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    features = []
    for row in root.findall(".//w:tr", namespace):
        cells = [" ".join(node.itertext()).strip() for node in row.findall("w:tc", namespace)]
        if len(cells) >= 4 and re.fullmatch(r"HP:\d{7}", cells[2]):
            features.append({"clinical_feature": cells[0], "label": cells[1], "hpo_id": cells[2], "notes": cells[3], "subject": "family_history" if "parent" in cells[0].lower() else "proband"})
    if not features:
        raise ValueError("No HPO table found in the supplied clinical document")
    return {"source_filename": path.name, "sha256": sha256(path), "features": features, "interpretation": "Family history is retained separately; no absent features were inferred. HPO features are recorded for review, not yet used for an automated similarity score."}


def inspect_data(data_dir: Path, output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    files = [{"name": p.name, "bytes": p.stat().st_size} for p in sorted(data_dir.iterdir()) if p.is_file() and not p.name.startswith(".")]
    vcfs = list(data_dir.glob("*.vcf.gz")) + list(data_dir.glob("*.vcf"))
    docs = list(data_dir.glob("*.docx"))
    if len(vcfs) != 1 or len(docs) != 1:
        raise ValueError("Expected exactly one VCF and one clinical DOCX in the data directory")
    fastqs = sorted(data_dir.glob("*.fastq.gz"))
    pairs = []
    for r1 in fastqs:
        if "_R1_" not in r1.name:
            continue
        r2 = r1.with_name(r1.name.replace("_R1_", "_R2_"))
        pair = {"read1": r1.name, "read2": r2.name, "mate_exists": r2.exists(), "sampled_pairs": 0, "valid_structure": r2.exists(), "paired_identifiers_match": r2.exists()}
        if r2.exists():
            with gzip.open(r1, "rt", encoding="ascii") as first, gzip.open(r2, "rt", encoding="ascii") as second:
                for _ in range(100):
                    one, two = [first.readline().rstrip() for _ in range(4)], [second.readline().rstrip() for _ in range(4)]
                    if not one[0] and not two[0]:
                        break
                    pair["sampled_pairs"] += 1
                    pair["valid_structure"] &= all(x[0].startswith("@") and x[2].startswith("+") and len(x[1]) == len(x[3]) and len(x[1]) > 0 for x in (one, two))
                    if not one[0] or not two[0]:
                        pair["paired_identifiers_match"] = False
                        break
                    pair["paired_identifiers_match"] &= one[0].split()[0].removesuffix("/1") == two[0].split()[0].removesuffix("/2")
        pairs.append(pair)
    phenotype = read_phenotype(docs[0])
    (output_dir / "phenotype.json").write_text(json.dumps(phenotype, indent=2) + "\n", encoding="utf-8")
    result = {"files": files, "total_bytes": sum(f["bytes"] for f in files), "vcf": str(vcfs[0].resolve()), "vcf_sha256": sha256(vcfs[0]), "phenotype_count": len(phenotype["features"]), "fastq_pairs": pairs, "fastq_scope": "First 100 read pairs per lane only. Not full FASTQ QC, CRC validation, alignment, variant confirmation or CNV analysis."}
    (output_dir / "inventory.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result
