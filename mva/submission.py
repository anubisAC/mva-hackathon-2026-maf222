"""Validate official Track 1 CSV shape locally without submitting anything."""

from __future__ import annotations

import csv
import math
from pathlib import Path

FIELDS = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"]


def validate(path: Path) -> dict:
    errors = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = set(FIELDS[:-1]) - set(reader.fieldnames or [])
        if missing:
            return {"valid": False, "errors": ["Missing columns: " + ", ".join(sorted(missing))], "rows": 0}
        rows = list(reader)
    if not 1 <= len(rows) <= 10:
        errors.append("The challenge requires 1 to 10 candidate rows")
    previous = 1.0
    seen = set()
    for index, row in enumerate(rows, 1):
        if row["proband_id"] != "PROBAND01":
            errors.append(f"Row {index}: expected challenge identifier PROBAND01")
        if row["finding_type"] not in {"primary", "secondary"}:
            errors.append(f"Row {index}: finding_type must be primary or secondary")
        try:
            epcr = float(row["epcr"])
            if not math.isfinite(epcr) or not 0 < epcr <= 1:
                raise ValueError
            if epcr > previous:
                errors.append(f"Row {index}: EPCR must be in descending order")
            previous = epcr
        except (ValueError, TypeError):
            errors.append(f"Row {index}: supply a reviewed EPCR in (0, 1]; heuristic scores are not probabilities")
        variants = []
        for suffix in ("1", "2"):
            values = [row.get(f"{key}_{suffix}", "") or "" for key in ("chrom", "pos", "ref", "alt")]
            if suffix == "2" and not any(values):
                continue
            if not all(values):
                errors.append(f"Row {index}: incomplete variant {suffix}")
                continue
            chrom, pos, ref, alt = values
            if chrom not in {f"chr{n}" for n in range(1, 23)} | {"chrX", "chrY", "chrM"}:
                errors.append(f"Row {index}: use a canonical chr-prefixed chromosome")
            if not pos.isdigit() or int(pos) < 1:
                errors.append(f"Row {index}: position must be a positive integer")
            if not set(ref) <= set("ACGT") or not set(alt) <= set("ACGT") or ref == alt:
                errors.append(f"Row {index}: require distinct concrete uppercase DNA alleles")
            variants.append(tuple(values))
        if len(variants) == 2 and variants[0] == variants[1]:
            errors.append(f"Row {index}: compound pair repeats the same variant")
        key = tuple(sorted(variants))
        if key in seen:
            errors.append(f"Row {index}: duplicate hypothesis")
        seen.add(key)
    return {"valid": not errors, "errors": errors, "rows": len(rows), "scope": "Format check only; does not verify GRCh38 mapping, causality, phase, probability calibration or competition score."}
