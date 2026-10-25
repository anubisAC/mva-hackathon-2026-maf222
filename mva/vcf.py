"""Small streaming, single-sample VCF reader used by the local analysis.

Coordinates remain 1-based inclusive. INFO/AF is deliberately not interpreted as
a population frequency: in this dataset it describes the called sample.
"""
from __future__ import annotations

import gzip
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


GRCH38_LENGTHS = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
    "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
    "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
    "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
    "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
}
GRCH37_LENGTHS = {"1": 249250621, "2": 243199373, "X": 155270560}


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.removeprefix("chr")
    return "MT" if chrom in {"M", "MT"} else chrom


def number(value: str | None, integer: bool = False) -> int | float | None:
    if value in {None, "", "."}:
        return None
    try:
        result = int(value) if integer else float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def parse_info(raw: str) -> dict[str, str | bool]:
    result: dict[str, str | bool] = {}
    for item in raw.split(";"):
        if not item or item == ".":
            continue
        key, separator, value = item.partition("=")
        result[key] = value if separator else True
    return result


def parse_gt(raw: str, alt_count: int) -> tuple[tuple[int | None, ...], bool]:
    alleles: list[int | None] = []
    valid = True
    for token in re.split(r"[/|]", raw):
        if token in {"", "."}:
            alleles.append(None)
        elif token.isdigit() and 0 <= int(token) <= alt_count:
            alleles.append(int(token))
        else:
            alleles.append(None)
            valid = False
    return tuple(alleles), valid


@dataclass(slots=True)
class VCFRecord:
    chrom: str
    pos: int
    identifier: str
    ref: str
    alts: tuple[str, ...]
    qual: float | None
    filter_raw: str
    info_raw: str
    sample: dict[str, str]
    gt: tuple[int | None, ...]
    gt_valid: bool
    line_number: int

    @property
    def end(self) -> int:
        if any(alt.startswith("<") or "[" in alt or "]" in alt for alt in self.alts):
            end = number(str(parse_info(self.info_raw).get("END", ".")), integer=True)
            if end is not None and end >= self.pos:
                return int(end)
        return self.pos + len(self.ref) - 1

    @property
    def called_alt_indices(self) -> tuple[int, ...]:
        return tuple(sorted({a for a in self.gt if a is not None and a > 0}))


class VCFReader:
    """One-pass reader; header and audit state remain accessible after iteration."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.samples: list[str] = []
        self.contigs: dict[str, int] = {}
        self.header_lines: list[str] = []
        self.records_seen = 0
        self.malformed_records = 0
        self.malformed_examples: list[dict] = []
        self.complete = False
        self.assembly = "unknown"
        self.assembly_evidence: list[str] = []

    def _open(self) -> TextIO:
        with self.path.open("rb") as stream:
            compressed = stream.read(2) == b"\x1f\x8b"
        return gzip.open(self.path, "rt", encoding="utf-8") if compressed else self.path.open(encoding="utf-8")

    def _infer_assembly(self) -> None:
        matches38 = [c for c, size in self.contigs.items() if GRCH38_LENGTHS.get(c) == size]
        mismatches38 = [c for c, size in self.contigs.items() if c in GRCH38_LENGTHS and GRCH38_LENGTHS[c] != size]
        matches37 = [c for c, size in self.contigs.items() if GRCH37_LENGTHS.get(c) == size]
        explicit38 = any(re.search(r"(?:GRCh38|hg38)", line, re.I) for line in self.header_lines)
        explicit37 = any(re.search(r"(?:GRCh37|hg19)", line, re.I) for line in self.header_lines)
        if matches37 or (explicit37 and not explicit38):
            self.assembly = "GRCh37"
            self.assembly_evidence.append("GRCh37 canonical contig length or explicit header reference")
        elif mismatches38:
            self.assembly = "conflicting"
            self.assembly_evidence.append("Canonical contig lengths conflict with GRCh38: " + ", ".join(mismatches38))
        elif len(matches38) >= 2 or explicit38:
            self.assembly = "GRCh38"
            if matches38:
                self.assembly_evidence.append(f"{len(matches38)} canonical contig lengths match GRCh38")
            if explicit38:
                self.assembly_evidence.append("GRCh38 or hg38 named in the VCF header")
        else:
            self.assembly_evidence.append("Insufficient header evidence to verify the reference assembly")

    def records(self, max_records: int | None = None) -> Iterator[VCFRecord]:
        if max_records is not None and max_records < 1:
            raise ValueError("max_records must be positive")
        header_found = False
        with self._open() as stream:
            for line_number, line in enumerate(stream, 1):
                if line.startswith("##"):
                    self.header_lines.append(line.rstrip())
                    if line.startswith("##contig=<"):
                        identifier = re.search(r"(?:<|,)ID=([^,>]+)", line)
                        length = re.search(r",length=(\d+)", line, re.I)
                        if identifier and length:
                            self.contigs[normalize_chrom(identifier.group(1))] = int(length.group(1))
                    continue
                if line.startswith("#CHROM"):
                    fields = line.rstrip("\r\n").split("\t")
                    self.samples = fields[9:]
                    if len(self.samples) != 1:
                        raise ValueError(f"Expected exactly one VCF sample; found {len(self.samples)}")
                    self._infer_assembly()
                    header_found = True
                    continue
                if line.startswith("#") or not line.strip():
                    continue
                if not header_found:
                    raise ValueError("VCF has no #CHROM header before data")
                if max_records is not None and self.records_seen >= max_records:
                    return
                self.records_seen += 1
                try:
                    fields = line.rstrip("\r\n").split("\t")
                    if len(fields) != 10:
                        raise ValueError(f"Expected 10 single-sample columns, found {len(fields)}")
                    chrom, position, identifier, ref, alt, qual, filter_raw, info_raw, formats, sample_raw = fields
                    pos = int(position)
                    if pos < 1 or not ref or ref == ".":
                        raise ValueError("Invalid genomic position or reference allele")
                    alts = tuple(alt.upper().split(",")) if alt != "." else ()
                    sample = dict(zip(formats.split(":"), sample_raw.split(":")))
                    genotype, gt_valid = parse_gt(sample.get("GT", "."), len(alts))
                    record = VCFRecord(normalize_chrom(chrom), pos, identifier, ref.upper(), alts,
                                       number(qual), filter_raw, info_raw, sample, genotype, gt_valid, line_number)
                except (ValueError, IndexError) as error:
                    self.malformed_records += 1
                    if len(self.malformed_examples) < 20:
                        self.malformed_examples.append({"line": line_number, "error": str(error)})
                    continue
                yield record
        if not header_found:
            raise ValueError("VCF has no #CHROM header")
        self.complete = True
