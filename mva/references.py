"""Download public reference information, never patient queries or coordinates."""

from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


def fetch(url: str, destination: Path) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    receipt = destination.with_suffix(destination.suffix + ".source.json")
    if destination.exists() and receipt.exists():
        prior = json.loads(receipt.read_text(encoding="utf-8"))
        if prior.get("url") == url and prior.get("sha256") == hashlib.sha256(destination.read_bytes()).hexdigest():
            return prior
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "MVA-local-research/0.1", "Accept": "application/json"})
            with urllib.request.urlopen(request, timeout=90) as response:
                content = response.read()
            break
        except (urllib.error.URLError, TimeoutError):
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
    destination.write_bytes(content)
    metadata = {"url": url, "retrieved_at": datetime.now(timezone.utc).isoformat(), "sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)}
    receipt.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata


def prepare_genes(config_path: Path, output_dir: Path) -> dict:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config["assembly"] != "GRCh38":
        raise ValueError("Current Ensembl reference preparation supports GRCh38 only")
    output_dir.mkdir(parents=True, exist_ok=True)
    result = {"assembly": "GRCh38", "genes": [], "sources": [], "scope": config["description"]}
    release_path = output_dir / "ensembl_data_release.json"
    result["sources"].append(fetch("https://rest.ensembl.org/info/data?content-type=application/json", release_path))
    result["ensembl_release"] = json.loads(release_path.read_text(encoding="utf-8"))
    for entry in config["genes"]:
        symbol = entry["symbol"]
        print(f"Preparing public reference for {symbol}", flush=True)
        lookup_path = output_dir / "raw" / f"{symbol}.json"
        url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/" + urllib.parse.quote(symbol) + "?expand=1;content-type=application/json"
        result["sources"].append(fetch(url, lookup_path))
        gene = json.loads(lookup_path.read_text(encoding="utf-8"))
        if gene.get("assembly_name") != "GRCh38":
            raise ValueError(f"Unexpected assembly for {symbol}")
        canonical = gene.get("canonical_transcript", "").split(".")[0]
        transcripts = [t for t in gene["Transcript"] if t["id"] == canonical and "Translation" in t]
        if len(transcripts) != 1:
            raise ValueError(f"Expected one protein-coding canonical transcript for {symbol}")
        transcript = transcripts[0]
        cds_path = output_dir / "raw" / f"{transcript['id']}.cds.json"
        url = f"https://rest.ensembl.org/sequence/id/{transcript['id']}?type=cds;content-type=application/json"
        result["sources"].append(fetch(url, cds_path))
        cds = json.loads(cds_path.read_text(encoding="utf-8"))["seq"].upper()
        translation = transcript["Translation"]
        start, end = sorted([translation["start"], translation["end"]])
        exons = [{"start": e["start"], "end": e["end"]} for e in transcript["Exon"]]
        expected = sum(max(0, min(e["end"], end) - max(e["start"], start) + 1) for e in exons)
        if len(cds) != expected:
            raise ValueError(f"CDS mapping mismatch for {symbol}: sequence {len(cds)}, intervals {expected}")
        result["genes"].append({**entry, "chrom": str(gene["seq_region_name"]), "start": gene["start"], "end": gene["end"], "strand": transcript["strand"], "transcript_id": transcript["id"], "transcript_version": transcript.get("version"), "transcript_selection": "Ensembl canonical, not independently MANE-validated", "exons": exons, "translation_start": start, "translation_end": end, "cds_sequence": cds})
    result["created_at"] = datetime.now(timezone.utc).isoformat()
    result_path = output_dir / "genes_grch38.json"
    result_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return {"path": str(result_path), "genes": len(result["genes"]), "ensembl_release": result["ensembl_release"]}


def prepare_enrichment(output_dir: Path, aliases_path: Path | None = None) -> dict:
    """Cache public human TF libraries locally; no gene lists are uploaded."""
    output_dir.mkdir(parents=True, exist_ok=True)
    aliases_path = aliases_path or Path("config/gene_aliases.json")
    alias_records = json.loads(aliases_path.read_text(encoding="utf-8"))["aliases"] if aliases_path.exists() else []
    aliases = {row["legacy"]: row["current"] for row in alias_records}
    libraries = []
    for library in ("ChEA_2016", "TRANSFAC_and_JASPAR_PWMs"):
        raw = output_dir / f"{library}.raw.gmt"
        source = fetch(f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={library}", raw)
        lines = raw.read_text(encoding="utf-8").splitlines()
        if library == "ChEA_2016":
            lines = [line for line in lines if "HUMAN" in line.split("\t")[0].upper()]
        normalized, universe = [], set()
        for line in lines:
            columns = line.split("\t")
            if len(columns) < 3:
                continue
            raw_genes = [g.split(",")[0].strip() for g in columns[2:] if g.strip()]
            genes = sorted({aliases.get(g, g) for g in raw_genes})
            universe.update(genes)
            normalized.append("\t".join([columns[0], columns[1], *genes]))
        if not normalized:
            raise ValueError(f"Empty gene set library: {library}")
        gmt = output_dir / f"{library}.gmt"
        background = output_dir / f"{library}.universe.txt"
        gmt.write_text("\n".join(normalized) + "\n", encoding="utf-8")
        background.write_text("\n".join(sorted(universe)) + "\n", encoding="utf-8")
        libraries.append({"name": library, "gmt": str(gmt), "universe": str(background), "sets": len(normalized), "genes": len(universe), "source": source, "background_caveat": "Union of genes annotated in this library, not a measured-expression background. Exploratory literature-set analysis only.", "context": "Human-only ChEA experiments" if library == "ChEA_2016" else "Motif-predicted targets; library mixes TF contexts and motifs, not patient binding evidence."})
    metadata = {"libraries": libraries, "symbol_aliases_applied": alias_records, "multiple_testing": "BH correction within each library; exploratory cross-library comparisons are not independent validation."}
    (output_dir / "libraries.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata
