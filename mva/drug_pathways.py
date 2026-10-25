"""Version-pinned public Reactome sets and local seed-module selection.

These are undirected co-membership hypotheses. They do not establish downstream
causality, drug effect direction, mutation rescue, or clinical benefit.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from .intake import sha256


RELEASE = 97
DEFAULT_CACHE = Path("references/drug_screen")
FILES = ("ReactomePathways.gmt.zip", "ReactomePathways.txt",
         "ReactomePathwaysRelation.txt", "pathway2summation.txt")
SOURCE_BASE = f"https://download.reactome.org/{RELEASE}/"
# Fixed before screening drug rankings. These labels exclude broad categories
# even if a later database release happens to shrink their gene counts.
BROAD_NAMES = frozenset({
    "cell cycle", "cell cycle, mitotic", "disease", "metabolism",
    "signal transduction", "immune system", "developmental biology",
    "programmed cell death", "dna repair", "metabolism of proteins",
    "gene expression (transcription)", "cellular responses to stimuli",
    "hemostasis",
})
DISEASE_ROOT = "R-HSA-1643685"
POLICY = {
    "min_genes": 5, "max_genes": 200,
    "selection": "leaf-most eligible seed-containing pathway for each seed",
    "exclude_top_level": True,
    "exclude_disease_root_and_descendants": DISEASE_ROOT,
    "broad_category_names": sorted(BROAD_NAMES),
    "size_definition": "unique GMT member symbols, before drug-target intersection",
    "deduplicate_gene_symbols": True,
    "semantics": "undirected pathway co-membership; no effect direction or benefit",
    "frozen_before_first_drug_screen": True,
}


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def download_file(cache_dir: Path, name: str) -> dict:
    """Fetch an entire public file; no gene or patient-dependent URL exists."""
    if name not in FILES:
        raise ValueError("Unknown bulk Reactome reference")
    directory = cache_dir / f"reactome_v{RELEASE}"
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / name
    receipt_path = destination.with_name(destination.name + ".receipt.json")
    url = SOURCE_BASE + name
    if destination.exists() and receipt_path.exists():
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if receipt["url"] == url and receipt["sha256"] == sha256(destination):
            return receipt
        raise ValueError(f"Cached reference fingerprint mismatch: {destination}")
    request = urllib.request.Request(url, headers={"User-Agent": "MVA-research-public-reference/1.0"})
    temporary = destination.with_name(destination.name + ".part")
    with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as handle:
        while block := response.read(1024 * 1024):
            handle.write(block)
        handle.flush()
        receipt = {
            "url": url, "resolved_url": response.url, "reactome_release": RELEASE,
            "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            "last_modified": response.headers.get("Last-Modified"),
            "etag": response.headers.get("ETag"),
            "bytes": temporary.stat().st_size, "sha256": sha256(temporary),
            "local_path": str(destination.resolve()),
            "request_scope": "complete public bulk file; no patient data or seed query",
        }
    temporary.replace(destination)
    write_json(receipt_path, receipt)
    return receipt


def parse_names(text: str) -> dict[str, dict]:
    pathways = {}
    for number, row in enumerate(csv.reader(io.StringIO(text), delimiter="\t"), 1):
        if not row:
            continue
        if len(row) != 3:
            raise ValueError(f"Pathway names line {number}: expected 3 columns")
        identifier, name, species = (value.strip() for value in row)
        if species != "Homo sapiens":
            continue
        if not identifier.startswith("R-HSA-"):
            raise ValueError("Human pathway has non-human identifier")
        if identifier in pathways:
            raise ValueError(f"Duplicate pathway identifier {identifier}")
        pathways[identifier] = {"id": identifier, "name": name, "species": species,
                                "gene_symbols": [], "parent_ids": [], "child_ids": [],
                                "description": "", "summations": []}
    return pathways


def normalize_reference(gmt_text: str, names_text: str, relations_text: str,
                        descriptions_text: str, receipts: list[dict] | None = None) -> dict:
    """Strictly parse official bulk formats, validate a DAG, retain full human sets."""
    pathways = parse_names(names_text)
    seen_gmt = set()
    for number, row in enumerate(csv.reader(io.StringIO(gmt_text), delimiter="\t"), 1):
        if not row:
            continue
        if len(row) < 3:
            raise ValueError(f"GMT line {number}: expected name, identifier and members")
        name, identifier, *genes = row
        name, identifier = name.strip(), identifier.strip()
        if not identifier.startswith("R-HSA-"):
            continue
        if identifier not in pathways:
            raise ValueError(f"GMT pathway absent from pathway names: {identifier}")
        if identifier in seen_gmt:
            raise ValueError(f"Duplicate GMT pathway {identifier}")
        canonical_name = pathways[identifier]["name"]
        # GMT disambiguates duplicate display names with the Reactome numeric ID.
        if name not in {canonical_name, canonical_name + "_" + identifier.rsplit("-", 1)[-1]}:
            raise ValueError(f"Pathway name mismatch for {identifier}")
        seen_gmt.add(identifier)
        pathways[identifier]["gmt_name"] = name
        pathways[identifier]["gene_symbols"] = sorted(set(g.strip() for g in genes if g.strip()))
    for number, row in enumerate(csv.reader(io.StringIO(relations_text), delimiter="\t"), 1):
        if not row:
            continue
        if len(row) != 2:
            raise ValueError(f"Hierarchy line {number}: expected parent and child")
        parent, child = row
        if not parent.startswith("R-HSA-") and not child.startswith("R-HSA-"):
            continue
        if parent not in pathways or child not in pathways:
            raise ValueError(f"Human hierarchy references unknown pathway: {parent}, {child}")
        pathways[parent]["child_ids"].append(child)
        pathways[child]["parent_ids"].append(parent)
    description_lines = descriptions_text.splitlines()
    if not description_lines or description_lines[0].split("\t") != ["Identifier", "Name", "Summation"]:
        raise ValueError("Unexpected Reactome summation header")
    for line in description_lines[1:]:
        # Summations are unquoted free text and can contain embedded tabs.
        columns = line.split("\t", 2)
        if len(columns) != 3:
            raise ValueError("Malformed Reactome summation row")
        identifier, _, description = columns
        if identifier not in pathways:
            continue
        if description not in pathways[identifier]["summations"]:
            pathways[identifier]["summations"].append(description)
    ancestors = {}

    def visit(identifier: str, visiting: set[str]) -> set[str]:
        if identifier in visiting:
            raise ValueError("Cycle in Reactome pathway hierarchy")
        if identifier in ancestors:
            return ancestors[identifier]
        result = set(pathways[identifier]["parent_ids"])
        for parent in pathways[identifier]["parent_ids"]:
            result.update(visit(parent, visiting | {identifier}))
        ancestors[identifier] = result
        return result

    for identifier, row in pathways.items():
        row["parent_ids"] = sorted(set(row["parent_ids"]))
        row["child_ids"] = sorted(set(row["child_ids"]))
        row["ancestor_ids"] = sorted(visit(identifier, set()))
        row["gene_count"] = len(row["gene_symbols"])
        row["gmt_available"] = identifier in seen_gmt
        row["description"] = "\n\n".join(row["summations"])
        row["source_url"] = f"https://reactome.org/content/detail/{identifier}"
    genes = sorted({gene for row in pathways.values() for gene in row["gene_symbols"]})
    return {
        "schema_version": 1, "database": "Reactome", "release": RELEASE,
        "pathway_species": "Homo sapiens", "source_receipts": receipts or [],
        "source_documentation": "https://reactome.org/download-data",
        "member_identifier_type": "Reactome GMT gene symbols, exact spelling; aliases not expanded",
        "member_species_caveat": "Human-host pathways can include pathogen symbols; intersect with human drug-target genes before testing.",
        "relationship_semantics": "Gene membership and parent-child pathway hierarchy; no signed causal drug-rescue relationships.",
        "pathway_count": len(pathways), "gmt_pathway_count": len(seen_gmt),
        "gene_universe": genes, "gene_universe_count": len(genes),
        "pathways": sorted(pathways.values(), key=lambda row: row["id"]),
    }


def build_reference(cache_dir: Path = DEFAULT_CACHE, download: bool = True) -> dict:
    directory = cache_dir / f"reactome_v{RELEASE}"
    if download:
        with ThreadPoolExecutor(max_workers=4) as pool:
            receipts = list(pool.map(lambda name: download_file(cache_dir, name), FILES))
    else:
        receipts = []
        for name in FILES:
            receipt = json.loads((directory / (name + ".receipt.json")).read_text(encoding="utf-8"))
            if sha256(directory / name) != receipt["sha256"]:
                raise ValueError(f"Cached reference fingerprint mismatch: {name}")
            receipts.append(receipt)
    with zipfile.ZipFile(directory / "ReactomePathways.gmt.zip") as archive:
        if archive.namelist() != ["ReactomePathways.gmt"]:
            raise ValueError("Unexpected Reactome GMT archive members")
        gmt_text = archive.read("ReactomePathways.gmt").decode("utf-8")
    value = normalize_reference(
        gmt_text, (directory / "ReactomePathways.txt").read_text(encoding="utf-8"),
        (directory / "ReactomePathwaysRelation.txt").read_text(encoding="utf-8"),
        (directory / "pathway2summation.txt").read_text(encoding="utf-8"), receipts,
    )
    write_json(cache_dir / f"reactome_human_v{RELEASE}.json", value)
    return value


def select_modules(reference: dict, seed_genes: list[str], min_genes: int = 5,
                   max_genes: int = 200) -> dict:
    """Select leaf-most eligible modules separately per seed, using only local data.

    An ancestor is discarded for a seed when an eligible descendant contains that
    same seed. This avoids pooling broad categories just because the seed occurs
    in both. A pathway needed for another seed can remain selected.
    """
    if min_genes < 1 or max_genes < min_genes:
        raise ValueError("Require 1 <= min_genes <= max_genes")
    seeds = sorted(set(gene.strip() for gene in seed_genes if gene.strip()))
    if not seeds:
        raise ValueError("At least one local seed gene is required")
    rows = []
    for pathway in reference["pathways"]:
        matching = sorted(set(pathway["gene_symbols"]) & set(seeds))
        if not matching:
            continue
        reasons = []
        if not pathway["parent_ids"]:
            reasons.append("top_level_category")
        if pathway["name"].casefold() in BROAD_NAMES:
            reasons.append("broad_category")
        if pathway["id"] == DISEASE_ROOT or DISEASE_ROOT in pathway["ancestor_ids"]:
            reasons.append("disease_category_or_descendant")
        if pathway["gene_count"] < min_genes:
            reasons.append("below_minimum_size")
        if pathway["gene_count"] > max_genes:
            reasons.append("above_maximum_size")
        rows.append({**pathway, "matching_seed_genes": matching,
                     "exclusion_reasons": reasons, "selected_for_seed_genes": []})
    for seed in seeds:
        eligible = [row for row in rows if not row["exclusion_reasons"] and seed in row["matching_seed_genes"]]
        for row in eligible:
            if not any(row["id"] in other["ancestor_ids"] for other in eligible):
                row["selected_for_seed_genes"].append(seed)
    for row in rows:
        if not row["selected_for_seed_genes"] and not row["exclusion_reasons"]:
            row["exclusion_reasons"].append("broader_ancestor_of_eligible_seed_module")
    selected = [row for row in rows if row["selected_for_seed_genes"]]
    excluded = [row for row in rows if not row["selected_for_seed_genes"]]
    genes = sorted({gene for row in selected for gene in row["gene_symbols"]})
    return {
        "schema_version": 1, "reference_database": reference["database"],
        "reference_release": reference["release"], "seed_genes": seeds,
        "policy": {**POLICY, "min_genes": min_genes, "max_genes": max_genes},
        "selected": selected, "excluded": excluded,
        "selected_module_count": len(selected), "excluded_module_count": len(excluded),
        "union_genes": genes, "union_gene_count": len(genes),
        "seed_coverage": {
            seed: {"reference_module_count": sum(seed in row["matching_seed_genes"] for row in rows),
                   "selected_module_count": sum(seed in row["selected_for_seed_genes"] for row in selected)}
            for seed in seeds
        },
        "interpretation": "Associated pathway membership only; not direct interaction, downstream direction, mutation correction, or evidence of therapeutic benefit.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--seed", action="append", default=[])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    reference = build_reference(args.cache, download=not args.offline)
    summary = {"release": reference["release"], "pathway_count": reference["pathway_count"],
               "gmt_pathway_count": reference["gmt_pathway_count"], "gene_universe_count": reference["gene_universe_count"]}
    if args.seed:
        modules = select_modules(reference, args.seed)
        write_json(args.output or args.cache / "reactome_seed_modules.json", modules)
        summary.update({"selected_module_count": modules["selected_module_count"], "union_gene_count": modules["union_gene_count"]})
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
