"""Acquire an unconditioned public ChEMBL medicine/mechanism catalog.

No patient inputs are accepted. ``max_phase=4`` is a historical development
annotation, not verification of current approval, availability, or suitability.
All source records and exact download receipts are retained for audit.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
import urllib.error
import urllib.parse
import urllib.request


BASE = "https://www.ebi.ac.uk/chembl/api/data/"
DOCS = "https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services"
COLLECTIONS = {
    "molecules": ("molecule", "molecules", "molecule_chembl_id", {"max_phase": "4"}),
    "mechanisms": ("mechanism", "mechanisms", "mec_id", {}),
    "human_targets": ("target", "targets", "target_chembl_id", {"tax_id": "9606"}),
    "warnings": ("drug_warning", "drug_warnings", "warning_id", {}),
}


def stamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def digest(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def save(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def request_json(url: str, path: Path, version: str | None, *, reuse: bool = True) -> tuple[dict, dict]:
    """Cache only a validated JSON response, with an exact URL and hash receipt."""
    receipt_path = path.with_name(path.name + ".source.json")
    if reuse and path.exists() and receipt_path.exists():
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if (receipt.get("url") == url and receipt.get("chembl_version") == version
                and receipt.get("sha256") == digest(path)):
            return json.loads(path.read_text(encoding="utf-8")), receipt
    path.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "MVA-public-catalog/0.1", "Accept": "application/json"})
            with urllib.request.urlopen(request, timeout=120) as response:
                raw = response.read()
                expected = response.headers.get("Content-Length")
                if expected is not None and len(raw) != int(expected):
                    raise ValueError("Truncated ChEMBL response")
                data = json.loads(raw.decode("utf-8"))
                if not isinstance(data, dict):
                    raise ValueError("Unexpected ChEMBL JSON response")
                resolved_url = response.geturl()
                last_modified = response.headers.get("Last-Modified")
            temporary = path.with_name(path.name + ".part")
            temporary.write_bytes(raw)
            temporary.replace(path)
            receipt = {"url": url, "resolved_url": resolved_url, "retrieved_at_utc": stamp(),
                       "chembl_version": version, "bytes": len(raw), "sha256": digest(path),
                       "last_modified": last_modified}
            save(receipt_path, receipt)
            return data, receipt
        except (urllib.error.URLError, TimeoutError):
            if attempt == 3:
                raise
            time.sleep(min(2 ** attempt, 8))
    raise AssertionError("Unreachable")


def collect(name: str, directory: Path, version: str, *, fetch=request_json) -> dict:
    """Follow every page; reject missing/duplicate IDs or unstable pagination."""
    endpoint, row_key, id_key, filters = COLLECTIONS[name]
    url = BASE + endpoint + ".json?" + urllib.parse.urlencode({**filters, "limit": "1000", "offset": "0", "order_by": id_key})
    rows, receipts, seen_ids, seen_urls = [], [], set(), set()
    expected_count = None
    while url:
        if url in seen_urls:
            raise ValueError(f"{name}: pagination loop")
        seen_urls.add(url)
        path = directory / "raw" / name / f"page_{len(receipts):04d}.json"
        data, receipt = fetch(url, path, version)
        page, metadata = data.get(row_key), data.get("page_meta")
        if not isinstance(page, list) or not isinstance(metadata, dict):
            raise ValueError(f"{name}: missing records/page metadata")
        total = metadata.get("total_count")
        if not isinstance(total, int) or total < 0:
            raise ValueError(f"{name}: invalid total count")
        if expected_count is None:
            expected_count = total
        if total != expected_count or metadata.get("offset") != len(rows):
            raise ValueError(f"{name}: changed total or discontinuous offset")
        for row in page:
            identifier = row.get(id_key)
            if identifier is None or identifier in seen_ids:
                raise ValueError(f"{name}: absent or repeated ID {identifier}")
            seen_ids.add(identifier)
        rows.extend(page)
        receipts.append({"path": path.relative_to(directory).as_posix(), "row_count": len(page), **receipt})
        next_page = metadata.get("next")
        if next_page:
            next_url = urllib.parse.urljoin(BASE, next_page)
            parsed = urllib.parse.urlparse(next_url)
            if parsed.scheme != "https" or parsed.netloc != "www.ebi.ac.uk" or parsed.path != f"/chembl/api/data/{endpoint}.json":
                raise ValueError(f"{name}: unexpected pagination destination")
            if not page:
                raise ValueError(f"{name}: empty nonterminal page")
            url = next_url
        else:
            url = None
        print(f"ChEMBL {name}: {len(rows):,}/{expected_count:,}", flush=True)
    if len(rows) != expected_count:
        raise ValueError(f"{name}: incomplete acquisition {len(rows)} != {expected_count}")
    result = {"collection": name, "chembl_version": version, "complete": True,
              "expected_count": expected_count, "record_count": len(rows), "pages": receipts, "rows": rows}
    save(directory / (name + ".json"), result)
    return result


def target_mapping(target: dict) -> dict:
    """A complex's components are candidates, not individually proven binders."""
    components = []
    for component in target.get("target_components") or []:
        symbols = sorted({item["component_synonym"] for item in component.get("target_component_synonyms", [])
                          if item.get("syn_type") == "GENE_SYMBOL" and item.get("component_synonym")})
        basis = "ChEMBL GENE_SYMBOL"
        if not symbols:
            symbols = sorted({item["xref_name"] for item in component.get("target_component_xrefs", [])
                              if item.get("xref_src_db") == "HGNC" and item.get("xref_name")})
            basis = "ChEMBL HGNC cross-reference" if symbols else "unmapped"
        components.append({"component_id": component.get("component_id"), "accession": component.get("accession"),
                           "component_type": component.get("component_type"), "relationship": component.get("relationship"),
                           "genes": symbols, "gene_mapping_basis": basis})
    return {"id": target["target_chembl_id"], "name": target.get("pref_name"), "type": target.get("target_type"),
            "tax_id": target.get("tax_id"), "species_group_flag": target.get("species_group_flag"),
            "components": components, "genes": sorted({gene for c in components for gene in c["genes"]})}


def normalize(molecules: list[dict], mechanisms: list[dict], targets: list[dict], warnings: list[dict], source: dict) -> dict:
    """Normalize all phase-four molecules, including those with zero mechanisms."""
    human = {item["target_chembl_id"]: target_mapping(item) for item in targets if str(item.get("tax_id")) == "9606"}
    mechanism_by_molecule, warnings_by_molecule = defaultdict(dict), defaultdict(dict)
    for rows, index, key in ((mechanisms, mechanism_by_molecule, "mec_id"), (warnings, warnings_by_molecule, "warning_id")):
        for row in rows:
            for identifier in {row.get("molecule_chembl_id"), row.get("parent_molecule_chembl_id")} - {None}:
                index[identifier][row[key]] = row
    drugs = []
    seen = set()
    for molecule in sorted(molecules, key=lambda item: item["molecule_chembl_id"]):
        identifier = molecule["molecule_chembl_id"]
        if identifier in seen or float(molecule.get("max_phase") or 0) != 4:
            raise ValueError("Duplicate molecule or non-phase-four member of catalog")
        seen.add(identifier)
        hierarchy = molecule.get("molecule_hierarchy") or {}
        parent = hierarchy.get("parent_chembl_id")
        family = {identifier, parent} - {None}
        found_mechanisms, found_warnings = {}, {}
        for member in family:
            found_mechanisms.update(mechanism_by_molecule.get(member, {}))
            found_warnings.update(warnings_by_molecule.get(member, {}))
        normalized_mechanisms = []
        for mechanism in sorted(found_mechanisms.values(), key=lambda item: item["mec_id"]):
            target_id = mechanism.get("target_chembl_id")
            target = human.get(target_id)
            normalized_mechanisms.append({
                "id": mechanism["mec_id"], "action_type": mechanism.get("action_type"),
                "mechanism_of_action": mechanism.get("mechanism_of_action"), "target_id": target_id,
                "target_name": target["name"] if target else None, "target_type": target["type"] if target else None,
                "target_tax_id": target["tax_id"] if target else None,
                "genes": target["genes"] if target else [], "target_components": target["components"] if target else [],
                "target_mapping_status": "human_component_genes" if target and target["genes"] else "human_no_gene_symbol" if target else "target_not_in_human_reference" if target_id else "no_target_record",
                "gene_specificity": "single_protein" if target and target["type"] == "SINGLE PROTEIN" else "component_membership_only" if target and target["genes"] else "unmapped",
                "direct_interaction": mechanism.get("direct_interaction"), "disease_efficacy": mechanism.get("disease_efficacy"),
                "mechanism_max_phase": mechanism.get("max_phase"), "molecular_mechanism": mechanism.get("molecular_mechanism"),
                "source_molecule_id": mechanism.get("molecule_chembl_id"), "source_parent_id": mechanism.get("parent_molecule_chembl_id"),
                "match_basis": "exact_molecule" if mechanism.get("molecule_chembl_id") == identifier else "shared_parent",
                "refs": mechanism.get("mechanism_refs") or [], "variant_sequence": mechanism.get("variant_sequence"),
                "mechanism_comment": mechanism.get("mechanism_comment"), "selectivity_comment": mechanism.get("selectivity_comment"),
            })
        warning_rows = sorted(found_warnings.values(), key=lambda item: item["warning_id"])
        withdrawal_warnings = [w for w in warning_rows if "withdraw" in str(w.get("warning_type", "")).lower()]
        structures = molecule.get("molecule_structures") or {}
        drugs.append({"id": identifier, "name": molecule.get("pref_name") or identifier, "molecule_type": molecule.get("molecule_type"),
                      "max_phase": float(molecule["max_phase"]), "parent_id": parent, "active_id": hierarchy.get("active_chembl_id"),
                      "withdrawn_flag": molecule.get("withdrawn_flag"), "withdrawal_warning_present": bool(withdrawal_warnings),
                      "black_box_warning": molecule.get("black_box_warning"), "veterinary": molecule.get("veterinary"),
                      "first_approval": molecule.get("first_approval"), "availability_type": molecule.get("availability_type"),
                      "current_regulatory_status": "unverified", "approval_scope": "ChEMBL maximum development phase 4; current jurisdiction, formulation, route and pediatric approval not verified",
                      "oral": molecule.get("oral"), "parenteral": molecule.get("parenteral"), "topical": molecule.get("topical"),
                      "canonical_smiles": structures.get("canonical_smiles"), "standard_inchi_key": structures.get("standard_inchi_key"),
                      "molecule_properties": molecule.get("molecule_properties"), "synonyms": molecule.get("molecule_synonyms") or [],
                      "mechanisms": normalized_mechanisms, "warnings": warning_rows,
                      "human_target_genes": sorted({g for m in normalized_mechanisms for g in m["genes"]}),
                      "source_url": BASE + "molecule/" + identifier + ".json"})
    summary = {"molecule_count": len(drugs), "unique_parent_count": len({d["parent_id"] or d["id"] for d in drugs}),
               "with_mechanisms": sum(bool(d["mechanisms"]) for d in drugs),
               "without_mechanisms": sum(not d["mechanisms"] for d in drugs),
               "with_human_target_genes": sum(bool(d["human_target_genes"]) for d in drugs),
               "without_human_target_genes": sum(not d["human_target_genes"] for d in drugs),
               "human_target_gene_count": len({g for d in drugs for g in d["human_target_genes"]}),
               "withdrawn_flag_true": sum(d["withdrawn_flag"] is True for d in drugs),
               "withdrawal_warning_present": sum(d["withdrawal_warning_present"] for d in drugs),
               "withdrawn_flag_missing": sum(d["withdrawn_flag"] is None for d in drugs),
               "veterinary_flagged": sum(bool(d["veterinary"]) for d in drugs),
               "black_box_warning_flagged": sum(bool(d["black_box_warning"]) for d in drugs),
               "molecule_types": dict(sorted(Counter(str(d["molecule_type"]) for d in drugs).items())),
               "currently_approved_verified": 0}
    return {"schema_version": 1, "source": source, "summary": summary, "limitations": [
        "Not a complete FDA/currently marketed drug list; max_phase is maximum development phase reached in ChEMBL.",
        "Warnings are source-reported and incomplete; no warning record does not establish safety.",
        "Known mechanisms do not enumerate all off-target binding, transcriptomic effects or mutation rescue mechanisms.",
        "Shared-parent joins may transfer mechanisms or warnings between related formulations; match_basis is explicit.",
        "Complex/family component membership does not prove direct binding to each gene product.",
        "A target absent from the human reference is unmapped, not assumed nonhuman or inactive.",
        "No patient data is used for acquisition and no patient data is sent to ChEMBL."], "drugs": drugs}


def build(directory: Path) -> dict:
    acquisition = json.loads((directory / "acquisition.json").read_text(encoding="utf-8"))
    if acquisition.get("status") != "complete":
        raise ValueError("Cannot normalize an incomplete acquisition")
    data = {}
    for name in COLLECTIONS:
        path = directory / (name + ".json")
        if digest(path) != acquisition["collections"][name]["sha256"]:
            raise ValueError(f"Changed cached collection: {name}")
        record = json.loads(path.read_text(encoding="utf-8"))
        if not record.get("complete") or record["chembl_version"] != acquisition["chembl_version"]:
            raise ValueError(f"Incomplete or mismatched release: {name}")
        data[name] = record["rows"]
    source = {"provider": "EMBL-EBI ChEMBL", "chembl_version": acquisition["chembl_version"],
              "chembl_release_date": acquisition["chembl_release_date"], "acquired_at_utc": acquisition["completed_at_utc"],
              "documentation_url": DOCS, "query_policy": "unconditioned phase-four molecule universe, all mechanisms, all human targets, all drug warnings",
              "acquisition_manifest": "acquisition.json", "acquisition_sha256": digest(directory / "acquisition.json"),
              "collections": acquisition["collections"]}
    catalog = normalize(data["molecules"], data["mechanisms"], data["human_targets"], data["warnings"], source)
    save(directory / "catalog.json", catalog)
    save(directory / "catalog_summary.json", {"source": source, "summary": catalog["summary"], "catalog_sha256": digest(directory / "catalog.json")})
    return catalog["summary"]


def acquire(directory: Path) -> dict:
    directory.mkdir(parents=True, exist_ok=True)
    before, _ = request_json(BASE + "status.json", directory / "status_before.json", None, reuse=False)
    if before.get("status") != "UP" or not before.get("chembl_db_version"):
        raise ValueError("ChEMBL service not ready")
    version = before["chembl_db_version"]
    manifest = {"status": "running", "started_at_utc": stamp(), "chembl_version": version,
                "chembl_release_date": before.get("chembl_release_date"), "collections": {}}
    save(directory / "acquisition.json", manifest)
    try:
        with ThreadPoolExecutor(max_workers=4) as pool:
            tasks = {name: pool.submit(collect, name, directory, version) for name in COLLECTIONS}
            for name, task in tasks.items():
                result = task.result()
                manifest["collections"][name] = {"record_count": result["record_count"], "page_count": len(result["pages"]),
                                                 "path": name + ".json", "sha256": digest(directory / (name + ".json"))}
        after, _ = request_json(BASE + "status.json", directory / "status_after.json", None, reuse=False)
        if after.get("chembl_db_version") != version or after.get("chembl_release_date") != before.get("chembl_release_date"):
            raise ValueError("ChEMBL release changed during acquisition")
        manifest.update(status="complete", completed_at_utc=stamp(), release_consistent=True)
        save(directory / "acquisition.json", manifest)
        return build(directory)
    except Exception as error:
        manifest.update(status="failed", error=str(error), failed_at_utc=stamp())
        save(directory / "acquisition.json", manifest)
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["download", "build"])
    parser.add_argument("--directory", type=Path, default=Path("references/drug_screen"))
    args = parser.parse_args()
    summary = acquire(args.directory) if args.action == "download" else build(args.directory)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
