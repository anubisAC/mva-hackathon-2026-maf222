"""Local, direction-aware overrepresentation analysis with an explicit universe.

No gene identifiers leave the machine. Gene symbols are case sensitive: only
surrounding whitespace and exact duplicates are removed. A literature-derived
signature can nominate a mechanism but cannot establish expression direction.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable


FDR_THRESHOLD = 0.05


def _integer(value: int, label: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{label} must be an integer")


def _log_combination(total: int, selected: int) -> float:
    return (
        math.lgamma(total + 1)
        - math.lgamma(selected + 1)
        - math.lgamma(total - selected + 1)
    )


def _log_add(left: float, right: float) -> float:
    if left == -math.inf:
        return right
    if right == -math.inf:
        return left
    high, low = max(left, right), min(left, right)
    return high + math.log1p(math.exp(low - high))


def hypergeometric_tail(
    overlap: int, population_size: int, set_size: int, signature_size: int
) -> float:
    """Return P[X >= overlap] without subtracting a CDF from one.

    X is hypergeometric(N=population_size, K=set_size, n=signature_size).
    Log-space accumulation and a recurrence avoid factorial overflow and loss
    of tiny upper-tail probabilities through catastrophic cancellation.
    """
    for label, value in (
        ("overlap", overlap),
        ("population_size", population_size),
        ("set_size", set_size),
        ("signature_size", signature_size),
    ):
        _integer(value, label)
    if population_size <= 0:
        raise ValueError("population_size must be positive")
    if not 0 <= set_size <= population_size:
        raise ValueError("set_size must be between zero and population_size")
    if not 0 <= signature_size <= population_size:
        raise ValueError("signature_size must be between zero and population_size")
    if overlap < 0:
        raise ValueError("overlap must be nonnegative")

    low = max(0, signature_size - (population_size - set_size))
    high = min(set_size, signature_size)
    if overlap <= low:
        return 1.0
    if overlap > high:
        return 0.0

    log_mass = (
        _log_combination(set_size, overlap)
        + _log_combination(population_size - set_size, signature_size - overlap)
        - _log_combination(population_size, signature_size)
    )
    log_tail = -math.inf
    for hits in range(overlap, high + 1):
        log_tail = _log_add(log_tail, log_mass)
        if hits < high:
            log_mass += (
                math.log(set_size - hits)
                + math.log(signature_size - hits)
                - math.log(hits + 1)
                - math.log(population_size - set_size - signature_size + hits + 1)
            )
    return min(1.0, math.exp(log_tail))


def benjamini_hochberg(p_values: Iterable[float]) -> list[float]:
    """Adjust every supplied test, retaining original order and tied values."""
    values = list(p_values)
    if any(
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or not 0 <= value <= 1
        for value in values
    ):
        raise ValueError("All p-values must be finite numbers between zero and one")
    adjusted = [1.0] * len(values)
    order = sorted(range(len(values)), key=values.__getitem__)
    running_minimum = 1.0
    for rank in range(len(order), 0, -1):
        index = order[rank - 1]
        running_minimum = min(running_minimum, values[index] * len(order) / rank)
        adjusted[index] = running_minimum
    return adjusted


def _symbols(values: Any, label: str) -> tuple[list[str], int]:
    if not isinstance(values, list):
        raise ValueError(f"{label} must be a list of gene-symbol strings")
    unique: list[str] = []
    seen: set[str] = set()
    duplicates = 0
    for value in values:
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{label} contains an empty or non-string gene symbol")
        symbol = value.strip()
        if any(character.isspace() for character in symbol):
            raise ValueError(f"{label} contains a gene symbol with internal whitespace: {symbol!r}")
        if symbol in seen:
            duplicates += 1
        else:
            seen.add(symbol)
            unique.append(symbol)
    return unique, duplicates


def _read_signature(path: Path) -> tuple[dict[str, Any], int]:
    try:
        signature = json.loads(path.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as error:
        raise ValueError(f"Invalid signature JSON: {error}") from error
    if not isinstance(signature, dict):
        raise ValueError("Signature must be a JSON object")
    for name in ("name", "species", "context", "source", "direction", "evidence_type"):
        if not isinstance(signature.get(name), str) or not signature[name].strip():
            raise ValueError(f"Signature requires a nonempty {name!r} string")
        signature[name] = signature[name].strip()
    if signature["direction"] not in {"up", "down", "undirected"}:
        raise ValueError("Signature direction must be up, down, or undirected")
    if signature["evidence_type"] not in {"measured_expression", "literature_gene_set"}:
        raise ValueError("Signature evidence_type must be measured_expression or literature_gene_set")
    signature["genes"], duplicates = _symbols(signature.get("genes"), "Signature genes")
    signature["direction_status"] = (
        "hypothesis_only_unvalidated_direction"
        if signature["evidence_type"] == "literature_gene_set" and signature["direction"] != "undirected"
        else "undirected"
        if signature["direction"] == "undirected"
        else "reported_measured_expression"
    )
    return signature, duplicates


def _read_universe(path: Path) -> tuple[list[str], int]:
    entries = [
        line.strip()
        for line in path.read_text(encoding="utf-8-sig").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    symbols, duplicates = _symbols(entries, "Universe")
    if not symbols:
        raise ValueError("Universe is empty; provide the genes eligible to enter this signature")
    return symbols, duplicates


def _read_gmt(path: Path) -> list[dict[str, Any]]:
    gene_sets: list[dict[str, Any]] = []
    names: set[str] = set()
    for line_number, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 3 or not fields[0].strip():
            raise ValueError(f"GMT line {line_number} needs a set name, description, and at least one gene")
        name = fields[0].strip()
        if name in names:
            raise ValueError(f"Duplicate GMT set name: {name!r}")
        genes, duplicates = _symbols(fields[2:], f"GMT set {name!r}")
        names.add(name)
        gene_sets.append({"set_name": name, "description": fields[1].strip(), "genes": genes, "duplicates_removed": duplicates})
    if not gene_sets:
        raise ValueError("GMT library contains no gene sets")
    return gene_sets


def _summary(report: dict[str, Any]) -> str:
    signature = report["signature"]
    lines = [
        f"Enrichment: {signature['name']}",
        f"Species: {signature['species']}",
        f"Context: {signature['context']}",
        f"Direction: {signature['direction']} ({signature['direction_status']})",
        f"Evidence: {signature['evidence_type']}; source: {signature['source']}",
        f"Signature genes: {len(signature['genes'])}; explicit universe: {report['universe']['gene_count']}",
        f"Tests: {report['library']['tested_set_count']}; BH FDR < 0.05: {report['significant_set_count']}",
        "",
        "Top sets (ranked by FDR, then raw p-value):",
    ]
    for result in report["results"][:10]:
        lines.append(
            f"{result['set_name']}: overlap={result['overlap']}/{result['set_size']}; "
            f"p={result['p_value']:.6g}; FDR={result['fdr']:.6g}; "
            f"genes={', '.join(result['overlap_genes']) or '(none)'}"
        )
    lines.extend(["", "Interpretation and limitations:"])
    lines.extend(f"- {warning}" for warning in report["warnings"])
    return "\n".join(lines) + "\n"


def run(signature_path: Path, gmt_path: Path, universe_path: Path, output_dir: Path) -> dict[str, Any]:
    """Validate local inputs, test all GMT sets, and write auditable results.

    The explicit universe must represent genes eligible for selection (for
    measured data, typically all reliably measured genes). Signature genes
    outside this universe are errors. GMT sets are intersected with the
    universe; empty and zero-overlap sets remain in the correction family.
    Species compatibility cannot be inferred from GMT and must be checked by
    the analyst. A source string is provenance supplied by the user, not proof
    that measurements or literature have been independently validated.
    """
    signature_path, gmt_path = Path(signature_path), Path(gmt_path)
    universe_path, output_dir = Path(universe_path), Path(output_dir)
    signature, signature_duplicates = _read_signature(signature_path)
    universe_genes, universe_duplicates = _read_universe(universe_path)
    gene_sets = _read_gmt(gmt_path)
    universe = set(universe_genes)
    selected = set(signature["genes"])
    outside = sorted(selected - universe)
    if outside:
        raise ValueError(
            f"Signature has {len(outside)} gene(s) outside the explicit universe: "
            + ", ".join(outside[:15])
            + ". Check identifier/species compatibility or supply the correct background."
        )

    warnings = [
        "Exploratory analysis: only BH FDR < 0.05 is flagged; a small raw p-value alone is not corrected significance. "
        "FDR is calculated across every set in this run, including empty and zero-overlap sets; it does not correct "
        "across other signatures, libraries, or repeated exploratory runs. BH control relies on its usual dependence assumptions.",
        "Overrepresentation is an association and does not establish regulator activity, causal mechanism, therapeutic benefit, "
        "or expression reversal. Disease-gene enrichment is not evidence of expression reversal.",
        "The universe must represent eligible genes, not just the selected hits. Results depend on background choice, "
        "library version, and selection method; literature selection can introduce publication bias.",
        "Symbols are case sensitive and only trimmed and deduplicated. Check species, identifier vocabulary, assay context, "
        "and library compatibility; no ortholog or alias mapping is performed.",
    ]
    if signature["evidence_type"] == "literature_gene_set":
        warnings.append(
            "This signature is a literature gene set, not measured expression. Any supplied up/down label is a "
            "hypothesis only and cannot support a claim about observed expression direction or reversal."
        )
    else:
        warnings.append(
            "Measured-expression status and direction are supplied in the signature metadata; this module does not "
            "verify the underlying experiment, differential-expression model, or source."
        )
    if not selected:
        warnings.append("The signature is empty: every set has p-value 1 and FDR 1; there is no enrichment evidence.")

    results: list[dict[str, Any]] = []
    for gene_set in gene_sets:
        members = set(gene_set["genes"]) & universe
        shared = sorted(selected & members)
        expected = len(selected) * len(members) / len(universe)
        results.append({
            "signature_name": signature["name"],
            "species": signature["species"],
            "context": signature["context"],
            "direction": signature["direction"],
            "direction_status": signature["direction_status"],
            "evidence_type": signature["evidence_type"],
            "signature_source": signature["source"],
            "set_name": gene_set["set_name"],
            "description": gene_set["description"],
            "set_size": len(members),
            "original_set_size": len(gene_set["genes"]),
            "genes_outside_universe": len(set(gene_set["genes"]) - universe),
            "duplicates_removed": gene_set["duplicates_removed"],
            "signature_size": len(selected),
            "universe_size": len(universe),
            "overlap": len(shared),
            "expected_overlap": expected,
            "fold_enrichment": len(shared) / expected if expected else None,
            "p_value": hypergeometric_tail(len(shared), len(universe), len(members), len(selected)),
            "overlap_genes": shared,
        })
    for result, fdr in zip(results, benjamini_hochberg(result["p_value"] for result in results)):
        result["fdr"] = fdr
        result["significant_fdr_0_05"] = fdr < FDR_THRESHOLD
    results.sort(key=lambda result: (result["fdr"], result["p_value"], result["set_name"]))
    if selected and not any(result["overlap"] for result in results):
        warnings.append("No selected gene overlaps any tested set. Check library coverage and identifiers before interpreting this null result.")
    excluded = sum(result["genes_outside_universe"] for result in results)
    if excluded:
        warnings.append(f"Excluded {excluded} GMT gene memberships outside the explicit universe; counts may include a gene in multiple sets.")
    if any(result["p_value"] == 0.0 and result["overlap"] for result in results):
        warnings.append("At least one p-value underflowed the floating-point range; a displayed zero is a numerical limit, not an exact probability of zero.")

    report: dict[str, Any] = {
        "schema_version": 1,
        "method": "One-sided hypergeometric overrepresentation with Benjamini-Hochberg FDR",
        "signature": signature,
        "signature_path": str(signature_path.resolve()),
        "signature_duplicates_removed": signature_duplicates,
        "universe": {"path": str(universe_path.resolve()), "gene_count": len(universe), "duplicates_removed": universe_duplicates},
        "library": {"path": str(gmt_path.resolve()), "tested_set_count": len(results), "excluded_gene_memberships": excluded},
        "significant_set_count": sum(result["significant_fdr_0_05"] for result in results),
        "fdr_threshold": FDR_THRESHOLD,
        "results": results,
        "warnings": warnings,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "json": output_dir / "enrichment_results.json",
        "csv": output_dir / "enrichment_results.csv",
        "summary": output_dir / "enrichment_summary.txt",
    }
    report["outputs"] = {name: str(path.resolve()) for name, path in paths.items()}
    paths["json"].write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    fields = list(results[0])
    with paths["csv"].open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for result in results:
            writer.writerow({**result, "overlap_genes": ";".join(result["overlap_genes"])})
    paths["summary"].write_text(_summary(report), encoding="utf-8")
    return report
