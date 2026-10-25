"""Compare public BUB1B-loss modules with cell-cycle masking and matched perturbations.

Offline analysis of thresholded signed gene lists. No patient data or network I/O.
Core statistics use the standard library; --plot uses the optional local matplotlib.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from .drug_signatures import read_pairs, file_hash
from .enrichment import hypergeometric_tail, benjamini_hochberg


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def table(path, rows):
    if not rows:
        Path(path).write_text("", encoding="utf-8")
        return
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v
                             for k, v in row.items()})


def symbols(mask, index):
    return sorted(g for g, i in index.items() if mask & (1 << i))


def bits(genes, index):
    return sum(1 << index[g] for g in set(genes) if g in index)


def metadata(pair, source):
    fields = pair["fields"] + [""] * 9
    return {"source": source, "gene": fields[0], "cell": fields[1], "cell_type": fields[2],
            "time": fields[3], "perturbation": fields[5], "method": fields[6],
            "batch": fields[7] if source == "LINCS_XPR" else "", "term": pair["term"]}


def stratum(meta):
    return tuple(meta[k] for k in ("source", "cell", "cell_type", "time", "perturbation", "method"))


def score(pair, module, universe, excluded=0, target=0, minimum=5):
    """Mask every part of the test; target removal prevents self-hit enrichment."""
    background = universe & ~excluded & ~target
    genes = module & background
    up, down = pair["up"] & background, pair["down"] & background
    n, m, nu, nd = background.bit_count(), genes.bit_count(), up.bit_count(), down.bit_count()
    ku, kd = (up & genes).bit_count(), (down & genes).bit_count()
    eligible = n > 0 and m >= minimum and nu > 0 and nd > 0
    return {"eligible": eligible, "background_genes": n, "module_genes": m,
            "up_list_size": nu, "down_list_size": nd, "up_hits": ku, "down_hits": kd,
            "score_pp": 100 * (ku / nu - kd / nd) if eligible else None,
            "up_mask": up & genes, "down_mask": down & genes}


def empirical_rank(value, reference, minimum=20):
    """Inclusive ties, plus-one adjustment, and explicit unassessable small pools."""
    valid = [x for x in reference if x is not None]
    if value is None or len(valid) < minimum:
        return {"n_reference_genes": len(valid), "p": None, "reference_median": None,
                "excess_over_reference_pp": None, "percentile_midrank": None}
    n = len(valid)
    upper = (1 + sum(x >= value for x in valid)) / (n + 1)
    lower = (1 + sum(x <= value for x in valid)) / (n + 1)
    center = statistics.median(valid)
    return {"n_reference_genes": n, "p": min(1.0, 2 * min(upper, lower)),
            "reference_median": center, "excess_over_reference_pp": value - center,
            "percentile_midrank": (sum(x < value for x in valid) + .5 * sum(x == value for x in valid)) / n}


def collapse_gene_scores(records):
    values = defaultdict(list)
    for gene, value in records:
        if value is not None:
            values[gene].append(value)
    return {gene: statistics.median(v) for gene, v in values.items()}


def adjust(rows):
    eligible = [row for row in rows if row.get("p") is not None]
    for row in rows:
        row["q"] = None
    for row, q in zip(eligible, benjamini_hochberg(r["p"] for r in eligible)):
        row["q"] = q
    return len(eligible)


def consistency(rows):
    """One score per cell; no pooling repeated contexts as independent replicates."""
    valid = [r for r in rows if r["score_pp"] is not None]
    scores = {r["cell"]: r["score_pp"] for r in valid}
    if len(scores) != len(valid):
        raise ValueError("Multiple strata per cell require explicit handling")
    loo = {cell: statistics.median([v for c, v in scores.items() if c != cell])
           for cell in scores if len(scores) > 1}
    positive = sum(v > 0 for v in scores.values())
    negative = sum(v < 0 for v in scores.values())
    sig_up = sum(r["score_pp"] > 0 and r.get("all_q") is not None and r["all_q"] <= .05 for r in valid)
    sig_down = sum(r["score_pp"] < 0 and r.get("all_q") is not None and r["all_q"] <= .05 for r in valid)
    needed = math.ceil(2 * len(valid) / 3)
    direction = ("up" if positive >= needed and sig_up >= 2 else
                 "down" if negative >= needed and sig_down >= 2 else "none")
    return {"assessable_cells": len(valid), "positive_cells": positive, "negative_cells": negative,
            "zero_cells": len(valid) - positive - negative,
            "median_score_pp": statistics.median(scores.values()) if scores else None,
            "per_cell_scores": scores, "all_benchmark_q05_up_cells": sig_up,
            "all_benchmark_q05_down_cells": sig_down,
            "cycle_benchmark_q05_cells": sum(r.get("cycle_q") is not None and r["cycle_q"] <= .05 for r in valid),
            "leave_one_cell_out_medians": loo,
            "leave_one_out_direction_consistent": bool(loo) and (all(v > 0 for v in loo.values()) or all(v < 0 for v in loo.values())),
            "exploratory_convergence": direction}


def run(config_path, reference_path, library_dir, output, plot=False):
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Choose a fresh output directory")
    output.mkdir(parents=True, exist_ok=True)
    policy = read(config_path)
    reference = read(reference_path)
    if reference["release"] != policy["reactome_release"]:
        raise ValueError("Wrong Reactome release")
    paths = {"policy": Path(config_path), "reactome": Path(reference_path),
             "LINCS_XPR": Path(library_dir) / "lincs-l1000-xpr.gmt.gz",
             "CM4AI": Path(library_dir) / "cm4ai.gmt.gz"}
    inputs = {k: {"path": str(p.resolve()), "sha256": file_hash(p)} for k, p in paths.items()}
    for source in ("LINCS_XPR", "CM4AI"):
        receipt = read(paths[source].with_name(paths[source].name + ".receipt.json"))
        if inputs[source]["sha256"] != receipt["sha256"]:
            raise ValueError("Public source differs from acquisition receipt")
        inputs[source]["source_url"] = receipt["url"]
    code_paths = [Path(__file__), Path(__file__).with_name("drug_signatures.py"), Path(__file__).with_name("enrichment.py")]
    code_hashes = {str(p.resolve()): file_hash(p) for p in code_paths}
    write(output / "policy.json", policy)
    write(output / "started.json", {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs,
                                    "code_sha256": code_hashes, "scope": "public perturbation libraries only"})
    index, groups, audits = {}, defaultdict(list), {}
    gene = policy["gene"]
    for source in ("LINCS_XPR", "CM4AI"):
        print(f"Reading public {source} library", flush=True)
        pairs, _, audits[source] = read_pairs(paths[source], index)
        for pair in pairs:
            pair["meta"] = metadata(pair, source)
        wanted = {stratum(p["meta"]) for p in pairs if p["meta"]["gene"] == gene}
        for pair in pairs:
            key = stratum(pair["meta"])
            if key in wanted:
                groups[key].append(pair)
    pathways = reference["pathways"]
    if not isinstance(pathways, dict):
        pathways = {p["id"]: p for p in pathways}
    root = policy["cycle_mask_root"]
    cycle_paths = [p for p in pathways.values() if p["id"] == root or root in p["ancestor_ids"]]
    cycle_genes = {g for p in cycle_paths for g in p["gene_symbols"]}
    cycle_mask = bits(cycle_genes, index)
    write(output / "cell_cycle_mask.json", {"root": root, "pathway_ids": sorted(p["id"] for p in cycle_paths),
                                          "genes": sorted(cycle_genes), "gene_count": len(cycle_genes),
                                          "observed_in_libraries": cycle_mask.bit_count()})
    modules = []
    for module in policy["modules"]:
        p = pathways[module["id"]]
        modules.append({**module, "name": p["name"], "genes": p["gene_symbols"],
                        "source_url": p["source_url"], "mask": bits(p["gene_symbols"], index)})
    write(output / "modules.json", [{k: v for k, v in m.items() if k != "mask"} for m in modules])
    contexts, coverage, overrep, cell_rows, benchmark_rows, batch_rows, qc, bg_rows = [], [], [], [], [], [], [], []
    controls_output = output / "matched_gene_scores.csv.gz"
    with gzip.open(controls_output, "wt", encoding="utf-8", newline="") as stream:
        fields = ["source", "cell", "time", "view", "module_id", "gene", "score_pp", "cell_cycle_gene"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for key, pairs in sorted(groups.items()):
            source, cell, cell_type, time, perturbation, method = key
            common = {"source": source, "cell": cell, "cell_type": cell_type, "time": time,
                      "perturbation": perturbation, "method": method}
            print(f"Comparing {source} {cell}: {len(pairs)} public contexts", flush=True)
            queries = [p for p in pairs if p["meta"]["gene"] == gene]
            universe = 0
            for pair in pairs:
                universe |= pair["up"] | pair["down"]
            target_masks = {p["meta"]["gene"]: bits([p["meta"]["gene"]], index) for p in pairs}
            bg_rows.append({**common, "genes": symbols(universe, index)})
            qc.append({**common, "query_contexts": len(queries), "all_contexts": len(pairs),
                       "distinct_perturbed_genes": len(target_masks), "background_genes": universe.bit_count(),
                       "query_batches": sorted({p["meta"]["batch"] for p in queries})})
            for view in policy["views"]:
                excluded = cycle_mask if view == "cell_cycle_removed" else 0
                for query in queries:
                    u, d = query["up"], query["down"]
                    contexts.append({**common, "view": view, "term": query["term"],
                                     "up_genes": symbols(u & ~excluded & ~target_masks[gene], index),
                                     "down_genes": symbols(d & ~excluded & ~target_masks[gene], index),
                                     "cycle_up_removed": (u & cycle_mask).bit_count(),
                                     "cycle_down_removed": (d & cycle_mask).bit_count(),
                                     "original_up_count": u.bit_count(), "original_down_count": d.bit_count()})
                for module in modules:
                    base = {**common, "view": view, "module_id": module["id"],
                            "module_label": module["label"], "domain": module["domain"]}
                    m = module["mask"]
                    stats = [score(p, m, universe, excluded, target_masks[p["meta"]["gene"]],
                                   policy["minimum_module_genes_in_background"]) for p in pairs]
                    lookup = {p["term"]: s for p, s in zip(pairs, stats)}
                    query_first = lookup[queries[0]["term"]]
                    coverage.append({**base, "total_module_genes": len(module["genes"]),
                                     "module_genes_in_background": query_first["module_genes"],
                                     "background_genes": query_first["background_genes"],
                                     "cycle_module_genes_removed": (m & universe & cycle_mask).bit_count(),
                                     "eligible": query_first["eligible"]})
                    for query in queries:
                        s = lookup[query["term"]]
                        for direction in ("up", "down"):
                            hits, size = s[direction + "_hits"], s[direction + "_list_size"]
                            expected = size * s["module_genes"] / s["background_genes"] if s["background_genes"] else 0
                            overrep.append({**base, "term": query["term"], "direction": direction,
                                            "eligible": s["eligible"], "hits": hits, "expected_hits": expected,
                                            "fold_enrichment": hits / expected if expected else None,
                                            "signature_size": size, "module_size": s["module_genes"],
                                            "background_size": s["background_genes"],
                                            "genes": symbols(s[direction + "_mask"], index),
                                            "p": hypergeometric_tail(hits, s["background_genes"], s["module_genes"], size) if s["eligible"] else None})
                    gene_scores = collapse_gene_scores([(p["meta"]["gene"], s["score_pp"]) for p, s in zip(pairs, stats)])
                    for target_gene, value in sorted(gene_scores.items()):
                        writer.writerow({"source": source, "cell": cell, "time": time, "view": view,
                                         "module_id": module["id"], "gene": target_gene, "score_pp": value,
                                         "cell_cycle_gene": target_gene in cycle_genes})
                    value = gene_scores.get(gene)
                    cell_row = {**base, "score_pp": value, "query_contexts": len(queries),
                                "query_context_scores": [lookup[p["term"]]["score_pp"] for p in queries]}
                    cell_rows.append(cell_row)
                    for kind in ("all", "cycle"):
                        reference_values = [v for g, v in gene_scores.items() if g != gene and (kind == "all" or g in cycle_genes)]
                        benchmark_rows.append({**base, "benchmark": kind, "score_pp": value,
                                               **empirical_rank(value, reference_values, policy["minimum_benchmark_genes"])})
                    if source == "LINCS_XPR":
                        for query in queries:
                            batch = query["meta"]["batch"]
                            batch_scores = collapse_gene_scores([(p["meta"]["gene"], s["score_pp"])
                                                                for p, s in zip(pairs, stats) if p["meta"]["batch"] == batch])
                            for kind in ("all", "cycle"):
                                values = [v for g, v in batch_scores.items() if g != gene and (kind == "all" or g in cycle_genes)]
                                batch_rows.append({**base, "term": query["term"], "batch": batch, "benchmark": kind,
                                                   "score_pp": lookup[query["term"]]["score_pp"],
                                                   **empirical_rank(lookup[query["term"]]["score_pp"], values,
                                                                    policy["minimum_benchmark_genes"])})
    family_sizes = {"context_overrepresentation": adjust(overrep), "cell_benchmarks": adjust(benchmark_rows),
                    "batch_benchmarks": adjust(batch_rows)}
    blookup = {(r["source"], r["cell"], r["view"], r["module_id"], r["benchmark"]): r for r in benchmark_rows}
    for row in cell_rows:
        for kind in ("all", "cycle"):
            b = blookup[(row["source"], row["cell"], row["view"], row["module_id"], kind)]
            for k in ("p", "q", "n_reference_genes", "excess_over_reference_pp", "reference_median", "percentile_midrank"):
                row[kind + "_" + k] = b[k]
    cross = []
    for module in modules:
        for view in policy["views"]:
            for source in ("LINCS_XPR", "CM4AI"):
                chosen = [r for r in cell_rows if r["module_id"] == module["id"] and r["view"] == view and r["source"] == source]
                cross.append({"source": source, "view": view, "module_id": module["id"],
                              "module_label": module["label"], "domain": module["domain"], **consistency(chosen)})
    outputs = {"query_contexts.json": contexts, "module_coverage.json": coverage,
               "context_enrichment.json": overrep, "cell_scores.json": cell_rows,
               "cell_benchmarks.json": benchmark_rows, "batch_benchmarks.json": batch_rows,
               "cross_cell_summary.json": cross, "observed_backgrounds.json": bg_rows,
               "library_qc.json": {"libraries": audits, "strata": qc}}
    for name, value in outputs.items():
        write(output / name, value)
    for name, rows in [("cell_scores.csv", cell_rows), ("context_enrichment.csv", overrep),
                       ("cross_cell_summary.csv", cross), ("batch_benchmarks.csv", batch_rows)]:
        table(output / name, rows)
    summary = {"status": "complete", "gene": gene, "completed_utc": datetime.now(timezone.utc).isoformat(),
               "query_contexts": sum(r["query_contexts"] for r in qc),
               "source_cells": {source: sorted({r["cell"] for r in qc if r["source"] == source}) for source in audits},
               "modules": len(modules), "cycle_mask_genes": len(cycle_genes), "test_family_sizes": family_sizes,
               "masked_convergence": [r for r in cross if r["source"] == "LINCS_XPR" and r["view"] == "cell_cycle_removed" and r["exploratory_convergence"] != "none"],
               "context_q05_tests": sum(r.get("q") is not None and r["q"] <= .05 for r in overrep),
               "cell_benchmark_q05_tests": sum(r.get("q") is not None and r["q"] <= .05 for r in benchmark_rows),
               "batch_benchmark_q05_tests": sum(r.get("q") is not None and r["q"] <= .05 for r in batch_rows)}
    write(output / "summary.json", summary)
    if plot:
        render_plot(output, modules, cell_rows)
    for path, expected in code_hashes.items():
        if file_hash(Path(path)) != expected:
            raise RuntimeError("Analysis code changed while running")
    write(output / "run_manifest.json", {"status": "complete", "inputs": inputs, "code_sha256": code_hashes,
                                        "completed_utc": datetime.now(timezone.utc).isoformat(),
                                        "outputs_sha256": {p.name: file_hash(p) for p in sorted(output.iterdir()) if p.is_file()}})
    return summary


def render_plot(output, modules, rows):
    import sys
    dependencies = Path(".analysis_dependencies").resolve()
    if dependencies.exists():
        sys.path.insert(0, str(dependencies))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    selected = [m for m in modules if m["domain"] != "control"]
    cells = sorted({r["cell"] for r in rows if r["source"] == "LINCS_XPR"}) + ["KOLF2.1J"]
    lookup = {(r["view"], r["module_id"], r["cell"]): r for r in rows}
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 8.2), sharey=True)
    values = [abs(r["score_pp"]) for r in rows if r["score_pp"] is not None and r["domain"] != "control"]
    limit = max(values) if values else 1
    for ax, view, title in zip(axes, ("full", "cell_cycle_removed"), ("All retained genes", "Cell-cycle genes removed")):
        data = np.array([[lookup[(view, m["id"], c)]["score_pp"] if lookup[(view, m["id"], c)]["score_pp"] is not None else np.nan for c in cells] for m in selected])
        im = ax.imshow(data, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
        for i, m in enumerate(selected):
            for j, c in enumerate(cells):
                r = lookup[(view, m["id"], c)]
                q = r["all_q"]
                if q is not None and q <= .05:
                    ax.text(j, i, "*", ha="center", va="center", color="black", fontsize=14, weight="bold")
        ax.set_xticks(range(len(cells)), cells, rotation=45, ha="right")
        ax.set_yticks(range(len(selected)), [m["label"] for m in selected], fontsize=9)
        ax.axvline(len(cells)-1.5, color="black", linewidth=1.5)
        ax.set_title(title, fontsize=13, pad=13)
        ax.set_xticks(np.arange(-.5, len(cells), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(selected), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=.6)
        ax.tick_params(which="minor", bottom=False, left=False)
    fig.subplots_adjust(left=.28, right=.92, bottom=.23, top=.84, wspace=.12)
    cax = fig.add_axes([.94, .28, .015, .48])
    fig.colorbar(im, cax=cax, label="Up-minus-down module membership (percentage points)")
    fig.suptitle("Public BUB1B-loss signatures vary across cell models", fontsize=16, x=.61, y=.96)
    fig.text(.28, .89, "Six LINCS cancer cell lines (96 h knockout); separate hiPSC model at right (144 h knockdown)", fontsize=10)
    fig.text(.28, .065, "Red: greater up-list representation. Blue: greater down-list representation.\n* q <= 0.05 versus other gene perturbations in the same cell/time/method.\nColors describe signed gene-list membership, not measured pathway activity or dysfunction.", fontsize=10, linespacing=1.6)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(output / ("module_comparison." + ext), dpi=180)
    write(output / "plot_environment.json", {"matplotlib": matplotlib.__version__, "numpy": np.__version__})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/bub1b_loss_comparison.json"))
    parser.add_argument("--reference", type=Path, default=Path("references/drug_screen/reactome_human_v97.json"))
    parser.add_argument("--libraries", type=Path, default=Path("references/drug_screen/signatures"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    result = run(args.config, args.reference, args.libraries, args.output, args.plot)
    print(json.dumps({k: v for k, v in result.items() if k != "masked_convergence"}, indent=2))


if __name__ == "__main__":
    main()
