"""Command-line entry point. Python 3.11+, standard library only."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .intake import inspect_data, sha256


def run_enrichment(signature: Path, references: Path, output: Path) -> list:
    from .enrichment import run
    libraries = json.loads((references / "libraries.json").read_text(encoding="utf-8"))
    result = []
    for library in libraries["libraries"]:
        name = library["name"]
        current = run(signature, references / f"{name}.gmt", references / f"{name}.universe.txt", output / name)
        result.append({"library": name, "significant_set_count": current["significant_set_count"], "sets": len(current["results"])})
    return result


def pipeline(args) -> dict:
    from .track1 import run as variants
    from .track2 import run as drugs
    from .report import build
    output = args.out or Path("outputs") / datetime.now(timezone.utc).strftime("run_%Y%m%dT%H%M%SZ")
    if output.exists():
        raise ValueError(f"Run output already exists: {output}. Choose a new --out directory to preserve earlier runs.")
    references = Path("references/ensembl/genes_grch38.json")
    if not references.exists() or not Path("references/enrichment/libraries.json").exists():
        raise ValueError("Public references are missing. Run python -m mva prepare first.")
    output.mkdir(parents=True)
    started = datetime.now(timezone.utc).isoformat()
    code_files = sorted(Path("mva").glob("*.py")) + sorted(Path("config").glob("*.json"))
    manifest = {"version": __version__, "started_at_utc": started, "python": platform.python_version(), "platform": platform.platform(), "command": sys.argv, "files": {p.as_posix(): sha256(p) for p in code_files}, "reference_sha256": sha256(references), "status": "running"}
    reference_inputs = [references, *Path("references/enrichment").glob("*.gmt"), *Path("references/enrichment").glob("*.universe.txt"), Path("references/enrichment/libraries.json")]
    manifest["reference_inputs_sha256"] = {p.as_posix(): sha256(p) for p in reference_inputs}
    manifest_path = output / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    try:
        print("Inspecting local inputs and sampling read-pair structure...", flush=True)
        intake = inspect_data(args.data, output / "intake")
        manifest["vcf_sha256"] = intake["vcf_sha256"]
        manifest["input_vcf"] = intake["vcf"]
        print("Auditing all VCF records and annotating panel variants locally...", flush=True)
        track1 = variants(Path(intake["vcf"]), references, output / "track1")
        print(f"Track 1 complete: {track1['records_seen']:,} records; {track1['hypothesis_count']} hypotheses.", flush=True)
        drugs(Path("config/track2_evidence.json"), output / "track2")
        (output / "track2/evidence_snapshot.json").write_bytes(Path("config/track2_evidence.json").read_bytes())
        run_enrichment(Path("config/mva_literature_signature.json"), Path("references/enrichment"), output / "enrichment")
        result = build(output, Path("config/track2_evidence.json"))
        manifest["status"] = "completed"
        manifest["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
        manifest["summary"] = result
        manifest["outputs_sha256"] = {p.relative_to(output).as_posix(): sha256(p) for p in sorted(output.rglob("*")) if p.is_file() and p != manifest_path}
    except Exception as error:
        manifest["status"] = "failed"
        manifest["error"] = str(error)
        raise
    finally:
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="MVA local research workflow")
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare", help="Download only public reference data; no patient uploads")
    prepare.add_argument("--genes", type=Path, default=Path("config/genes.json"))
    runner = commands.add_parser("run", help="Run both tracks and enrichment into a new output directory")
    runner.add_argument("--data", type=Path, default=Path("data"))
    runner.add_argument("--out", type=Path)
    inventory = commands.add_parser("inventory")
    inventory.add_argument("--data", type=Path, default=Path("data"))
    inventory.add_argument("--out", type=Path, default=Path("outputs/intake"))
    track1 = commands.add_parser("track1")
    track1.add_argument("--vcf", type=Path, required=True)
    track1.add_argument("--reference", type=Path, default=Path("references/ensembl/genes_grch38.json"))
    track1.add_argument("--out", type=Path, default=Path("outputs/track1"))
    track1.add_argument("--max-records", type=int, help="Smoke test only; output is explicitly incomplete")
    track2 = commands.add_parser("track2")
    track2.add_argument("--config", type=Path, default=Path("config/track2_evidence.json"))
    track2.add_argument("--out", type=Path, default=Path("outputs/track2"))
    enrichment = commands.add_parser("enrich")
    enrichment.add_argument("--signature", type=Path, default=Path("config/mva_literature_signature.json"))
    enrichment.add_argument("--references", type=Path, default=Path("references/enrichment"))
    enrichment.add_argument("--out", type=Path, default=Path("outputs/enrichment"))
    review = commands.add_parser("review")
    review.add_argument("--out", type=Path, required=True)
    check = commands.add_parser("validate-submission")
    check.add_argument("path", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "prepare":
            from .references import prepare_genes, prepare_enrichment
            result = {"genes": prepare_genes(args.genes, Path("references/ensembl")), "enrichment": prepare_enrichment(Path("references/enrichment"))}
        elif args.command == "run":
            result = pipeline(args)
        elif args.command == "inventory":
            result = inspect_data(args.data, args.out)
        elif args.command == "track1":
            from .track1 import run
            result = run(args.vcf, args.reference, args.out, args.max_records)
            result.pop("qc", None)
        elif args.command == "track2":
            from .track2 import run
            full = run(args.config, args.out)
            result = {"counts": full["counts"], "outputs": full["outputs"]}
        elif args.command == "enrich":
            result = run_enrichment(args.signature, args.references, args.out)
        elif args.command == "review":
            from .report import build
            result = build(args.out, Path("config/track2_evidence.json"))
        else:
            from .submission import validate
            result = validate(args.path)
        print(json.dumps(result, indent=2, allow_nan=False))
        return 1 if isinstance(result, dict) and result.get("valid") is False else 0
    except (ValueError, OSError, KeyError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
