"""Offline, auditable Track 2 hypothesis triage; no efficacy prediction or prescribing.

The JSON evidence ledger is deliberately separate from the decision rules. Nothing
in this module makes network requests or reads individual-level genomic data.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from urllib.parse import urlparse


REQUIRED_GATES = (
    "variant_mechanism",
    "pathway_activation",
    "functional_rescue",
    "genome_safety",
    "exposure_feasibility",
)
STATES = ("accepted", "deferred", "rejected")
DISCLAIMER = (
    "Research hypotheses only. Accepted means eligible for further experimental "
    "prioritization, not effective, safe, or recommended for a person. No clinical "
    "efficacy probabilities or medication doses are generated."
)


def _require(value: bool, message: str) -> None:
    if not value:
        raise ValueError(message)


def _references(refs: list, sources: dict, location: str) -> None:
    _require(isinstance(refs, list) and bool(refs), f"{location}: evidence references required")
    for ref in refs:
        _require(ref in sources, f"{location}: unknown source {ref!r}")


def validate_config(config: dict) -> None:
    """Reject incomplete evidence provenance and malformed triage rules."""
    _require(config.get("schema_version") == 1, "Unsupported Track 2 schema_version")
    sources = config.get("sources", {})
    _require(isinstance(sources, dict) and bool(sources), "sources must be a nonempty object")
    for source_id, source in sources.items():
        _require(bool(source.get("title")), f"{source_id}: title required")
        _require(source.get("kind") in {"primary_study", "regulatory", "official_requirements"},
                 f"{source_id}: primary, regulatory, or official source required")
        parsed = urlparse(source.get("url", ""))
        _require(parsed.scheme == "https" and bool(parsed.netloc), f"{source_id}: HTTPS URL required")
        _require(bool(source.get("accessed_on")), f"{source_id}: access date required")
    mechanisms = config.get("mechanisms", [])
    _require(isinstance(mechanisms, list) and bool(mechanisms), "mechanisms must be a nonempty list")
    for mechanism in mechanisms:
        _references(mechanism.get("source_ids"), sources, mechanism.get("gene", "mechanism"))
    candidates = config.get("candidates", [])
    _require(isinstance(candidates, list), "candidates must be a list")
    seen = set()
    for candidate in candidates:
        candidate_id = candidate.get("id")
        _require(bool(candidate_id) and candidate_id not in seen, "Candidate ids must be unique and nonempty")
        seen.add(candidate_id)
        for field in ("name", "rationale", "uncertainties", "validation", "mechanism_action"):
            _require(bool(candidate.get(field)), f"{candidate_id}: {field} required")
        _require(candidate.get("scope") in {"disease_modifying_hypothesis", "downstream_stress_hypothesis", "symptom_directed", "excluded_control"},
                 f"{candidate_id}: invalid scope")
        _require(candidate.get("direction") in {"compatible", "uncertain", "opposed"},
                 f"{candidate_id}: invalid mechanism direction")
        _references(candidate.get("source_ids"), sources, candidate_id)
        approval = candidate.get("approval", {})
        _require(approval.get("status") in {"approved", "unverified", "not_approved"},
                 f"{candidate_id}: invalid approval status")
        _require(bool(approval.get("jurisdiction")), f"{candidate_id}: approval jurisdiction required")
        if approval["status"] == "approved":
            refs = approval.get("source_ids")
            _references(refs, sources, f"{candidate_id} approval")
            _require(any(sources[ref]["kind"] == "regulatory" for ref in refs),
                     f"{candidate_id}: approved status requires regulatory provenance")
            _require(bool(approval.get("indication")) and bool(approval.get("formulation")),
                     f"{candidate_id}: approved indication and formulation required")
        gates = candidate.get("gates", {})
        for gate in REQUIRED_GATES:
            record = gates.get(gate, {})
            _require(record.get("status") in {"pending", "pass", "fail"}, f"{candidate_id}: invalid or missing gate {gate}")
            _require(bool(record.get("reason")), f"{candidate_id}: gate {gate} needs a reason")
            if record["status"] in {"pass", "fail"}:
                _references(record.get("source_ids"), sources, f"{candidate_id} gate {gate}")
        for experiment in candidate["validation"]:
            for field in ("question", "design", "readouts", "falsification"):
                _require(bool(experiment.get(field)), f"{candidate_id}: validation {field} required")
    disclosure = config.get("ai_disclosure", {})
    _require(isinstance(disclosure, dict), "ai_disclosure must be an object")


def triage_candidate(candidate: dict) -> dict:
    """Hard exclusions precede provenance/functional gates; never rank efficacy."""
    reasons = []
    if candidate.get("exclusion_reason"):
        reasons.append(candidate["exclusion_reason"])
    if candidate["direction"] == "opposed":
        reasons.append("Drug action conflicts with the proposed rescue direction.")
    failed = [name for name, gate in candidate["gates"].items() if gate["status"] == "fail"]
    reasons.extend(f"Failed evidence gate: {gate}." for gate in failed)
    if reasons:
        state = "rejected"
    else:
        if candidate["approval"]["status"] != "approved":
            reasons.append("Approved medication status is not established by this ledger.")
        if candidate["direction"] != "compatible":
            reasons.append("Mechanistic direction remains uncertain in the relevant model.")
        pending = [gate for gate in REQUIRED_GATES if candidate["gates"][gate]["status"] != "pass"]
        reasons.extend(f"Unresolved evidence gate: {gate}." for gate in pending)
        state = "deferred" if reasons else "accepted"
        if not reasons:
            reasons.append("All experimental-prioritization gates passed; this is not evidence of clinical benefit.")
    return {**candidate, "state": state, "decision_reasons": reasons}


def _link_sources(refs: list, sources: dict) -> str:
    return "; ".join(f"[{sources[ref]['title']}]({sources[ref]['url']})" for ref in refs)


def _cell(value: str) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def render_report(result: dict, config: dict) -> str:
    sources = config["sources"]
    counts = result["counts"]
    lines = [
        "# Track 2: mechanism-to-drug research draft", "", DISCLAIMER, "",
        f"Evidence review: {config['reviewed_on']}. Config SHA-256: `{result['config_sha256']}`.", "",
        f"Triage: {counts['accepted']} accepted for further experiments; {counts['deferred']} deferred; "
        f"{counts['rejected']} rejected. These are categorical review decisions, not efficacy scores.", "",
        "## Case linkage and scope", "", config["case_linkage"], "",
        "No patient-level data were used by this module. Gene-level literature is conditional background. "
        "A high variant rank does not establish loss or gain of function. The causal allele, transcript, "
        "inheritance, protein abundance, and cellular phenotype must be checked before proposing a mechanism.", "",
        "## Conditional disease mechanisms", "",
    ]
    for mechanism in config["mechanisms"]:
        lines += [f"### {mechanism['gene']}", "", mechanism["description"], "",
                  f"Required direction: {mechanism['required_direction']}", "",
                  _link_sources(mechanism["source_ids"], sources), ""]
    lines += ["## Candidate decision table", "",
              "| Candidate | State | Role | Direction | Approval evidence |",
              "|---|---|---|---|---|"]
    for candidate in result["candidates"]:
        approval = candidate["approval"]
        lines.append("| " + " | ".join(_cell(value) for value in (
            candidate["name"], candidate["state"], candidate["scope"], candidate["direction"],
            f"{approval['jurisdiction']}: {approval['status']}")) + " |")
    lines += ["", "Symptom-directed hypotheses, if added, must remain distinct from disease modification. "
              "Changing a stress marker alone cannot establish corrected chromosome segregation.", ""]
    for candidate in result["candidates"]:
        approval = candidate["approval"]
        lines += [f"## {candidate['name']} — {candidate['state']}", "",
                  f"**Action:** {candidate['mechanism_action']}", "", candidate["rationale"], "",
                  f"**Regulatory scope:** {approval['jurisdiction']}; {approval['status']}. "
                  f"{approval.get('indication', 'No approved indication verified')}. "
                  f"Formulation: {approval.get('formulation', 'not established')}. "
                  "An approval for another indication does not establish MVA approval or suitability.", ""]
        for reason in candidate["decision_reasons"]:
            lines.append(f"- {reason}")
        lines += ["", "**Evidence limitations and possible harm**", ""]
        lines.extend(f"- {item}" for item in candidate["uncertainties"])
        lines += ["", "**Falsifiable validation plan**", ""]
        for experiment in candidate["validation"]:
            lines += [f"- Question: {experiment['question']}",
                      f"  Design: {experiment['design']}",
                      f"  Readouts: {experiment['readouts']}",
                      f"  Falsification: {experiment['falsification']}"]
        lines += ["", "**Evidence gates**", ""]
        for gate in REQUIRED_GATES:
            record = candidate["gates"][gate]
            lines.append(f"- {gate}: {record['status']} — {record['reason']}")
        refs = list(dict.fromkeys(candidate["source_ids"] + approval.get("source_ids", [])))
        lines += ["", _link_sources(refs, sources), ""]
    lines += ["## Methods and interpretation", "",
              "The local standard-library pipeline validates source provenance, records approval jurisdiction "
              "and formulation, enforces direction-of-effect exclusions, and requires five documented gates "
              "before experimental acceptance. It never fetches data, infers efficacy from citation counts, "
              "or converts a genetic rank into a treatment recommendation. The initial ledger is a focused "
              "literature seed, not an exhaustive systematic review. Selection bias is possible.", "",
              "Experiments should use blinded analysis, independent biological replicates, matched vehicle "
              "and isogenic correction controls, and prespecified endpoints. Reduced aneuploid fraction can "
              "reflect selective killing or cytostasis; live-cell lineage tracking and viable-cell counts "
              "must distinguish correction from depletion. Do not equate a stronger checkpoint signal, "
              "lower cytokines, or improved viability with restored organismal development.", "",
              "## Submission readiness", ""]
    lines.extend(f"- {requirement}" for requirement in config["submission_requirements"])
    lines += ["", _link_sources(config["requirement_source_ids"], sources), "", "**Outstanding items**", ""]
    lines.extend(f"- {blocker}" for blocker in result["submission_blockers"])
    lines += ["", "**AI methods disclosure**", ""]
    for field in ("provider", "plan_or_tier", "data_handling_setting"):
        lines.append(f"- {field}: {config['ai_disclosure'].get(field) or 'UNCONFIRMED — fill from actual account/settings'}")
    lines += ["", "## Source register", ""]
    for source_id, source in sources.items():
        lines += [f"- **{source_id}**: [{source['title']}]({source['url']}); "
                  f"{source['kind']}; accessed {source['accessed_on']}. {source['supports']}"]
    return "\n".join(lines) + "\n"


def run(config_path: Path, output_dir: Path) -> dict:
    """Write reproducible local JSON, CSV and Markdown evidence-review artifacts."""
    config_path, output_dir = Path(config_path), Path(output_dir)
    raw = config_path.read_bytes()
    config = json.loads(raw.decode("utf-8-sig"))
    validate_config(config)
    candidates = [triage_candidate(candidate) for candidate in config["candidates"]]
    counts = dict.fromkeys(STATES, 0)
    counts.update(Counter(candidate["state"] for candidate in candidates))
    blockers = [
        "Link the locally validated causal allele(s) to the mechanism and document alternatives.",
        "Complete missing functional, genome-safety, and exposure evidence before final candidate claims.",
        "Review the report with a domain expert and complete the official methods template.",
        "Set the team/username report filename, repository URL, and three-minute pitch URL before submission.",
    ]
    for field in ("provider", "plan_or_tier", "data_handling_setting"):
        if not config.get("ai_disclosure", {}).get(field):
            blockers.append(f"Confirm AI disclosure: {field}.")
    result = {
        "schema_version": 1,
        "reviewed_on": config["reviewed_on"],
        "config_sha256": hashlib.sha256(raw).hexdigest(),
        "purpose": "experimental_hypothesis_triage",
        "disclaimer": DISCLAIMER,
        "case_linkage": config["case_linkage"],
        "submission_ready": False,
        "counts": counts,
        "candidates": candidates,
        "sources": config["sources"],
        "submission_blockers": blockers,
        "outputs": {"json": "track2_candidates.json", "csv": "track2_candidates.csv", "report": "track2_research_draft.md"},
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / result["outputs"]["json"]).write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with (output_dir / result["outputs"]["csv"]).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["id", "name", "state", "scope", "direction", "approval_jurisdiction", "approval_status", "decision_reasons", "source_urls"])
        writer.writeheader()
        for candidate in candidates:
            refs = list(dict.fromkeys(candidate["source_ids"] + candidate["approval"].get("source_ids", [])))
            writer.writerow({
                **{key: candidate[key] for key in ("id", "name", "state", "scope", "direction")},
                "approval_jurisdiction": candidate["approval"]["jurisdiction"],
                "approval_status": candidate["approval"]["status"],
                "decision_reasons": " | ".join(candidate["decision_reasons"]),
                "source_urls": " | ".join(config["sources"][ref]["url"] for ref in refs),
            })
    (output_dir / result["outputs"]["report"]).write_text(render_report(result, config), encoding="utf-8")
    return result
