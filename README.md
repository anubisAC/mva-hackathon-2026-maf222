# MVA Hackathon 2026 - maf222

This repository contains local variant triage, public drug-evidence comparisons, and an unfitted target-first simulation prototype. Healthy-tissue rescue and tumor effects are separate objectives. No drug-specific winner, binding prediction, clinical response or tumor-survival estimate is established.

## Quick start without challenge data

Use Python 3.11 or newer (the saved analyses used 3.12.10). Core calculations use only the standard library. From the repository root:

```powershell
python -m unittest discover -s tests
python tools/run_target_demo.py --output outputs/public_target_demo
```

The demonstration produces the same 76,800 virtual molecular evaluations as the saved target pilot, without reading phenotype/variant files, downloading databases or simulating a drug. These are model sensitivity evaluations, including repeated baselines, not independent experiments. All 186 tests passed during packaging, including 12 independent-frequency and 19 predictor/read-evidence checks; the complete Linux environment ran the optional compiled-recruiter and BAM tests. See docs/MODEL.md for equations, assumptions and omitted biology.

## Reproduce the complete analyses

See REPRODUCING.md for the ordered workflow. Authorized participants must obtain the challenge data independently; the source bundle contains no VCF, sequencing reads, clinical documents or private patient-context export. Public reference downloads are also omitted to keep the repository small. Analyses with gated inputs run locally. Do not place those inputs into Git history.

The source snapshot includes historical configurations and evidence ledgers. The active decision is config/contest_objective.json and the active mechanistic policy is config/target_first_framework.json. Earlier fixed drug roles are preserved as history and do not determine the current target model.

Code automates annotation, ranking, enrichment, catalog matching and numerical experiments. Literature appraisal, target selection, medicine proposals and final CSV selection involved AI-assisted curation. The report's EPCR is a subjective contest estimate, not a computed posterior.

## AI use

OpenAI, GPT-6 Astra through Codex; ChatGPT Pro plan; model training disabled, as reported by the participant (maf222) on October 20, 2026. AI assisted programming, literature retrieval, analysis design, interpretation and report drafting. This declaration records the participant's stated settings; it is not an independent audit of provider retention or service terms.

## Source record and data handling

docs/SOURCES.json lists public resources and primary studies. docs/PROJECT_HISTORY.md records the wider project history; its output links refer to locally generated directories, not files shipped in this source bundle. Exact-byte historical replay requires the original permitted reference caches and manifests. Fresh downloads can change results. The organizers' rules govern data retention, publication and submission reuse, including the required post-contest deletion confirmation. No data deletion or email is performed by these tools.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.
