# MVA Hackathon research workflow

This project implements and compares two approaches to both hackathon tracks, a broader drug-discovery screen, and a target-first simulation prototype for Track 2. **Strategy A** follows the research structure in **DJ-1 Modulated Genes in Parkinson's Disease**: collect evidence, build context-specific gene sets, test regulatory enrichment, interpret mechanisms, and propose experiments. **Strategy B** tests a different starting point: genome-wide coding-variant evidence and allele-specific functional rescue. **Strategy C** searches large public drug-target and perturbation-signature libraries under the conditional BUB1B mechanism. Core analyses run locally with Python 3.11 or newer and have no third-party Python dependencies; optional scientific figures use Matplotlib and NumPy.

**The active Track 2 approach is the [target-first framework](outputs/target_first_2026-10-20_v2/FRAMEWORK.md), adopted October 20.** Select intervention targets from the variant-to-phenotype model, test competing mechanisms, then identify drugs capable of the desired action. Healthy-tissue rescue and tumor control remain separate. The initial protein-turnover and attachment-balance prototype has run, but its parameters and phenotype outcomes are not fitted; it is not a drug efficacy ranking.

The [October 19 four-drug comparison](outputs/main_comparison_2026-10-19/COMPARISON.md) remains historical supporting work. It selected sodium 4-phenylbutyrate and niacin as leads, with pioglitazone and gentamicin as comparators. Those drugs have no reserved lead positions in the current framework.

The earlier [October 20 BUB1B-loss comparison](outputs/bub1b_loss_comparison_2026-10-20/ASSESSMENT.md) evaluates nine public perturbation profiles in seven cell models, before and after removing cell-cycle genes and benchmarking against other perturbations. It finds a suggestive unfolded-protein-response pattern in cancer cells, with the opposite direction in the separate hiPSC model; no pathway comparison passes correction for multiple testing. Its provisional phenylbutyrate assessment and [separate niacin review](outputs/bub1b_loss_comparison_2026-10-20/NIACIN_ASSESSMENT.md) remain evidence for comparison, not a constraint on target discovery.

## Run the target-first prototype

The target hypotheses, phenotype endpoints, separate objective arms and explicit parameter ranges are in `config/target_first_framework.json`. With the existing public catalog and local allele/phenotype context available:

```powershell
python -m mva.target_first --output outputs\target_first_replay
python -m unittest tests.test_target_first -v
```

The runner evaluates virtual molecular interventions across four allele scenarios and hypothetical tumor states, then audits the complete cached phase-four catalog for target/action matches. Neither a virtual-intervention result nor a catalog match is a drug nomination. Exact gene matches do not recover all indirect or complex-mediated pharmacology. The saved [model specification](outputs/target_first_2026-10-20_v2/MODEL.md) explains the equations, units and omitted biology.

To model a measured functional drug response, copy and complete `config/target_first_drug_profiles.template.json`, then supply it with `--drug-profiles`. Missing potency, maximum effect, functional direction or free-exposure inputs remain unresolved. Binding affinity alone is not accepted as a functional response. Each current profile models one specified action; broader drug effects and phenotype outcomes require additional calibration.

```powershell
python -m mva.target_first --drug-profiles path\to\completed_profiles.json --output outputs\target_first_drug_replay
```

The completed pilot is in `outputs/target_first_2026-10-20_v2/`; the preceding local prototype run is retained in `outputs/target_first_2026-10-20/`. The active objective is recorded in `config/contest_objective.json`, with its prior version preserved beside the new framework. Replays require a fresh output directory. The analysis is local and does not transmit phenotype or variant data.

## Reproduce the BUB1B loss comparison

With the Strategy C public references cached, run the fixed module comparison into a fresh directory:

```powershell
python -m mva.bub1b_loss --output outputs\bub1b_loss_replay
```

The rules are in `config/bub1b_loss_comparison.json`. Outputs include exact gene hits, per-cell scores, complete matched-target reference scores, cell-cycle-mask coverage, XPR-identifier sensitivity comparisons, multiple-testing correction and input/output hashes. LINCS and CM4AI are kept separate; repeated contexts do not count as independent cell models. Gene masking and matched-target comparisons only partially address proliferation and cell-state confounding. These are analyses of public thresholded lists, not raw-expression regression or measurements of pathway activity.

To include PNG, SVG and PDF figures, install the optional local plotting dependencies, then add `--plot` to the analysis command:

```powershell
python -m pip install --target .analysis_dependencies -r requirements-analysis-plot.txt
python -m mva.bub1b_loss --output outputs\bub1b_loss_replay_with_plot --plot
```

The saved assessment and niacin review interpret the completed run; they are not regenerated automatically by the analysis command. Replays require a new evidence review if inputs or results change.

## Run Strategy A: DJ1-style baseline

From this project directory:

```powershell
python -m mva prepare
python -m mva run
python -m unittest discover -s tests -v
```

`prepare` downloads public Ensembl canonical transcripts and Enrichr libraries with source receipts and checksums. It never uploads patient genes, variants or files. Once references are cached, `run` needs no network access and creates a new dated directory under `outputs/`. Open its `REVIEW.md` for both tracks, limitations and next experiments. Use `python -m mva --help` for individual stages.

The saved comparison baseline is `outputs/review_ready/`. Its manifest and output hashes are checked before and after the Strategy B run.

## What Strategy A does

- **Intake:** inventory the supplied files, extract the HPO table, preserve family history separately, fingerprint the VCF, and check the first 100 read pairs per lane. This is not full read-level QC.
- **Track 1:** audit every VCF record, map called alternate alleles to a literature-selected GRCh38 starter panel, annotate supported coding consequences locally, apply quality checks and retain auditable recessive hypotheses. Genotype frequency is never mistaken for population frequency. Phase and unknown annotations remain explicit.
- **DJ1 enrichment:** run hypergeometric overrepresentation with Benjamini–Hochberg FDR using local GMT files and explicit backgrounds. Human ChEA experiments and motif predictions remain distinct. CCDC84 is explicitly mapped to CENATAC using documented primary evidence.
- **Track 2:** evaluate public drug evidence by mechanism direction, approved indication/formulation, uncertainty, counter-evidence and experimental gates. The initial literature set does not imply patient expression or drug efficacy.
- **Review:** connect the molecular hypotheses to the mechanism ledger and generate local draft submission rows. EPCR is deliberately blank pending scientific review; heuristic scores are not probabilities.

## Run Strategy B: genome-wide evidence and allele-specific rescue

Prepare the public reference cache and the local phenotype input, then choose a fresh output directory:

```powershell
python -m mva.alternative_references download
python -m mva.alternative_references models
python -m mva.alternative_references phenotypes
python -m mva.strategy_b --baseline outputs\review_ready --output outputs\strategy_b_replay
```

`download` acquires the full public Ensembl 116 GTF/CDS references, HPO release `v2026-09-01`, and the GRCh38 ClinVar VCF under `references/alternative/`, with download receipts and SHA-256 checksums. `models` builds `references/alternative/genes_grch38.jsonl.gz`; rejected transcript models and their reasons are recorded in `gene_models_summary.json`. The reference commands accept `--directory`, but the combined Strategy B runner currently uses the default `references/alternative/` location.

`phenotypes` reads the existing `outputs/review_ready/intake/phenotype.json` and combines its features with the public HPO associations in `outputs/strategy_b_inputs/phenotype_reference.json`. This combined file contains patient-derived information and stays under local outputs. Family history is retained separately. Neither this preparation nor the analysis sends patient queries or files to external services.

The initial alternative output directory is `outputs/strategy_b/`. The replay command above uses a new directory because the runner refuses to overwrite an existing result. It verifies that the VCF is identical to the baseline input and preserves the baseline artifacts. Once the reference files and local phenotype input exist, the Strategy B runner works offline.

## What Strategy B does

- **Track 1 scope:** scan the same complete VCF using an exon/splice interval index across eligible canonical nuclear protein-coding transcripts. Reference selection is independent of disease-gene names. Nonstandard, incomplete or ambiguous transcript models are excluded with recorded reasons.
- **Evidence:** retain exact local ClinVar allele matches, review status, conflicting classifications and condition labels. Available ExAC/ESP/1000 Genomes frequency fields are source-labeled; absent frequencies remain unknown. Decompose each research score into genomic, genotype-pattern, patient phenotype, family-history and clinical evidence components.
- **Inheritance and robustness:** show biallelic hypotheses separately from exploratory single alleles. Use catalogued gene inheritance as context, retain unknown phase, and exclude same-block supplied cis pairs. Rebuild rankings after removing phenotype, all ClinVar-derived evidence, or one patient HPO feature at a time. Full ClinVar removal also removes its frequency evidence and rebuilds eligibility.
- **Track 2:** conditionally investigate rescue of the baseline-nominated mechanism through stop-codon context, predicted nonsense-mediated decay, readthrough-product identity and restored checkpoint function. Compare medicinal probes using current approval evidence, conflicting results, exposure gaps, normal-versus-tumor effects and explicit failure criteria. The mechanism branch is disclosed in `config/strategy_b.json`; it does not constrain the genome-wide gene screen.
- **Comparison:** keep both approaches available for review. Assess coverage, evidence quality, ranking sensitivity and experiments that distinguish plausible mechanisms; no contest accuracy or winning probability is inferred.

Open `outputs/strategy_b/COMPARISON.md` for the combined review. Each Strategy B output directory also contains:

| File | Purpose |
| --- | --- |
| `track1/genes.csv` | Gene priorities with component scores, inheritance context and the best available biallelic hypothesis |
| `track1/biallelic_genes.csv` | Independently ranked homozygous and possible compound-heterozygous hypotheses |
| `track1/candidates.csv` and `track1/candidates.json.gz` | Auditable variant consequences, quality flags, phase and source evidence |
| `track1/sensitivity.csv` | Evidence-removal ranks, scores and explicit score ties |
| `track1/qc.json` | Scan coverage, reference exclusions, scoring policy and evidence limitations |
| `track2/RESCUE_REVIEW.md` | Conditional rescue rationale, source citations, counter-evidence and experimental gates |
| `comparison.json` and `run_manifest.json` | Machine-readable comparison and reproducibility fingerprints |

## Run Strategy C: broad computational drug screening

The historical selection is recorded in `config/main_comparison.json`, the [four-drug comparison](outputs/main_comparison_2026-10-19/COMPARISON.md), and its [structured evidence](outputs/main_comparison_2026-10-19/decision.json). The [earlier eight-nomination review](outputs/contest_hypotheses/REVIEW.md) and [identity/context audit](outputs/contest_hypotheses/screen_evidence.json) remain available as supporting research. The current objective in `config/contest_objective.json` selects targets before drugs and evaluates theoretical rescue and tumor effects separately. Experimental rescue in the patient and direct allele binding are not nomination prerequisites.

The saved Strategy C results are in [the combined review](outputs/drug_screen_c/REVIEW.md) and [the evidence matrix](outputs/drug_screen_c/evidence_matrix.csv). To download the public references and reproduce the analysis in a fresh directory:

```powershell
python -m mva.drug_catalog download
python -m mva.drug_pathways --seed BUB1B
python -m mva.drug_screen --output outputs\drug_screen_c_replay\targets
python -m mva.drug_signatures --download --output outputs\drug_screen_c_replay\signatures
python -m mva.drug_screen_review --root outputs\drug_screen_c_replay
```

`drug_catalog download` acquires the complete public ChEMBL phase-four molecule collection and supporting mechanism, human-target and warning collections; it also builds the normalized catalog. `drug_pathways` downloads the complete Reactome v97 GMT, names, hierarchy and descriptions, then selects the requested gene's modules locally. `drug_signatures --download` acquires complete public drug, knockout and knockdown libraries before local selection and comparison. References, source receipts and file hashes are retained under `references/drug_screen/`. Patient files, mutation coordinates and gene-specific queries are not sent to these services.

With references already cached, replace `drug_catalog download` with `drug_catalog build`, run `python -m mva.drug_pathways --offline --seed BUB1B`, and omit `--download` from the signature command. The analysis then works offline. Choose another fresh output root for each replay. The target runner defaults to `outputs/drug_screen_c/targets`; the signature runner requires `--output` to run the analysis; the combined-review runner defaults to `outputs/drug_screen_c` and requires both completed branches beneath it. It also verifies the saved Strategy A and B outputs, so those local results must remain available.

## What Strategy C does

- **Catalog coverage:** retain 4,225 ChEMBL phase-four molecule records, consolidated into 2,996 parent families for target testing. Phase four is a database development-stage annotation; it does not verify current human approval, jurisdiction, formulation or market availability. Withdrawals and other source warnings remain visible.
- **Pathway relevance:** select eight specific BUB1B-containing Reactome modules, comprising 216 unique genes. The fixed 5–200-gene size limit applies to each module, not their union. Broad categories and Disease descendants are excluded; eligible ancestors are replaced by their closest eligible descendants. These are undirected co-membership sets, not demonstrated downstream causal relationships.
- **Target tests:** test unambiguous direct human single-protein annotations, keeping complexes and protein families as contextual evidence. Correct all assessable parent-family tests using Benjamini–Hochberg, with both an ordinary overlap test and a separate test conditioned on target-annotation popularity. In the saved run, 1,312 families have eligible target data, only two module genes are covered, and no family passes the corrected threshold. Missing coverage is not evidence of lack of efficacy.
- **Expression proxies:** compare 78,463 paired drug up/down contexts representing 5,888 compound names against public BUB1B knockout and knockdown signatures. Preserve cell, dose, duration, source identifier and disagreement across contexts. Same-cell comparisons and a separate hiPSC cross-cell comparison remain distinct. These are thresholded gene-list comparisons, not full expression measurements or this patient's mutation-specific response. Unlisted genes are assigned zero for list membership; that does not establish unchanged expression.
- **Interpretation:** signature reversal scores are descriptive and have no calibrated significance or multiple-testing correction. Cytostasis, cell death and suppression of protective compensation can look like reversal. Exact drug-name matches do not establish structural identity or current approval. The combined matrix preserves separate evidence layers without inventing a composite efficacy score. No molecular docking, treatment simulation or functional rescue experiment has been performed.

The screen produces hypotheses for subsequent mechanism, exposure and functional review. None of its rankings establishes a drug's benefit or safety for this case. Its saved historical review records a known adverse mechanistic comparator that reversal alone failed to reject. Later reviews advanced theoretical drug nominations. The current target-first framework tests the intervention mechanism before promoting a drug; all prior rankings remain available as supporting evidence.

## Research boundaries

All three strategies are exploratory research, not completed diagnoses or treatment recommendations. Strategy A's seven-gene panel is documented in `config/genes.json` and is not exhaustive. Strategy B broadens the coding screen and adds explicit phenotype/clinical evidence, but does not exhaustively assess alternative transcripts, regulatory variants, CNV/structural events, repeat expansions, mitochondrial mechanisms or mosaicism. Independent consequence annotation, read-based confirmation and segregation remain necessary. No expression experiment was supplied, so Strategy A's five-gene literature signature is undirected and Strategy C uses public perturbation proxies.

ClinVar's legacy population-frequency fields cover only a subset of alleles and do not replace a comprehensive current population reference. A missing exact ClinVar match is not evidence of rarity, novelty or benignity; equivalent indels are not normalized for matching. HPO membership uses the source's propagated annotations, and correlated features are not independent likelihood evidence. Gene inheritance is a union across associated disease records, not confirmation of the causal condition. Unphased pairs may be in cis, and supplied trans phase still needs validation. A single allele can represent an incomplete recessive or carrier finding.

Computational NMD predictions and readthrough literature do not establish rescue in this sample. More protein alone is insufficient: allele-specific RNA, protein identity, chromosome-segregation function, exposure and normal-versus-tumor effects require experiments. No drug efficacy, clinical benefit or improved contest score has been demonstrated.

The initial requirements search surfaced other teams' reported findings, and the Strategy A result was visible during Strategy B design. Consequently this work must not be described as blinded independent recovery. The analysis code contains no fixed answer coordinates.

## Project layout

```text
mva/          Analysis, enrichment, reference preparation and report generation
config/       Public evidence, gene panel, aliases and literature signature
tests/        Synthetic correctness tests without challenge data
data/         Existing restricted challenge files, excluded from source control
references/   Cached public reference downloads and source receipts
outputs/      Local derived results and run manifests, excluded from source control
```

Complete Strategy A and B runs record code/configuration hashes, the VCF hash, reference hashes and output hashes in `run_manifest.json`. Strategy B additionally records the scoring policy, stage-specific code hashes, baseline fingerprints and whether the baseline stayed unchanged. Strategy C records public-reference and code fingerprints in its branch manifests and a combined manifest; it does not rescan the VCF. Retain the original code, configuration, local inputs, reference bytes and receipts for exact replay; rerun into a fresh output directory. Ensembl 116, HPO `v2026-09-01` and Reactome v97 are pinned. ClinVar, ChEMBL service collections and public signature download URLs can change over time, so their cached bytes, release metadata and receipts are required to reproduce the same evidence.

## Submission preparation

Official Track 1 requires a maximum of ten GRCh38 prediction rows, reviewed EPCR estimates, a methods report and a GitHub repository. Track 2 requires a mechanism-grounded report, repository and three-minute pitch. The official pages also request disclosure of the AI provider, actual account plan/tier and data-handling settings. Unknown account settings are not invented in these drafts.

```powershell
python -m mva validate-submission outputs\YOUR_RUN\submission_review\track1_UNREVIEWED.csv
```

The unreviewed file is expected to fail validation until EPCR is completed. Validation checks format only. Nothing here submits, publishes, emails or uploads results. Review derived data before publishing a source repository.

The official rules require deletion of challenge data and intermediate/derived datasets within 30 days of hackathon close and a deletion-confirmation email. No deletion or email is performed by this project. Keep research outputs local through review.

Requirements were checked against the [official challenge](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026). The submission deadline shown there is October 24, 2026 at 23:59 UTC. Source instructions and current rules are cached under `references/challenge/` for inspection.

## Independent population frequency audit

The October 20 follow-up adds gnomAD v4 references for 13 selected genes and matches 508 local candidate annotations. It separates population evidence from ClinVar classifications. Of these annotations, 332 have usable frequency matches and 176 remain unknown. This does not replace the genome-wide population annotation gap outside the selected genes.

```powershell
python -m mva.gnomad_review download
python -m mva.gnomad_review run --output outputs/gnomad_frequency_review_replay
```

The download sends the public gene IDs in `config/gnomad_frequency_review.json`; the analysis-informed selection is visible to the service. Patient allele matching runs locally. The fixed query, independent-frequency rules, six comparisons and missing-data checks are included in the source. The local audit is in `outputs/gnomad_frequency_review_2026-10-20`; the revised submission package is `outputs/submission_package_2026-10-20_v2`. The original package is retained. All 167 unit tests passed during this revision.

## Predictor and complete-read follow-up

The subsequent Track 1 revision adds REVEL v1.3, CADD v1.7, AlphaMissense and Ensembl MANE v1.4 SpliceAI evidence from original public datasets, matched locally. The scores remain supporting evidence and do not change the frozen ranking or subjective EPCR. The complete eight-file FASTQ scan examined 535,917,552 pairs and competitively mapped the 2,436 recruited pairs against GRCh38 plus decoys. Duplicate-aware fragment counts and stricter quality/retrieval comparisons support both site calls. Phase and measured protein function remain unresolved.

See [reproduction details](docs/READ_AND_PREDICTOR_METHODS.md). All 186 tests pass in the complete Linux environment, including synthetic recruitment and BAM tests. The local reports and supporting summaries are in `outputs/submission_package_2026-10-20_v3`; prior packages are preserved. Raw/recruited reads and read names are not included in the submission/source archives.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.
