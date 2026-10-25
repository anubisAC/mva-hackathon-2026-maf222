# Reproduction instructions

Run from the repository root with Python 3.11+ (original run: 3.12.10). Core analysis uses no third-party Python dependencies. Optional figures use requirements-analysis-plot.txt. The supplied reports and final submission CSV were curated from these outputs; automatic replay does not recreate their prose or subjective EPCR.

## Data-independent check

```powershell
python -m unittest discover -s tests
python tools/run_target_demo.py --output outputs/public_target_demo
```

The public target demonstration uses the same numerical core, seed, parameters and iteration order as the completed pilot. It intentionally omits the separate clinical-context packaging and medicine-catalog audit. Its decompressed grid was compared byte-for-byte with the saved pilot during package validation.

## Full local workflow (authorized challenge access required)

Obtain the original challenge inputs through the official gated access process and place them in data/. This folder must contain the supplied VCF, phenotype DOCX and read files expected by the intake module. See the module's help/error messages for inventory checks. Do not distribute the inputs or include them in this repository. Use fresh output directories; analysis commands refuse to overwrite results.

The following paths reproduce the expected upstream dependency layout in a new checkout. If outputs/review_ready or outputs/strategy_b already exist, use a clean checkout or preserve them and adapt all dependent paths consistently.

```powershell
python -m mva prepare
python -m mva run --data data --out outputs/review_ready
python -m mva.alternative_references download
python -m mva.alternative_references models
python -m mva.alternative_references phenotypes
python -m mva.strategy_b --baseline outputs/review_ready --output outputs/strategy_b
python -m mva.gnomad_review download
python -m mva.gnomad_review run --output outputs/gnomad_frequency_review_replay
python -m mva.drug_catalog download
python -m mva.drug_pathways --seed BUB1B
python -m mva.drug_screen --output outputs/drug_screen_c/targets
python -m mva.drug_signatures --download --output outputs/drug_screen_c/signatures
python -m mva.drug_screen_review --root outputs/drug_screen_c
python -m mva.bub1b_loss --output outputs/bub1b_loss_replay
python -m mva.target_first --output outputs/target_first_replay
```

The gnomAD step downloads complete public references for the 13 gene IDs in config/gnomad_frequency_review.json. Its fixed query is embedded in the source; downloads require no patient inputs. The analysis-informed gene selection is visible to the external service. Local matching uses Strategy B candidates and phenotype weights. This is a scoped frequency audit with six ranking scenarios, not genome-wide gnomAD annotation. Unknown records remain unknown; public caches are omitted from the source archive.

The target run reads the Strategy B phenotype, candidate and allele-context outputs, plus the complete public ChEMBL catalog. It does not estimate kinetic parameters from those clinical files. Measured pharmacology can be supplied later through config/target_first_drug_profiles.template.json and --drug-profiles, with a new output path. Missing evidence stays unresolved; synthetic controls must remain identified as synthetic.

The baseline generator can produce historical literature reviews with earlier drug priorities. Those outputs are historical analyses, not the final Track 2 proposal. The frozen active model is config/target_first_framework.json; the final report distinguishes completed analysis from conditional proposals and future calibration.

## Predictor and actual-read follow-up

See docs/READ_AND_PREDICTOR_METHODS.md for original-dataset predictor lookup, the compiled allele-blind recruiter, full-stream scanning, competitive genome alignment, duplicate handling and separate quality sensitivities. These steps require pysam and the specified external software; core workflow dependencies above remain unchanged.

## Version and provenance limits

Ensembl 116, HPO v2026-09-01 and Reactome v97 are pinned. The original catalog metadata recorded ChEMBL 37. gnomAD API responses, ClinVar, ChEMBL service collections and public signature URLs may change: exact historical outputs require the original cached bytes and their source receipts, not merely the same command. Manifests fingerprint inputs, references, source and outputs. The private supporting-results bundle preserves the relevant original manifests and a package hash manifest.

The packaged source reflects October 20. Earlier runs may fingerprint earlier configurations; this package does not claim those hashes are identical. The exact target-pilot source snapshot is supplied with the model results. A data-independent virtual replay was checked against the historical grid. The entire private-data workflow was not rerun during packaging.

Runtime evidence: the original Strategy B execution took about 81 seconds with prepared caches; this excludes input generation, downloads and literature work. Hardware-normalized performance and monetary cost were not measured.
