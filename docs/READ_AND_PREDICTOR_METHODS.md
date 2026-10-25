# Site-read and predictor reproduction

These follow-ups add evidence without changing the frozen gene/pair scores or subjective EPCR. The source contains no patient alleles or reads. Obtain authorized inputs separately. All allele matching, recruitment and mapping occur locally. Download scripts contact the listed public resources; whole-gene requests reveal the investigated gene to the provider. The original participant authorized the CADD and Ensembl BUB1B interval requests.

## Environment and public resources

The original execution used Ubuntu 24.04 under WSL, Python 3.12, minimap2 2.31-r1302, pysam 0.24.1 and its bundled samtools 1.24. Install a C compiler and zlib development headers. Use native Linux temporary storage for the approximately 12 GB index; very small writes over a Windows mount can be slow. Allow approximately 24 GB RAM and adequate disk space. The source bundle does not contain software binaries or public reference tables.

```bash
python3 tools/download_read_predictor_resources.py
python3 -m pip install pysam==0.24.1
mkdir -p tools/read_validation_runtime
tar -xjf references/read_validation/minimap2-2.31_x64-linux.tar.bz2 -C tools/read_validation_runtime
gcc -O3 -std=c11 -Itools tools/recruit_reads.c -lz -o tools/read_validation_runtime/recruit_reads
gzip -dk references/read_validation/GRCh38_no_alt_plus_hs38d1.fna.gz
python3 -m unittest discover -s tests
```

The downloader also fetches a CPython 3.12 Linux wheel and verifies its PyPI SHA256. Original runs loaded that unpacked wheel from tools/read_validation_runtime/python. A normal compatible installation is sufficient for replay. Core tests need only the standard library; five recruiter tests skip without the compiled Linux executable and one BAM test skips without pysam. All 186 tests ran and passed in the original complete environment. Tests use synthetic inputs.

Public reference: NCBI GRCh38 no-alt plus hs38d1 analysis set, GCA_000001405.15. Check its compressed MD5 against the downloaded NCBI checksum file. This is the same assembly as the VCF but not its exact maskedGRC-exclusions-v2 reference file. Reference alleles are checked locally before retrieval. Decoy sequences are included and alternative loci are excluded.

## Full FASTQ scan and competitive mapping

Run after the Strategy B workflow has created outputs/strategy_b/context/candidates.json. The supplied four paired lanes follow *_L001_R1_*.fastq.gz through L004 with matching R2 files. Use fresh output directories. Do not run the two steps concurrently: the first must finish writing and copying the index before alignment.

```bash
python3 tools/prepare_site_reads.py --reference references/read_validation/GRCh38_no_alt_plus_hs38d1.fna --candidates outputs/strategy_b/context/candidates.json --data data --output outputs/read_replay --mapper tools/read_validation_runtime/minimap2-2.31_x64-linux/minimap2
python3 tools/align_site_reads.py --reference references/read_validation/GRCh38_no_alt_plus_hs38d1.fna --candidates outputs/strategy_b/context/candidates.json --output outputs/read_replay --mapper tools/read_validation_runtime/minimap2-2.31_x64-linux/minimap2
```

An optional --index path on the alignment step can point to a byte-identical index in faster native Linux storage. No original FASTQ has to be copied there. Policy: 750-base flanks, target bases masked to N; unique canonical 17-mers recruit either mate and retain both; 21-mer membership provides a stricter sensitivity comparison. All ideal 150-base reads spanning each target are checked equally for reference and alternate retrieval. Four complete stream scans validate paired names, sequence/quality lengths, matched record counts and gzip integrity. This is targeted read retrieval followed by whole-reference mapping, not full-genome read alignment or variant recall. Exact-kmer selection may miss highly divergent or error-rich reads.

Mapping uses minimap2 short-read preset, secondary alignments retained, up to ten secondary records. Samtools fixmate/markdup marks coordinate duplicates across lanes treated as one library (no UMI evidence). Primary evidence requires a proper pair, MAPQ >=30, base quality >=20, at least five bases from a read end, no duplicate/secondary/supplementary/QC-fail flags and an actual aligned base. Overlapping mates count once per fragment; discordant qualifying mates are excluded. Counts are repeated at base quality 30 and with 21-mer recruitment. A genuine shared fragment would be required for direct phasing.

Source snapshots preserve the exact initial execution scripts. Its slow index-writing attempt was interrupted after all lane scans were independently completed or continuing; the same public reference and mapper preset were used to rebuild the index in native temporary storage. Completed scans were not repeated. The portable preparation script incorporates that storage improvement.

## Original-dataset predictors

```bash
python3 tools/lookup_predictors.py --local-only
python3 tools/lookup_predictors.py --remote-only
python3 -m mva.predictor_review
```

Local-only scans the complete public chromosome-15 REVEL v1.3 archive and global AlphaMissense hg38 table. Remote-only fetches the entire public BUB1B interval from CADD v1.7 and Ensembl MANE v1.4 SpliceAI raw/masked tables by indexed HTTP ranges. Candidate coordinates and alleles are not used as remote query arguments. Final matching uses exact GRCh38 chromosome, position, reference and alternate. Missense matching additionally checks amino acids and transcript base ID. REVEL and AlphaMissense are not applicable to stop-gained alleles. AlphaMissense transcript version differs from the local model; the version difference is recorded. SpliceAI header version is 1.3.1 and the downloaded dataset uses a 50-base maximum splice-site distance. Identical duplicate source annotations are deduplicated, while their original match counts are retained.

Scores are correlated model evidence, not independent votes, measured residual protein function, a clinical classification or drug response. Public caches, source receipts, exact matching rows and manifests are required for historical byte-level replay. Current remote URLs can change. The local supporting package includes result summaries and receipts; original patient reads, recruited FASTQs, BAMs and read-name tables remain outside the submission/source archives.
