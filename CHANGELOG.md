# core-unit-bioinformatics/workflow-nxf-somatic-variantcomp: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Dates follow ISO 8601 (`YYYY-MM-DD`).

## v0.7.0-dev - [2026-10-07]

### `Added`

- Benchmark subworkflow `BENCHMARK`: `rtg vcfeval` (squash-ploidy, decompose, all records, tumor VAF as score) of every technology against the truth technology per sample; `--comparison_mode all_pairs` adds pairwise comparisons of the other technologies
- `bin/benchmark_metrics.py` (`BENCHMARK_METRICS`): recall, precision and F1 overall, by variant type, by truth/query VAF bin (`--vaf_bins`) and for exact/position-only matching; per-variant TP/FN/FP table
- `bin/aggregate_metrics.py` (`AGGREGATE_METRICS`): micro (summed counts) and macro (mean ± sd over samples) aggregation
- Outputs `benchmark/`
- Test profiles: rtg-tools workaround for aarch64 (local conda tests on a Raspberry Pi)

### `Fixed`

- SEQC2 test-data script: works when the NCBI md5 list is unreachable (cached copy or warning)

## v0.6.0-dev - [2026-10-07]

### `Added`

- Coverage subworkflow `COVERAGE`: reference check of alignment headers (`check_reference.py --type bam`, incl. translation of the common region to the alignment's contig names), indexing of alignments without index, `mosdepth` per alignment (duplicates kept for amplicon libraries), per-sample callable regions (`CALLABLE_REGIONS`; `--callable_mode base|region`, `--min_coverage`, `--require_normal_coverage`)
- The comparison region of a sample is its callable region if alignments are given for all its technologies, else the common region
- `coverage/` outputs and mosdepth results in MultiQC

### `Changed`

- SEQC2 small test data: BAM slices from two 1-Mb windows of chr22 with many exome calls (23-24 Mb, 50-50.8 Mb), subsampled to 50% with a fixed seed

## v0.5.0-dev - [2026-10-07]

### `Added`

- VCF preprocessing subworkflow `PREPROCESS_VCF`: reference check incl. sample names, tumor/normal selection and BGZF recompression, contig renaming, removal of non-reference contigs and gVCF reference blocks, contig header from the `.fai`, `bcftools norm` (split, left-align, atomize, REF check, deduplicate), harmonization, uniform filters and restriction to the comparison region
- `bin/harmonize_vcf.py` (local module `HARMONIZE_VCF`): caller-independent `VC_*` INFO tags from FORMAT/AD, variant type and indel length, tumor GT `0/1`, removal of symbolic alleles
- Filter parameters `--filter_pass_only`, `--min_vaf`, `--min_depth`, `--min_alt_reads`
- Reference FASTA is decompressed once if gzip-compressed; genome BED of all reference contigs
- Outputs `vcf_preprocessed/`

## v0.4.0-dev - [2026-10-07]

### `Added`

- Reference preparation (`--fasta`, `--fasta_fai`; `.fai`/`.gzi` created with `samtools faidx` when needed)
- Reference consistency check `bin/check_reference.py` (local module `CHECK_REFERENCE`): contig renaming (`chr` prefix, `MT`/`chrM`), header length check (nuclear mismatch = build error, mitochondrial = header fix), REF allele check (`--max_ref_mismatch`), dropping of unknown contigs, intervals/records beyond contig ends
- Region preparation: per-technology BED check, sort and merge, N-way intersection of all target BEDs (`COMMON_REGIONS`), optional `--regions_bed` / `--exclude_bed`, region summary
- nf-test module/subworkflow tests and pytest unit tests

## v0.3.0-dev - [2026-10-07]

### `Added`

- `tests/seqc2/seqc2_testdata.py` (+ `seqc2_testdata.yaml`): reproducible download (resumable, md5-verified), preparation and run of the SEQC2 HCC1395/HCC1395BL WES + WGS test dataset in `small` (chr22, Pi-sized) and `full` (genome-wide, HPC) mode, with provenance README and checksums
- Profiles `test_seqc2` (small SEQC2) and `test_full` (full SEQC2)
- nf-test pipeline test for `-profile test_seqc2`; pytest for the SEQC2 script
- `docs/testing.md`

## v0.2.0-dev - [2026-10-07]

### `Added`

- Tab-separated samplesheet (`sample`, `technology`, `library_type`, `tumor_id`, `normal_id`, `caller`, `vcf`, `vcf_index`, `tumor_bam`, `tumor_bai`, `normal_bam`, `normal_bai`, `target_bed`) with `n/a` for empty fields
- Relative samplesheet paths are resolved against the samplesheet's directory
- Cross-row samplesheet checks (duplicates, truth technology per sample, one BED and library type per technology, BAM/index consistency, file existence) reported all at once
- Parameters `--truth_technology` and `--comparison_mode`
- `pipeline_info/samplesheet_validated.tsv` with resolved paths
- Synthetic edge-case test dataset generator `tests/synthetic/make_synthetic_data.py`; `-profile test` uses its output from the test-data repo
- nf-test function tests for the samplesheet checks

## v0.1.0-dev - [2026-10-07]

Initial development version, created with the [nf-core](https://nf-co.re/) template (nf-core/tools 4.1.0).

### `Added`

- nf-core pipeline template (internal name `variantcomp`, repo `workflow-nxf-somatic-variantcomp`), default branch `main`
- Development conda environment `envs/variantcomp-dev.yaml`
- CUBI metadata: `CITATION.cff`, `pyproject.toml` (version + ruff/pytest config)
- Design and roadmap: `docs/pipeline-plan.md`, developer instructions: `CLAUDE.md`

### `Fixed`

### `Dependencies`

### `Deprecated`
