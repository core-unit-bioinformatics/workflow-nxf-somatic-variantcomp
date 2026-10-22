# core-unit-bioinformatics/workflow-nxf-somatic-variantcomp: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Dates follow ISO 8601 (`YYYY-MM-DD`).

## Unreleased

### `Added`

- Narrative report `report/variantcomp_report.html` (`bin/report_html.py`, `REPORT_HTML`): Summary (key numbers, overall metrics, recall vs. precision per sample, truth variants found, discordance classes), Introduction, Results (compared region step by step, input variant sets and their filtering, metrics, truth variants, every discordance class explained with its variants, annotation per technology and set, known variants) and Methods (tools, parameters, coverage, consistency checks); offline, inline plots (decision D33)
- datavzrd table browser `report/tables/` (nf-core `datavzrd`): truth variants, all compared variants, one table per discordance class, removed records, accounting, region steps, coverage, enrichment, genes, required variants, checks; linked by `variant_id`
- Variant accounting per VCF (`variant_accounting.py`, `VARIANT_ACCOUNTING`): records kept, removed (with the reason) and added per preprocessing step, every removed record, and checks (`vcf_preprocessed/accounting/`); MultiQC section
- Region steps (`region_summary.py steps`, `REGION_STEPS`): size of the comparison region after every step (`regions/common_targets.steps.tsv`); MultiQC section
- Coverage in the common region and the target region of every alignment (`coverage_summary.py`, second mosdepth run on the target BED, `--coverage_thresholds`, default `10,20,30,50,100`): `coverage/coverage.coverage.tsv` and MultiQC General Statistics columns
- Report tables `variantcomp.truth_variants.tsv` (every truth variant, found / class per query technology) and `variantcomp.annotation_counts.tsv` (annotation per sample, technology and set: truth variants, query calls, FN, FP)
- Consistency checks at the end of the run (`consistency_checks.py`, `CONSISTENCY_CHECKS`, `report/consistency_checks.tsv`): accounting, metrics vs. variant tables, aggregation, discordance classes, regions, required variants; a failed check stops the run after the reports were written (`CHECKS_GATE`, decision D34)
- Overview metro map in the style of the nf-core/rnaseq map (`scripts/make_metro_map_overview.py`, `docs/images/variantcomp_metro_map_overview.{svg,png}`): seven stages, main tools only, one line per input with its branches; shown in the README above the detailed map
- `scripts/check_syntax.py`: strict-syntax check before commits: `nextflow inspect` with the profiles of `nf-core pipelines download` (minimum and current Nextflow), `nextflow lint` without errors or unaccepted warnings

### `Changed`

- Dev env: Nextflow 26.10.2 (was 26.04.4)
- SEQC2 test data are prepared with the separate **seqc2-toolkit** (dataset configs `variantcomp_small` / `variantcomp_full`, workflow profile `variantcomp`) instead of `tests/seqc2/seqc2_testdata.py`, which is removed (decision D31). The test data layout changed: the samplesheet and known-variant sheet are in `seqc2/small/variantcomp/`, VCFs are named `<pair>_<technology>.<aligner>.<caller>.vcf.gz`, the GTF `gencode.v46.basic.annotation.gtf.gz`; the file contents are unchanged
- GIAB genome stratifications: one file parameter per stratification (`--giab_homopolymer`, `--giab_tandem_repeat`, `--giab_segmental_duplication`, `--giab_low_mappability`, `--giab_all_difficult`); a given file switches its annotation on, no file (the default) means no GIAB annotation. Replaces `--giab_stratifications`, `--giab_dir` and `--giab_files` (decision D30). `download_resources.py` and the SEQC2 dataset script write one parameter per file.

### `Fixed`

- `region_summary.py`: overlapping intervals of a BED are merged (and counted in the log) instead of stopping the run
- Known-variant VCFs: `known_variants.py prepare` copied the ID into `INFO/KV_LABEL`; an ID list like `COSM50383;COSM50383` then broke the INFO column and `bcftools norm` failed (`Could not add dummy header for INFO`). The records are now kept as they are (ID, QUAL, FILTER and all INFO fields with their header definitions; multiallelic records are split by `bcftools norm`); only `INFO/VC_T_VAF` (now Number=A) and, for CSV/TSV lists, the ID and `INFO/KV_GENE` are added. `KV_LABEL` is removed: the label is the ID. The check report counts records (`records_read`, `records_kept` instead of `variants_read`, `variants_kept`)
- The pipeline didn't compile after the GIAB file parameters (stray brace in the utils subworkflow; `nextflow inspect` / `nf-core pipelines download`: "Script compilation failed")

## v0.15.0-dev - [2026-10-14]

Implements the feedback of 2026-10-13 (`docs/plan-feedback-20261013.md`, decisions D23–D27).

### `Added`

- `CHECK_INPUTS`: all inputs are checked against the reference before any other step runs: every VCF (contigs, samples, REF alleles against the FASTA), alignment header, BED (target, regions, exclude, repeats, strata, known-variant regions), the stratification and known-variant sheets and the known-variant lists; toggles `--check_vcf_contigs`, `--check_ref_alleles`, `--check_bam_contigs`, `--check_bed_contigs`, `--check_known_variants`
- `--keep_mito` (default `false`): chrM/MT is dropped from all inputs and reported
- `--min_region_bp` (default 10000): minimum size of the common and the callable regions
- `--callable_regions` and `--read_evidence` (default `true`): use the alignments for callable regions / read evidence; both need BAMs (checked at the start); without them, BAMs are ignored
- `--min_coverage_tumor` / `--min_coverage_normal` (tumor/normal thresholds); `--require_normal_coverage` requires normal BAMs for paired rows
- Early PASS split (`--filter_pass_only`): non-PASS records are normalized in a side branch (`NORMALIZE_VCF`); new discordance class `filtered_non_pass` with the FILTER value
- `--vaf_source auto|af|ad` (default `auto`: FORMAT/AF, else AD) and `--max_missing_vaf` (default 0.05: warning and NA VAF strata instead of an error)
- Parameters for formerly fixed thresholds: `--evidence_min_mapq`, `--evidence_min_baseq`, `--evidence_max_depth`, `--amplicon_keep_duplicates`, `--enrichment_fdr`
- `ANNOTATE_VARIANTS`: Ensembl VEP (nf-core `ensemblvep/vep`, offline cache or `--vep_gtf`) + `bcftools +split-vep` (consequence, gene, impact, variant class, biotype); msisensor-pro repeat regions (nf-core `msisensorpro/scan`, or `--repeats_bed`); new columns in `variants.annotated.tsv`; params `--vep`, `--vep_cache`, `--vep_gtf`, `--vep_cache_version`, `--vep_genome`, `--vep_species`, `--vep_args`, `--vep_cache_args`, `--annotate_repeats`, `--repeats_bed`
- Report: discordance by consequence group, repeat context, variant class and impact (`variantcomp.strata.tsv`, `discordance_by_annotation.png`); recall of the required known variants per technology and over all samples
- MultiQC tables of every step summary (`bin/mqc_tables.py`, `MQC_TABLES`): input checks, `bcftools norm` logs, harmonization, regions, callable regions
- `docs/output.md` lists every output file with its columns and types; `docs/roadmap.md` collects further work
- Synthetic GTF (test data) and a GENCODE chr22 GTF for SEQC2 (`tests/seqc2/`)
- GIAB genome stratifications v3.5 in the annotation: `--giab_stratifications` (default `true`), `--giab_dir`, `--giab_files`; the BEDs are checked in `CHECK_INPUTS` and give one `giab_<name>` column per file
- `scripts/download_resources.py`: downloads the VEP cache (Ensembl CHECKSUMS), the GIAB stratifications (md5) and optionally a GENCODE GTF for offline runs; writes `resources.config` / `resources.params.yaml` with the parameters

### `Changed`

- `check_reference.py` rewritten on pysam (subcommands `vcf`, `bam`, `bed`, `rename-bed`); contig lengths are never fixed: a nuclear length mismatch is an error
- `CHECK_REFERENCE` takes the input file itself; new module `RENAME_BED`
- `known_variants.py prepare` checks REF alleles against the FASTA; no maximum number of required variants
- Python code formatted with black (line length 100), linted with ruff
- `conf/modules.config` split into `conf/modules/<subworkflow>.config`
- `BED_INTERSECT_ALL` sorts its input BEDs
- Benjamini–Hochberg from `scipy.stats.false_discovery_control` (SciPy ≥ 1.11)
- `present_but_filtered` now means a PASS call removed by the uniform filters or the region
- `vcf_preprocessed/reference_check/` renamed to `vcf_preprocessed/norm_logs/`
- Metro map redrawn: a central VCF line; input lines merge into the steps that use them
- Docstrings with INPUT/OUTPUT columns in every script, a docstring for every function

### `Fixed`

- `tests/seqc2/seqc2_testdata.py`: downloads in Python (`urllib`, resumable) with an explicit CA bundle, also exported for samtools/bcftools (`CURL_CA_BUNDLE`); fixes `curl: (77) error adding trust anchors` with conda curl builds; `--ca-bundle`
- SEQC2 small dataset: alignments for all pairs (callable regions and read evidence need them for every row); GIAB stratifications in the release layout for `--giab_dir` instead of `stratification.csv`

- Tests on small machines: tasks ran one at a time (the 3 GB test limit capped every 6 GB label request to the whole machine) and every task spent about 7 s in `conda info` + `source activate`; the test profiles now request small resources per label and nf-test activates the environments with micromamba (10x faster `PREPROCESS_VCF` test on the Pi)
- `mqc_tables.py`: row order of the MultiQC tables no longer depends on the order in which the inputs arrive

### `Removed`

- `--gene_bed` and `BEDTOOLS_MAP_GENES` (genes come from VEP)
- `--min_coverage` (replaced by `--min_coverage_tumor` / `--min_coverage_normal`)
- Header-extraction steps `BCFTOOLS_VIEW_HEADER`, `BCFTOOLS_QUERY_CONTIGS`, `SAMTOOLS_IDXSTATS` and the `check_ref_alleles` module (replaced by `CHECK_INPUTS`)

## v0.14.0-dev - [2026-10-13]

### `Changed`

- The samplesheet is read as in the nf-core template: `samplesheetToList` (nf-schema) in `PIPELINE_INITIALISATION`; the checks across rows (`check_sheets.py samplesheet`) run as `CHECK_SAMPLESHEET`, the first step of `VARIANTCOMP`, and everything else waits for them
- `--stratification` and `--known_variants` are plain file parameters: `PIPELINE_INITIALISATION` no longer reads them; each sheet is checked by `check_sheets.py` in the subworkflow that uses it (`DISCORDANCE`, `PREPARE_KNOWN_VARIANTS`)
- Metro map and README step table: the sheet checks are shown where they run (1.1, 1.2, 7.4, 8.1)

### `Fixed`

- Required known-variant tables list the technologies in name order (the order depended on the order the inputs arrived)
- `discordance/variants.annotated.tsv` and `coverage/callable_summary.tsv` are concatenated in file-name order (the row order depended on the work directory paths)

## v0.13.1-dev - [2026-10-08]

### `Changed`

- Metro map redrawn: one box per subworkflow, one numbered row per step with its tools, one line per use case (always; target BEDs; alignments; genes/stratifications; required known variants; external truth set)
- README: step table with a detailed explanation of every step, generated from the same data as the map (`scripts/make_metro_map.py --markdown`)

## v0.13.0-dev - [2026-10-08]

### `Changed`

- Tools run in nf-core modules, logic in `bin/` scripts (one tool per module, no tool pipes in local modules, all options in `conf/modules.config`):
  - reference check: facts from `bcftools view -h` / `bcftools query` (VCF), the BED itself and `samtools idxstats` (BAM); decisions in `check_reference.py contigs`; REF alleles from the `bcftools norm --check-ref x` log (patched nf-core module emits the log) and `check_reference.py ref-alleles`
  - target regions: `bedtools multiinter` + `gawk` + `bedtools merge` (subworkflow `BED_INTERSECT_ALL`), `bedtools intersect` / `subtract`; summary by `region_summary.py`
  - callable regions: nf-core `mosdepth`, `mosdepth_callable.py`, `bedtools sort`, `BED_INTERSECT_ALL`, `bedtools intersect`
  - read evidence: `bcftools query` (sites), nf-core `bcftools mpileup` (+ `call --keep-alts`), `bcftools norm`
  - annotation: `annotate_variants.py to-bed`, `bedtools sort`, `bedtools map` (genes), `bedtools intersect -c` (strata), `annotate_variants.py merge`
  - exact and position matches: nf-core `bcftools isec` (`--collapse none` / `all`)
  - known variants: the required-variant summary and MultiQC table are written by `known_summary.py` instead of Groovy code; `known_variants.py prepare` no longer reads the FASTA (REF alleles are checked by `bcftools norm`)
  - removed local modules `read_evidence`, `annotate_variants`, `common_regions`, `callable_regions`
- `harmonize_vcf.py`: documented field selection (AD, else AF + DP), new INFO tags `VC_T_SRC` / `VC_N_SRC`, a per-VCF summary (`*.harmonize_summary.tsv`); errors instead of silent skipping for unsplit records, invalid AD/AF/DP and tumor records without VAF (`--max-missing-vaf`, default 0); records without reads (AD 0,0) get an undefined VAF (`zero_depth`, counted) and `FORMAT/AF = nan` counts as missing (both occur in the SEQC2 MuTect2 3.8 VCFs)
- Stricter scripts with inline comments and `ERROR:`/`WARNING:`/`INFO:` messages naming file, record and expected value: `classify_discordance.py` (comparison identity, status values, harmonized and pileup VCF headers), `aggregate_metrics.py` (duplicate rows, value ranges), `known_variants.py` (TSV values, VAF tag, technologies of `--evidence`), `variantcomp_report.py` (columns, status and class values), `check_reference.py` (mixed naming styles, inputs without any reference contig)
- `known_variants/required_variants.tsv` and `required_summary.tsv` start with the list name (`list` column)
- Stratification columns of `variants.annotated.tsv` are in name order

### `Fixed`

- Position-only matches (`match` / `position`) are counted per position: an extra allele at a shared position is no longer an FP
- A single target BED (or none: WGS only) no longer fails the region intersection (`bedtools multiinter` needs two files; a single BED is only merged)
- `bedtools intersect` of the callable and `--regions_bed` regions uses `-sorted` with the reference index, so the output stays in reference order

### `Added`

- Unit tests for `mosdepth_callable.py`, `region_summary.py`, `annotate_variants.py`, `known_summary.py` and the new checks (98 pytest tests)
- nf-core-style metro map (`docs/images/variantcomp_metro_map.svg/png`, `scripts/make_metro_map.py`) and per-step tool table in the README
- Tool and dataset citations in `CITATIONS.md` and the MultiQC methods description
- Documentation of the VCF requirements, the tested inputs and the harmonization rules (`docs/usage.md`), tools per output folder (`docs/output.md`)

## v0.12.0-dev - [2026-10-08]

### `Changed`

- **Breaking:** samplesheet, stratification sheet and known-variant sheet are CSV (comma-separated); known-variant lists may be CSV or TSV; several variant files of a list are separated by `;`
- **Breaking:** the samplesheet has no index columns any more (`vcf_index`, `tumor_bai`, `normal_bai`): index files next to the data files are used (`.tbi`/`.csi`, `.bai`/`.bam.bai`/`.csi`, `.crai`); missing alignment indices are created once with `samtools index` (`PREPARE_ALIGNMENTS`) for coverage and read evidence
- Sheet validation moved from Nextflow code to `bin/check_sheets.py` (`CHECK_SHEETS`): all problems are reported at once with line numbers; Nextflow only stages the files (missing files stop the run with the path)
- Validated sheets are published as `pipeline_info/*.validated.csv`
- ruff: long lines are left to `ruff format` (E501 ignored)

## v0.11.0-dev - [2026-10-08]

### `Added`

- `scripts/make_offline_bundle.py`: pipeline archive at a git tag, list (and optionally pull) of all Singularity images named like Nextflow's cache, Nextflow plugins, README
- `docs/hpc.md`: offline installation, PBS Pro example configuration, full SEQC2 test (T2)

### `Changed`

- README with the implemented workflow steps; plan status and implementation decisions

## v0.10.0-dev - [2026-10-08]

### `Added`

- Known variants (`--known_variants` sheet: name, sample or `all`, mode, variants, regions): `bin/known_variants.py` converts VCF/TSV lists to normalized VCFs on the reference
  - mode `required`: PASS/FAIL per sample, technology and variant (`FAIL_filtered`, `FAIL_outside_region`, `FAIL_not_called`) with read evidence; warnings, MultiQC table and `--fail_on_missing_required`
  - mode `benchmark`: every technology against an external truth set (rtg vcfeval in ALT mode) within its regions and the comparison region, with the main benchmark metrics
- Test data: known-variant lists in the synthetic set; SEQC2 truth set v1.2.1 and a demo required list in the SEQC2 set
- nf-test: shared conda cache for tests (`NFT_CONDA_CACHE`)

### `Changed`

- `rtgtools/vcfeval` module patched: ROC outputs per variant type are optional (not written in ALT mode)
- `benchmark_metrics.py`: `--regions` restricts the match strata; works with VCFs without the `VC_*` tags

## v0.9.0-dev - [2026-10-08]

### `Added`

- `bin/variantcomp_report.py` (`VARIANTCOMP_REPORT`): enrichment of discordant vs. concordant variants (indel, low VAF, stratifications; Fisher's exact test, BH correction), gene summary with per-gene enrichment, discordance class summary, plots, and MultiQC custom-content sections (recall/precision per comparison and per sample, recall by VAF, discordance classes)
- Parameter `--low_vaf`; outputs `report/`
- Container for the report built with Seqera Containers (python, pandas, scipy, matplotlib, upsetplot)

## v0.8.0-dev - [2026-10-07]

### `Added`

- Discordance subworkflow `DISCORDANCE`: `bcftools mpileup` read evidence of every alignment at all compared positions (`READ_EVIDENCE`; duplicates kept for amplicon data), classification of every FN/FP (`present_but_filtered`, `germline_in_normal`, `low_coverage`, `subthreshold_evidence`, `no_evidence`, `absent`) and power-adjusted recall (`bin/classify_discordance.py`), annotation with genes and stratifications (`ANNOTATE_VARIANTS`, bedtools)
- Parameters `--gene_bed`, `--stratification` (sheet with name and BED), `--germline_vaf`, `--evidence_min_alt_reads`, `--power_min_alt_reads`, `--detection_power`
- Outputs `discordance/`, `benchmark/power.*.tsv`
- Test data: gene and stratification BEDs in the synthetic set; GENCODE genes and GIAB v3.5 stratifications in the SEQC2 set

### `Changed`

- Reference check keeps all columns of harmonized BEDs (gene names)

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
