# core-unit-bioinformatics/workflow-nxf-somatic-variantcomp: Output

## Introduction

This document describes every file the pipeline writes, and for every table its columns:
the name, the type of the value and what it means. All paths are relative to `--outdir`.

Conventions used in the column tables:

| Type | Meaning |
| --- | --- |
| contig | contig name of the reference (after renaming, e.g. `chr1`) |
| int | integer; `pos` is 1-based (VCF), BED `start` is 0-based and `end` exclusive |
| float | decimal number (4 digits); `NA` if undefined (e.g. a ratio with denominator 0) |
| yes/no | flag |
| text | free text, a name or a label |

BED files are tab-separated without header; TSV files are tab-separated with one header
line. Missing values are `NA` (numbers) or `.` (names, alleles).

## Pipeline overview

- [Input checks](#input-checks) - reference checks of every input at the start
- [Regions](#regions) - the common comparison region
- [Coverage](#coverage) - depth of every alignment and per-sample callable regions
- [Preprocessed VCFs](#preprocessed-vcfs) - normalized, harmonized and filtered VCFs
- [Benchmark](#benchmark) - recall, precision and F1 of every technology against the truth
- [Discordance](#discordance) - why variants differ, with read evidence and annotation
- [Annotation](#annotation) - VEP consequences and repeat regions
- [Known variants](#known-variants) - required variants and external truth sets
- [Report](#report) - enrichment statistics, gene and annotation tables, plots
- [MultiQC](#multiqc) - the report of all results and step summaries
- [Pipeline information](#pipeline-information) - run reports and checked sheets

Tools that produce each output folder (citations in [`CITATIONS.md`](../CITATIONS.md)):

| Folder | Tools |
| --- | --- |
| `input_checks/` | `check_reference.py` (pysam), `known_variants.py` (pysam) |
| `regions/` | `bedtools sort`, `bedtools multiinter`, `gawk`, `bedtools merge`, `bedtools intersect`, `bedtools subtract`, `region_summary.py` |
| `coverage/` | `check_reference.py`, `mosdepth`, `mosdepth_callable.py`, `bedtools sort`, `bedtools multiinter`, `bedtools intersect`, `region_summary.py` |
| `vcf_preprocessed/` | `bcftools view`, `bcftools annotate`, `bcftools reheader`, `bcftools norm`, `harmonize_vcf.py`, `bcftools concat` |
| `benchmark/` | `rtg vcfeval`, `bcftools isec`, `benchmark_metrics.py`, `aggregate_metrics.py` |
| `discordance/` | `bcftools query`, `bcftools mpileup` + `bcftools call`, `bcftools norm`, `classify_discordance.py`, `annotate_variants.py` |
| `annotation/` | Ensembl VEP, `bcftools +split-vep`, msisensor-pro, `bedtools map`, `bedtools intersect` |
| `known_variants/` | `known_variants.py`, `bcftools norm`, `bcftools isec`, `rtg vcfeval`, `benchmark_metrics.py`, `known_summary.py` |
| `report/` | `variantcomp_report.py` (pandas, SciPy, Matplotlib) |
| `multiqc/` | MultiQC, `mqc_tables.py` |

### Input checks

All inputs are checked against the reference before any other step (subworkflow
`CHECK_INPUTS`); a failed check stops the run. Each check can be switched off
(`--check_vcf_contigs`, `--check_ref_alleles`, `--check_bam_contigs`,
`--check_bed_contigs`, `--check_known_variants`): its problems are then warnings.

<details markdown="1">
<summary>Output files</summary>

- `input_checks/vcf/<sample>_<technology>.vcf.reference_check.tsv`
- `input_checks/alignments/<sample>_<technology>_<tumor|normal>.bam.reference_check.tsv`
- `input_checks/bed/<role>_<name>.bed.reference_check.tsv` (role: `target`, `regions`,
  `exclude`, `repeats`, `stratum`, `known_regions`)
- `input_checks/known_variants/<list>.known.reference_check.tsv`

</details>

Every report has two columns, `key` and `value`. The keys (only those that apply to the
input type are written):

| Key | Type | Meaning |
| --- | --- | --- |
| `status` | text | `ok`, `fixed` (contigs renamed or dropped, warnings) or `error` |
| `input`, `reference_contigs` | text, int | name of the input; number of reference contigs |
| `header_contigs`, `records`, `intervals` | int | contigs in the header (VCF, BAM); records (VCF); intervals (BED) |
| `samples` | text | sample columns of the VCF (comma-separated) |
| `renamed_contigs` | text | renames `input->reference`, e.g. `1->chr1,MT->chrM`, or `none` |
| `unknown_contigs` | text | contigs absent from the reference (their records/intervals are dropped) |
| `dropped_mito` | text / int | the mitochondrial contig (dropped unless `--keep_mito`) |
| `length_mismatches` | text | contigs with another length than in the reference (an error: another build) |
| `records_on_unknown_contigs`, `records_on_mito`, `intervals_dropped` | int | dropped records or intervals |
| `intervals_beyond_contig_end` | int | BED intervals ending after their contig (an error) |
| `ref_alleles_checked`, `ref_mismatches`, `ref_mismatch_rate` | int, int, float | REF alleles compared with the FASTA; mismatches; their fraction (error above `--max_ref_mismatch`) |
| `variants_read`, `variants_kept` | int | known-variant lists: variants read and usable |
| `warning`, `error` | text | one row per message |

### Regions

<details markdown="1">
<summary>Output files</summary>

- `regions/common_targets.bed`: the comparison region, BED3 (contig, start, end): the
  bases targeted by every technology (WGS counts as genome-wide), intersected with
  `--regions_bed` and minus `--exclude_bed`; sorted, merged.
- `regions/common_targets.summary.tsv`: size of every target BED and of the common region.
- `regions/targets/<technology>.bed`: the target BED of each technology after the input
  check (contigs renamed, unknown and mitochondrial contigs removed), sorted; all columns of
  the input BED are kept.

</details>

`common_targets.summary.tsv`:

| Column | Type | Meaning |
| --- | --- | --- |
| `region` | text | technology, or `common` for the common region |
| `type` | text | `target` (technology with a BED), `genome` (no BED, e.g. WGS) or `common` |
| `intervals` | int | number of intervals |
| `bp` | int | size of the region |
| `bp_in_common` | int | bases of the region inside the common region |
| `fraction_in_common` | float | `bp_in_common / bp` |

The run stops if the common region is smaller than `--min_region_bp`.

### Coverage

Only with `--callable_regions` (default; needs alignments for every samplesheet row).

<details markdown="1">
<summary>Output files</summary>

- `coverage/<sample>.callable.bed`: BED3 of the bases of the common region callable in all
  alignments of the sample (tumor depth ≥ `--min_coverage_tumor`, normal depth ≥
  `--min_coverage_normal` with `--require_normal_coverage`). The comparison of this sample
  is restricted to it.
- `coverage/callable_summary.tsv`: one row per sample (columns below).
- `coverage/<sample>/<sample>_<technology>_<tumor|normal>.*`: mosdepth outputs per
  alignment: `*.mosdepth.summary.txt` (mean depth per contig), `*.mosdepth.global.dist.txt`
  and `*.mosdepth.region.dist.txt` (cumulative depth distributions), `*.regions.bed.gz`
  (contig, start, end, mean depth of every common-region interval), `*.quantized.bed.gz`
  (contig, start, end, `LOW_COVERAGE` or `CALLABLE`).

</details>

`callable_summary.tsv`:

| Column | Type | Meaning |
| --- | --- | --- |
| `sample` | text | sample |
| `alignments` | int | alignments (tumor and normal of every technology) intersected |
| `mode` | text | `base` (per-base depth) or `region` (mean depth per interval) |
| `min_coverage` | text | thresholds `tumor/normal` |
| `common_bp` | int | size of the common region |
| `callable_bp` | int | callable bases of the sample |
| `fraction_callable` | float | `callable_bp / common_bp` |

The run stops if a callable region is smaller than `--min_region_bp`.

### Preprocessed VCFs

<details markdown="1">
<summary>Output files</summary>

- `vcf_preprocessed/<sample>/`
  - `<sample>_<technology>.normalized.vcf.gz(.tbi)`: all records (PASS and non-PASS)
    after normalization and harmonization, before the uniform filters and the region
    restriction.
  - `<sample>_<technology>.compare.vcf.gz(.tbi)`: the records that are compared: PASS (with
    `--filter_pass_only`), uniform filters applied, inside the comparison region.
  - `<sample>_<technology>.<pass|non_pass>.normalized.harmonize_summary.tsv` (without
    `--filter_pass_only`: `<sample>_<technology>.normalized.harmonize_summary.tsv`):
    harmonization counts per branch (keys below).
- `vcf_preprocessed/norm_logs/<sample>_<technology>.<branch>.norm.norm.log`: the log of
  `bcftools norm` with the line `Lines total/split/joined/realigned/mismatch_removed/
  dup_removed/skipped`.

</details>

Every VCF goes through the same steps (see CLAUDE.md, "Normalization chain"):

1. Reference check (at the start, see [Input checks](#input-checks)).
2. Tumor and normal columns selected in this order; recompressed with BGZF.
   2b. With `--filter_pass_only` (default): split into PASS records (FILTER `PASS` or `.`)
   and the other records; both go through steps 3–7.
3. Contigs renamed to the reference (`chr` prefix, `MT`/`chrM`).
4. Records on other contigs (and on the mitochondrial contig) and pure gVCF reference blocks
   removed.
5. Contig header lines rewritten from the reference `.fai`.
6. `bcftools norm`: multiallelic records split, indels left-aligned, MNVs atomized into
   SNVs, records with a REF mismatch removed, duplicates removed.
7. Harmonization (`harmonize_vcf.py`, see [usage](usage.md#harmonization-of-vaf-depth-and-alt-reads)).
8. Uniform filters (`--min_vaf`, `--min_depth`, `--min_alt_reads`) and restriction to the
   comparison region (records overlapping it are kept).

INFO tags added by the harmonization (in `*.normalized.vcf.gz` and `*.compare.vcf.gz`):

| Tag | Type | Meaning |
| --- | --- | --- |
| `VC_T_VAF`, `VC_N_VAF` | float | VAF of the tumor / normal (FORMAT/AF or AD[1]/sum(AD), `--vaf_source`) |
| `VC_T_DP`, `VC_N_DP` | int | depth (FORMAT/DP, else sum(AD)) |
| `VC_T_ALT`, `VC_N_ALT` | int | alt reads (AD[1], else round(AF × DP)) |
| `VC_T_SRC`, `VC_N_SRC` | text | source of the VAF: `AF` or `AD` |
| `VC_TYPE` | text | `SNV`, `MNV`, `INS`, `DEL`, `COMPLEX` |
| `VC_LEN` | int | `len(ALT) - len(REF)` |

A header line `##variantcomp_vaf=unavailable` marks a VCF with too many records without
VAF (`--max_missing_vaf`): its VAF strata and detection power are `NA`.

Keys of `*.harmonize_summary.tsv` (`key`, `value`): `records_read`, `records_written`,
`records_dropped_symbolic`, `vaf_available` (yes/no), and per sample (prefix `T_`, `N_`):
`source_AF`, `source_AD`, `source_zero_depth` (AD 0,0: no reads), `source_none` (no
VAF), `depth_DP`, `depth_AD_sum`, `depth_none`, `af_nan` (FORMAT/AF was NaN).

### Benchmark

<details markdown="1">
<summary>Output files</summary>

- `benchmark/metrics.aggregated.tsv`: **main result**: per truth/query technology and
  stratum, over all samples (columns below).
- `benchmark/metrics.per_pair.tsv`: the metrics of every sample and technology pair.
- `benchmark/<sample>/<sample>_<truth>_vs_<query>.metrics.tsv`: the metrics of one pair.
- `benchmark/<sample>/<sample>_<truth>_vs_<query>.variants.tsv`: every compared variant of
  the pair (columns below).
- `benchmark/<sample>/vcfeval/`: rtg vcfeval output per pair (`*.tp.vcf.gz`,
  `*.tp-baseline.vcf.gz`, `*.fn.vcf.gz`, `*.fp.vcf.gz`, `*.summary.txt`, ROC tables over
  the tumor VAF).
- `benchmark/power.aggregated.tsv`, `benchmark/power.per_pair.tsv`: power-adjusted recall (same columns,
  stratum type `power_adjusted`); `NA` without read evidence.

</details>

The truth technology's calls are the baseline and every other technology is compared with
it. With `--comparison_mode all_pairs`, all other technologies are also compared pairwise
(the alphabetically first technology is the baseline). Matching uses `rtg vcfeval`
(haplotype-aware, `--squash-ploidy`, `--decompose`), so different representations of the
same variant still match.

`*.metrics.tsv` and `metrics.per_pair.tsv`:

| Column | Type | Meaning |
| --- | --- | --- |
| `sample`, `truth`, `query` | text | sample, truth (baseline) and query technology |
| `stratum_type`, `stratum` | text | see the strata below |
| `tp_baseline` | int | truth variants found by the query (baseline representation) |
| `fn` | int | truth variants missed by the query |
| `tp_call` | int | query calls that match a truth variant (call representation) |
| `fp` | int | query calls without a truth variant |
| `recall` | float | `tp_baseline / (tp_baseline + fn)` |
| `precision` | float | `tp_call / (tp_call + fp)` |
| `f1` | float | harmonic mean of recall and precision |

`metrics.aggregated.tsv` has `truth`, `query`, `stratum_type`, `stratum`, `n_samples`
(int), the four counts summed over the samples, `micro_recall`, `micro_precision`,
`micro_f1` (float, from the summed counts) and `macro_recall_mean`, `macro_recall_sd`,
`macro_precision_mean`, `macro_precision_sd`, `macro_f1_mean`, `macro_f1_sd` (float, mean
and standard deviation of the per-sample values; `NA` with fewer than 2 samples for sd).

| Stratum type | Strata | Meaning |
| --- | --- | --- |
| `overall` | `all` | all compared variants |
| `variant_type` | `SNV`, `INDEL` | SNV includes atomized MNVs; INDEL includes insertions, deletions and complex variants |
| `truth_vaf` | VAF bins (`--vaf_bins`), `unknown` | recall by the VAF of the truth call: the detection limit of the query |
| `query_vaf` | VAF bins, `unknown` | precision by the VAF of the query call: low-VAF false positives |
| `match` | `exact`, `position` | `bcftools isec`: same CHROM/POS/REF/ALT (per allele), same CHROM/POS (per position) |

`*.variants.tsv`:

| Column | Type | Meaning |
| --- | --- | --- |
| `sample`, `truth`, `query` | text | the pair |
| `status` | text | `TP`, `FN` (truth only) or `FP` (query only) |
| `chrom`, `pos`, `ref`, `alt` | contig, int, text, text | the variant (truth representation for TP/FN, query for FP) |
| `type`, `length` | text, int | `SNV`, `MNV`, `INS`, `DEL`, `COMPLEX`; `len(ALT) - len(REF)` |
| `truth_vaf`, `truth_dp`, `truth_alt` | float, int, int | values of the truth call (`NA` for FP) |
| `query_vaf`, `query_dp`, `query_alt` | float, int, int | values of the query call (`NA` for FN) |
| `query_filter` | text | FILTER of the query call |

### Discordance

<details markdown="1">
<summary>Output files</summary>

- `discordance/variants.annotated.tsv`: **every compared variant of every comparison**:
  the columns of `*.variants.tsv`, the discordance columns and the annotation columns below.
- `discordance/<sample>/<sample>_<truth>_vs_<query>.discordance.tsv`: the same per
  comparison, without annotation.
- `discordance/<sample>/evidence/<sample>_<technology>_<tumor|normal>.evidence.vcf.gz`:
  allele counts (FORMAT/AD, FORMAT/DP) of every alignment at all compared positions of the
  sample (only with `--read_evidence`).

</details>

Discordance columns:

| Column | Type | Meaning |
| --- | --- | --- |
| `class` | text | discordance class (below); `concordant` for TP |
| `other_in_normalized` | yes/no | the technology that misses the variant has it in its normalized VCF |
| `other_filter` | text | its FILTER (`NA` if not called) |
| `other_vaf_normalized` | float | its VAF |
| `other_tumor_dp`, `other_tumor_alt` | int | depth and alt reads of its tumor alignment (read evidence) |
| `normal_dp`, `normal_alt`, `normal_vaf` | int, int, float | summed over the normal alignments of both technologies |
| `query_power` | float | detection power of the query for truth variants |
| `detectable` | yes/no | `query_power >= --detection_power` |

The read-evidence columns (`other_tumor_*`, `normal_*`, `query_power`, `detectable`) are
left out with `--read_evidence false`.

A variant found in only one technology is explained from the point of view of the
technology that misses it (the query for an FN, the truth for an FP):

| Class | Meaning |
| --- | --- |
| `filtered_non_pass` | The other technology called it, but its FILTER is not PASS (`other_filter`); with `--filter_pass_only` such calls are not compared. |
| `present_but_filtered` | The other technology called it as PASS, but the call was removed by the uniform filters (`--min_vaf`, `--min_depth`, `--min_alt_reads`) or lies outside the comparison region. |
| `germline_in_normal` | The normal alignments show the alt allele (VAF ≥ `--germline_vaf`, ≥ `--evidence_min_alt_reads` reads): likely germline. |
| `low_coverage` | The other technology's tumor depth is below `--min_coverage_tumor`. |
| `subthreshold_evidence` | The other tumor alignment has alt reads, but the variant was not called. |
| `no_evidence` | The other tumor alignment has no alt read (MQ ≥ `--evidence_min_mapq`, BQ ≥ `--evidence_min_baseq`). |
| `absent` | No call and no read evidence to look further. |

Power-adjusted recall: for each truth variant, the detection power of the query is
`P(X ≥ --power_min_alt_reads)` with `X ~ Binomial(query depth, truth VAF)`; only variants
with power ≥ `--detection_power` count.

There is no separate `outside_callable` class: both technologies of a sample are compared
in the same region, so a call outside it is reported as `present_but_filtered`.

### Annotation

<details markdown="1">
<summary>Output files</summary>

- `annotation/vep/variants.vep.vcf.gz(.tbi)`: all compared variants (sites only) annotated
  by Ensembl VEP (INFO/CSQ), and `*.html`: the VEP summary (only with `--vep`).
- `annotation/repeats/<reference>.repeats.bed`: repeat (microsatellite) regions of the
  reference from msisensor-pro: contig, start, end, repeat unit (text), repeat times (int).
  Can be reused with `--repeats_bed` (only with `--annotate_repeats` without
  `--repeats_bed`).

</details>

Annotation columns of `discordance/variants.annotated.tsv`:

| Column | Type | Meaning |
| --- | --- | --- |
| `consequence` | text | most severe VEP consequence (SO terms joined by `&`) |
| `consequence_group` | text | `non_synonymous` (missense, nonsense, frameshift, in-frame, splice site, ...), `synonymous` or `other` (intronic, UTR, intergenic, ...) |
| `gene` | text | gene symbol (VEP SYMBOL), `.` if none |
| `impact` | text | `HIGH`, `MODERATE`, `LOW`, `MODIFIER` |
| `variant_class` | text | Sequence Ontology class (`SNV`, `insertion`, `deletion`, ...) |
| `biotype` | text | biotype of the transcript (e.g. `protein_coding`) |
| `in_repeat` | yes/no | the variant overlaps a repeat region |
| `repeat_unit`, `repeat_times` | text, int | repeat unit and number of repeats (`.` if none) |
| `giab_<name>` | yes/no | one column per GIAB file given with `--giab_<name>` (overlap with the variant's REF span) |
| `<stratum>` | yes/no | one column per `--stratification` BED |

### Known variants

<details markdown="1">
<summary>Output files</summary>

- `known_variants/required_summary.tsv`: per list, sample and technology (plus sample `all`
  per list and technology):

  | Column | Type | Meaning |
  | --- | --- | --- |
  | `list`, `sample`, `technology` | text | the list, the sample (`all` = summed over the samples), the technology |
  | `required`, `passed`, `failed` | int | number of listed variants; found (PASS); not found |
  | `recall` | float | `passed / required` |
  | `all_passed` | yes/no | every listed variant was found |

- `known_variants/required_variants.tsv`: every required variant per sample and technology:
  `list`, `sample`, `technology`, `label`, `gene`, `chrom` (contig), `pos` (int), `ref`,
  `alt`, `status` (`PASS`, `FAIL_filtered`, `FAIL_outside_region`, `FAIL_not_called`),
  `in_region` (yes/no), `called` (yes/no), `filter`, `vaf` (float), `depth` (int),
  `alt_reads` (int) of the call, `evidence_depth`, `evidence_alt_reads` (int, tumor
  alignment), `same_position_call` (other alleles called at the position, `.`).
- `known_variants/required/<list>/`: the same per list and sample.
- `known_variants/benchmark/known_benchmark.aggregated.tsv`, `.per_pair.tsv`, and
  `benchmark/<sample>/`: every technology against each external truth set, with the columns
  of the main benchmark (`truth` is the list name).
- `known_variants/lists/`: the normalized known-variant VCFs.

</details>

### Report

<details markdown="1">
<summary>Output files</summary>

- `report/variantcomp.enrichment.tsv`: are discordant variants different from concordant
  ones? Per comparison and for FN, FP and FN+FP:

  | Column | Type | Meaning |
  | --- | --- | --- |
  | `truth`, `query`, `discordant` | text | comparison; `FN`, `FP` or `FN+FP` |
  | `feature` | text | `indel`, `vaf_lt_<--low_vaf>`, `non_synonymous`, `synonymous`, `high_or_moderate_impact`, `in_repeat`, every stratification |
  | `discordant_with`, `discordant_without`, `concordant_with`, `concordant_without` | int | the 2 × 2 table |
  | `fraction_discordant`, `fraction_concordant` | float | share with the feature |
  | `odds_ratio`, `p_value`, `q_value` | float | Fisher's exact test; Benjamini-Hochberg q over all tests (significant: q < `--enrichment_fdr`) |

- `report/variantcomp.strata.tsv`: per comparison, annotation (`consequence_group`,
  `impact`, `variant_class`, `in_repeat`, every stratification) and value: `tp`, `fn`, `fp`
  (int), `recall`, `precision`, `discordance_rate` = (FN + FP) / (TP + FN + FP) (float).
- `report/variantcomp.genes.tsv`: per comparison and gene: `tp`, `fn`, `fp`,
  `discordance_rate`, `samples`, `odds_ratio`, `p_value`, `q_value` (discordant vs.
  concordant in the gene vs. all other genes); sorted by p-value.
- `report/variantcomp.classes.tsv`: per comparison and status (FN, FP): number of variants
  per discordance class.
- `report/plots/`: `recall_precision.png` (bars: all samples, dots: per sample),
  `recall_by_vaf.png`, `discordance_classes.png`, `discordance_by_annotation.png`,
  `enrichment.png` (log2 odds ratios).

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

The report contains: recall and precision per comparison and per sample, recall by truth
VAF, the discordance classes, discordance by annotation, the required variants, and one
table per step summary (input checks, normalization by `bcftools norm`, harmonization,
comparison region, callable regions), the mosdepth and VEP summaries and the software
versions.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Checked sheets with resolved paths: `samplesheet.validated.csv` (columns of the
    samplesheet plus `is_truth`), and if given `stratification.validated.csv` and
    `known.validated.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
