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
- [Report](#report) - the MultiQC report (main report), the table browser, the zip to send
  both, the consistency checks, enrichment statistics, gene and annotation tables, plots
- [MultiQC](#multiqc) - output files and configuration of the MultiQC report
- [Pipeline information](#pipeline-information) - run reports and checked sheets

Tools that produce each output folder (citations in [`CITATIONS.md`](../CITATIONS.md)):

| Folder | Tools |
| --- | --- |
| `input_checks/` | `check_reference.py` (pysam), `known_variants.py` (pysam) |
| `regions/` | `bedtools sort`, `bedtools multiinter`, `gawk`, `bedtools merge`, `bedtools intersect`, `bedtools subtract`, `region_summary.py` |
| `coverage/` | `check_reference.py`, `mosdepth`, `mosdepth_callable.py`, `bedtools sort`, `bedtools multiinter`, `bedtools intersect`, `region_summary.py`, `coverage_summary.py` |
| `vcf_preprocessed/` | `bcftools view`, `bcftools annotate`, `bcftools reheader`, `bcftools norm`, `harmonize_vcf.py`, `bcftools concat`, `variant_accounting.py` (pysam) |
| `benchmark/` | `rtg vcfeval`, `bcftools isec`, `benchmark_metrics.py`, `aggregate_metrics.py` |
| `discordance/` | `bcftools query`, `read_evidence.py` (pysam), `classify_discordance.py`, `annotate_variants.py` |
| `annotation/` | Ensembl VEP, `bcftools +split-vep`, msisensor-pro, `bedtools map`, `bedtools intersect` |
| `known_variants/` | `known_variants.py`, `bcftools norm`, `bcftools isec`, `rtg vcfeval`, `benchmark_metrics.py`, `known_summary.py` |
| `report/` | `variantcomp_report.py` (pandas, SciPy, Matplotlib), `consistency_checks.py`, `report_content.py` (pandas), datavzrd, `report_bundle.py` |
| `multiqc/` | MultiQC (with the sections of `report_content.py`; plots with the Plotly library bundled in MultiQC), `mqc_tables.py` |

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
| `records_read`, `records_kept` | int | known-variant lists: VCF records / list rows read and usable (a multiallelic record counts once) |
| `warning`, `error` | text | one row per message |

### Regions

<details markdown="1">
<summary>Output files</summary>

- `regions/common_targets.bed`: the comparison region, BED3 (contig, start, end): the
  bases targeted by every technology (WGS counts as genome-wide), intersected with
  `--regions_bed` and minus `--exclude_bed`; sorted, merged.
- `regions/common_targets.summary.tsv`: size of every target BED and of the common region.
- `regions/common_targets.steps.tsv`: size of the comparison region after every step of
  its construction (columns below).
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

`common_targets.steps.tsv` (one row per step, in order: `1_genome` (the compared contigs of
the reference), `2_targets_intersection` (intersection of all target BEDs; WGS has none and
doesn't restrict it), `3_regions_bed` (only with `--regions_bed`), `4_exclude_bed` (only with
`--exclude_bed`); overlapping intervals of a step are merged first):

| Column | Type | Meaning |
| --- | --- | --- |
| `step` | int | order of the step |
| `name` | text | step name (above) |
| `intervals` | int | intervals after merging |
| `bp` | int | size of the region after the step |
| `fraction_of_genome` | float | `bp / bp of 1_genome` |
| `fraction_of_previous` | float | `bp / bp of the previous step` (`NA` for the first) |
| `removed_bp` | int | bases removed by the step (`NA` for the first) |

A step can't add bases: `region_summary.py steps` stops if one does.

The run stops if the common region is smaller than `--min_region_bp`.

### Coverage

Only with `--callable_regions` (default; needs alignments for every samplesheet row).

<details markdown="1">
<summary>Output files</summary>

- `coverage/<sample>.callable.bed`: BED3 of the bases of the common region callable in all
  alignments of the sample (tumor depth ≥ `--min_coverage_tumor`, normal depth ≥
  `--min_coverage_normal` with `--require_normal_coverage`). The comparison of this sample
  is restricted to it. Built from the callable bins of every alignment clipped to the
  common region (`mosdepth_callable.py --clip`), intersected with `bedtools multiinter` and
  merged (see [usage](usage.md#coverage-and-callable-regions)).
- `coverage/callable_summary.tsv`: one row per sample (columns below).
- `coverage/coverage.coverage.tsv`: depth of every alignment in the common region and in
  the target region of its technology (columns below). The same values are columns of the
  MultiQC General Statistics table (`Depth common`, `Depth target`, `<t>X common`, `<t>X
  target`; the third threshold is shown, the others can be switched on).
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

`coverage.coverage.tsv` (mosdepth `--thresholds` of `--coverage_thresholds`; a second
mosdepth run per alignment on the target BED of its technology):

| Column | Type | Meaning |
| --- | --- | --- |
| `alignment` | text | `<sample>_<technology>_<tumor\|normal>` |
| `sample`, `technology`, `role` | text | samplesheet row and `tumor` / `normal` |
| `region` | text | `common` (common comparison region), `target` (the target BED of the technology) or `genome` (technology without target BED: the whole genome, from mosdepth's global distribution, 2-digit precision) |
| `bp` | int | size of the region |
| `mean_depth` | float | mean depth, weighted by interval length |
| `pct_ge_<t>X` | float | percent of the bases with depth ≥ t, one column per threshold |

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
- `vcf_preprocessed/accounting/<sample>_<technology>.accounting.tsv`: records of the input
  VCF kept, removed and added per step, with the reason (`variant_accounting.py`, columns
  below).
- `vcf_preprocessed/accounting/<sample>_<technology>.removed.tsv`: every removed record:
  `sample`, `technology`, `step_name`, `reason`, `chrom` (contig as in the input),
  `pos` (1-based), `ref`, `alt`, `filter`.
- `vcf_preprocessed/accounting/<sample>_<technology>.checks.tsv`: checks that the counts
  add up (`sample`, `technology`, `check`, `expected`, `observed`, `status` `ok`/`FAIL`);
  part of the [consistency checks](#report).

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
6. `bcftools norm`: multiallelic records split, indels left-aligned, MNVs and complex
   calls atomized (`--atom-overlaps .`: one ALT per atomized record, see
   [usage](usage.md#overlapping-alleles)), records with a REF mismatch removed, duplicates
   removed.
7. Harmonization (`harmonize_vcf.py`, see [usage](usage.md#harmonization-of-vaf-depth-and-alt-reads));
   records with a symbolic ALT or the overlapping allele `*` are dropped and counted.
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

`*.accounting.tsv`, one row per step, outcome and reason (the input VCF is read with pysam;
the steps are reproduced from the step outputs and logs):

| Column | Type | Meaning |
| --- | --- | --- |
| `sample`, `technology` | text | samplesheet row |
| `step`, `step_name` | int, text | 1 `input`, 2 `pass_split`, 3 `contigs`, 4 `gvcf_blocks`, 5 `normalization`, 6 `harmonization`, 7 `filters` |
| `outcome` | text | `input`, `removed`, `added` (records created by splitting) or `kept` (step 7: the compared records) |
| `reason` | text | `filter_not_pass` (FILTER not PASS/`.`), `mitochondrial`, `contig_not_compared`, `gvcf_block` (only symbolic `<NON_REF>`/`<*>` ALT), `ref_mismatch`, `duplicate`, `split` (added), `symbolic_allele`, `overlap_allele` (ALT `*` left after splitting: the position is covered by a variant of another record, e.g. a spanning deletion), `outside_region`, `vaf_below_min`, `depth_below_min`, `alt_reads_below_min`; `.` for input/kept |
| `records` | int | records |
| `snv`, `indel`, `other` | int | records by type (`other`: symbolic or mixed) |

A record failing several filters of step 7 is counted once, under the first reason in the
order above. Checks (`*.checks.tsv`): the records entering normalization equal the total of
the `bcftools norm` log; the records written by the harmonization equal the harmonized
VCF; the compared VCF holds exactly the records passing the filters and the region; input =
the removed records of steps 2–4 + the records entering normalization.

Keys of `*.harmonize_summary.tsv` (`key`, `value`): `records_read`, `records_written`,
`records_dropped_symbolic`, `records_dropped_overlap` (ALT `*`), `vaf_available` (yes/no), and per sample (prefix `T_`, `N_`):
`source_AF`, `source_AD`, `source_zero_depth` (AD 0,0: no reads), `source_none` (no
VAF), `depth_DP`, `depth_AD_sum`, `depth_none`, `af_nan` (FORMAT/AF was NaN).

### Benchmark

<details markdown="1">
<summary>Output files</summary>

Every matching method has its own folder, named by the tool:

- `benchmark/rtgtools_vcfeval/<sample>/<pair>/`: rtg vcfeval output (`*.tp-baseline.vcf.gz`:
  truth variants found, `*.fn.vcf.gz`: truth variants missed, `*.tp.vcf.gz`: query calls
  matching a truth variant, `*.fp.vcf.gz`: query calls without, `*.summary.txt`, ROC tables
  over the tumor VAF). `<pair>` = `<sample>_<truth>_vs_<query>`.
- `benchmark/bcftools_isec_exact/<sample>/<pair>/`: `bcftools isec --collapse none` (same
  CHROM, POS, REF, ALT): `sites.txt` (CHROM, POS, REF, ALT, flags truth/query: `10` truth
  only, `01` query only, `11` both), `0000.vcf.gz` (truth only), `0001.vcf.gz` (query
  only), `0002.vcf.gz`/`0003.vcf.gz` (in both, truth and query records), `README.txt`.
- `benchmark/bcftools_isec_position/<sample>/<pair>/`: the same with `--collapse all`
  (same CHROM and POS, any allele).
- `benchmark/metrics/`: recall, precision and F1 of every method:
  - `metrics.aggregated.tsv`: **main result**, per method, truth/query technology and
    stratum, over all samples (columns below); `metrics.per_pair.tsv`: every sample and
    technology pair;
  - `metrics.rtgtools_vcfeval.{aggregated,per_pair}.tsv`,
    `metrics.bcftools_isec_exact.{aggregated,per_pair}.tsv`,
    `metrics.bcftools_isec_position.{aggregated,per_pair}.tsv`: the same, one method each;
  - `power.aggregated.tsv`, `power.per_pair.tsv` (and `power.rtgtools_vcfeval.*`):
    power-adjusted recall (stratum type `power_adjusted`); `NA` without read evidence;
  - `<sample>/<pair>.metrics.tsv`: the metrics of one pair; `<sample>/<pair>.variants.tsv`:
    every compared variant of the pair with its status per method (columns below).

</details>

The truth technology's calls are the baseline and every other technology is compared with
it. With `--comparison_mode all_pairs`, all other technologies are also compared pairwise
(the alphabetically first technology is the baseline). Three matching methods are
computed from the same compared VCFs (same records, same region):

- **RTG vcfeval** (main result; `--squash-ploidy`, `--decompose`): haplotype-aware, so
  different representations of the same haplotype match (e.g. two SNVs against a
  deletion + insertion that give the same sequence). All strata (variant type, VAF) and the
  discordance analysis use this method.
- **bcftools isec, same allele** (`--collapse none`): a match needs the same CHROM, POS, REF
  and ALT, counted per record.
- **bcftools isec, same position** (`--collapse all`): a match needs the same CHROM and POS,
  counted per position (several alleles at one position count once).

All VCFs are normalized the same way before (split, left-aligned, MNVs atomized, complex
variants decomposed by `bcftools norm --atomize`), so most representation differences are
already removed and the methods often give the same numbers. They differ where a haplotype
is called as different events (vcfeval matches, same-allele isec doesn't), and where the
technologies call different alleles at the same position (only same-position isec
matches). The column `isec_exact` / `isec_position` of `*.variants.tsv` shows for every
variant where the methods disagree. Tested with the synthetic events E06 and E15
(`design.tsv` of the test data).

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

`metrics.per_pair.tsv` and `metrics.aggregated.tsv` start with the column `method`
(`rtgtools_vcfeval`, `bcftools_isec_exact`, `bcftools_isec_position`); the per-pair
`<pair>.metrics.tsv` has no method column (the stratum shows it).

`metrics.aggregated.tsv` has `method`, `truth`, `query`, `stratum_type`, `stratum`, `n_samples`
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
| `match` | `exact`, `position` | `bcftools isec` (method `bcftools_isec_exact` / `bcftools_isec_position`): same CHROM/POS/REF/ALT (per allele), same CHROM/POS (per position); all other strata are RTG vcfeval |

`*.variants.tsv`:

| Column | Type | Meaning |
| --- | --- | --- |
| `sample`, `truth`, `query` | text | the pair |
| `status` | text | RTG vcfeval: `TP`, `FN` (truth only) or `FP` (query only) |
| `isec_exact`, `isec_position` | text | the status of the same record with bcftools isec, same allele / same position (`TP`, `FN`, `FP`; `NA` for a record vcfeval decomposed, which isec doesn't know) |
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
- `discordance/<sample>/read_evidence/<sample>_<technology>_<tumor|normal>.evidence.vcf.gz`:
  read counts of every alignment at every compared variant of the sample (only with
  `--read_evidence`), one record per variant: FORMAT/DP (reads counted), FORMAT/AD (REF
  reads, ALT reads), FORMAT/OT (reads with another allele); see Read evidence below.

</details>

Discordance columns:

| Column | Type | Meaning |
| --- | --- | --- |
| `class` | text | the main reason (below): exactly one per FN/FP; `concordant` for TP. Not shown in the table browser: it is the first entry of `reasons` |
| `reasons` | text | all reasons, comma-separated: the main reason first, then `germline_in_normal` and the context reasons (below); `concordant` for TP |
| `other_in_normalized` | yes/no | the technology that misses the variant has it in its normalized VCF (called) |
| `other_filter` | text | its FILTER (`NA` if not called) |
| `other_vaf_normalized` | float | its VAF |
| `other_tumor_dp`, `other_tumor_alt`, `other_tumor_other` | int | reads counted, alt reads and reads with another allele at the variant in its tumor alignment (read evidence) |
| `normal_dp`, `normal_alt`, `normal_vaf` | int, int, float | summed over the normal alignments of both technologies |
| `query_power` | float | detection power of the query for truth variants |
| `detectable` | yes/no | `query_power >= --detection_power` |

The read-evidence columns (`other_tumor_*`, `normal_*`, `query_power`, `detectable`) are
left out with `--read_evidence false`.

A variant found in only one technology is explained from the point of view of the
technology that misses it (the query for an FN, the truth for an FP).

**Main reason** (`class`, first entry of `reasons`): exactly one, they exclude each other.
It describes the call and the reads of the technology that misses the variant:

| Main reason | Meaning |
| --- | --- |
| `filtered_non_pass` | Called, but its FILTER is not PASS (`other_filter`); with `--filter_pass_only` such calls are not compared. |
| `present_but_filtered` | Called with PASS, but the call was removed by the uniform filters (`--min_vaf`, `--min_depth`, `--min_alt_reads`) or lies outside the comparison region. |
| `absent` | Not called, and no tumor alignment to look further (none given, or `--read_evidence false`). |
| `low_coverage` | Not called; the tumor alignment has fewer than `--min_coverage_tumor` reads at the variant. |
| `subthreshold_evidence` | Not called, although the tumor alignment has enough reads and at least one read with the alt allele. |
| `no_evidence` | Not called; the tumor alignment has enough reads but none with the alt allele. |

The read-evidence reasons (`low_coverage`, `subthreshold_evidence`, `no_evidence`) are only
used for variants that the other technology did not call: a call (PASS or not) is
explained by its call. A variant can't be `filtered_non_pass` and `subthreshold_evidence`
at the same time.

**Further reasons** (after the main reason, any number):

| Reason | Meaning |
| --- | --- |
| `germline_in_normal` | The normal alignments show the alt allele (alt reads ≥ `--evidence_min_alt_reads`, VAF ≥ `--germline_vaf`): likely germline (or a contaminated normal). |
| `indel` | The variant is an insertion, deletion or complex variant. |
| `low_vaf` | The VAF of the call that exists (truth call of an FN, query call of an FP) is below `--low_vaf`. |
| `in_repeat` | The variant lies in a repeat (microsatellite) region of the reference (column `in_repeat`). |
| `giab_<name>`, `<stratum>` | The variant overlaps this GIAB stratification (`--giab_<name>`) or `--stratification` BED (column of the same name). |

The context reasons (`indel` and below) come from the annotation (`annotate_variants.py
merge --context-reasons`). An insertion counts as overlapping a region if its anchor base
or the base after it does (the inserted bases go between them; a left-aligned insertion's
anchor lies just before a repeat). The statistical test of these features is the
enrichment table of the report.

**Read evidence** (`read_evidence.py`, one run per alignment). For every compared variant
of the sample (and every known variant), the reads of the alignment are counted (mapping
quality ≥ `--evidence_min_mapq`, no secondary, supplementary, QC-failed or duplicate
reads; duplicates are kept for amplicon libraries with `--amplicon_keep_duplicates`; at
most `--evidence_max_depth` reads):

- SNV / MNV: the read's bases at the variant (base quality ≥ `--evidence_min_baseq`) are
  compared with ALT and REF.
- Indels and complex variants: a read supports the variant if its sequence (soft-clipped
  bases included) contains the haplotype of the window around the variant with the variant
  applied; the window spans every position the indel can be shifted to in a repeat plus
  `--evidence_indel_flank` bases on each side. It supports REF if it contains the reference
  window; a read aligned across the window that contains neither counts as other (another
  indel, or a base error in the window). So the count doesn't depend on where or whether
  the aligner placed the indel. Reads that end inside the window aren't counted, so the
  depth of an indel is lower than the depth of an SNV at the same position.

`bcftools mpileup`, used before, discovers indel alleles itself: in repeats it reports
other lengths than the caller and it misses long deletions, so the alt reads of called
indels were counted as 0 and such variants were wrongly explained as `no_evidence`. An
end-of-run check now requires that most TP variants (`--evidence_check_min_fraction`, SNVs
and indels separately) have alt reads in the query tumor alignment.

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
| `transcript` | text | transcript of the most severe consequence (VEP Feature; `.` without `--hgvs`) |
| `hgvsc`, `hgvsp` | text | coding and protein change (HGVS, VEP `--hgvs`, default) without the transcript/protein ID, e.g. `c.35G>T`, `p.Gly12Val`; `.` if none (non-coding, intergenic) or without `--hgvs` |
| `dbsnp`, `cosmic` | text | IDs of known variants at the allele (VEP `--check_existing`, only with `--vep_cache`): rs IDs; COSV/COSM IDs; comma-separated, `.` if none |
| `clin_sig` | text | ClinVar significance (only with `--vep_cache`), comma-separated, `.` if none |
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
- `known_variants/benchmark/`: every technology against each external truth set, in the
  layout of `benchmark/` (`rtgtools_vcfeval/`, `bcftools_isec_exact/`,
  `bcftools_isec_position/`, `metrics/known_benchmark.*.tsv`, `metrics/<sample>/`). All three
  methods are restricted to the same region: the confident regions of the list intersected
  with the comparison region of the sample (vcfeval `--bed-regions`, bcftools isec `-R`),
  so a call that overlaps the region only partly is counted by every method; with the columns
  of the main benchmark (`truth` is the list name).
- `known_variants/lists/`: the normalized known-variant VCFs: the records of the input lists
  (ID, QUAL, FILTER, INFO kept), on the reference contigs, split and left-aligned, with
  `INFO/VC_T_VAF` (truth VAF, if known) and, for CSV/TSV lists, `INFO/KV_GENE`.

</details>

### Report

**Start here: `main.html`** at the top of the output folder. It opens the table browser
on its first page, which shows the MultiQC report; the menu at the top of the page
switches to the variant tables, and the links of the MultiQC report open them in a new
tab. `main.html` also links to the MultiQC report and the tables on their own.

The **MultiQC report** `multiqc/multiqc_report.html` is the main report of a run (decision
D37): one self-contained HTML file, interactive (hover, zoom, export of every plot; the
toolbox hides, highlights and renames samples in all plots and tables). It starts with the
results and ends with the details and the quality control:

1. **Summary**: banner of the consistency checks; cards with recall, precision, F1, FN and
   FP of every query technology (summed over the samples), with buttons for the matching
   method (RTG vcfeval, bcftools isec same allele, bcftools isec same position); cards of
   the run (samples, technologies, compared region, truth variants, required variants);
   key findings as a list; links to the variant tables.
2. **Introduction**: this run (question, samples, callers, compared region, filters) and
   how the technologies are compared, as lists; one row per technology (role, library type,
   caller, target BED and its size).
3. **Recall and precision per sample**: MultiQC scatter plot, both axes 0 to 1, one dot per
   sample and query technology (hover: sample name), one dataset per matching method.
4. **Overlap of the technologies per sample**: a Venn diagram per sample and comparison:
   variants only in the truth technology (FN), truth variants in both (TP), calls only in
   the query technology (FP), for every matching method (buttons); with vcfeval the number
   of query calls is added in brackets if it differs from the truth records.
   **Truth variants found**: vertical bars per sample, label `found/total`, per method.
5. **Known variants**: required variants as cards (found / required per technology) and a
   PASS/FAIL matrix (variant × technology × sample, details on hover); the benchmark
   against an external truth set as cards and bars per sample.
6. **Why variants differ**: discordant variants per sample (vertical bars) by main reason
   (exactly one per variant, stacked); views FN + FP, FN only (in the truth technology
   only), FP only (in the query technology only); and the number of discordant variants
   per sample with each further reason (germline, indel, low VAF, repeat, strata; grouped,
   a variant can have several).
7. **Discordant variants in detail** (per comparison): what is found in one technology only
   and why, as a list (main reasons, variants with several reasons, indel share, median VAF,
   protein-altering variants, genes); share of every category of type, consequence,
   impact, VAF, repeat, GIAB strata and main reason among TP, FN and FP; VAF of TP, FN and
   FP; the combinations of reasons (main reason + further reasons); the protein-altering discordant variants, the genes and the
   enrichment tests (top 15 rows each, the rest in the table browser).
8. **Comparison metrics per sample**: recall, precision and F1 by variant type, recall by
   truth VAF (one line per sample), power-adjusted recall; **Metrics per sample**: table
   of counts and the metrics of all three matching methods.
9. **Compared genomic region**: cards (genome, every target BED with its technology and
   role, common region, callable per sample) in Mbp; the steps as a flow; vertical bars
   (log or linear scale); the callable region per sample.
10. **Input variant sets and their filtering**: records per VCF, compared and removed per
    reason (counts or percent).
11. **Coverage**: mean depth and percent of bases ≥ each threshold per alignment, in the
    common and the target region.
12. The QC of the steps: input checks, normalization (`bcftools norm`), harmonization,
    mosdepth, VEP; **Consistency checks** (passed/failed per category, failed checks listed);
    methods description, software versions and parameters.

No table with variants in rows has more than 15 rows in MultiQC; the sections link to the
table browser. The bar plots are vertical with the samples on the x-axis: MultiQC's own bar
plots are horizontal only, so these plots are drawn with the Plotly library that MultiQC
bundles; they follow the toolbox like the MultiQC plots. The General Statistics table is
switched off (`skip_generalstats`), because MultiQC always puts it above all sections;
coverage is in the section Coverage.

The **table browser** `report/tables/index.html` (datavzrd) holds the variant tables. Every
view starts with a description (what the table shows, how to read it, links to the MultiQC
report and related tables, the column glossary). Views: the MultiQC report (first page),
truth variants, discordant variants (all FN and FP, with the class
explanations), all compared variants, genes, enrichment, metrics per sample, required
variants, known-variant benchmark, removed records, preprocessing accounting, region steps,
coverage, consistency checks, samplesheet. In the variant tables:

- the first columns are `variant` (gene and protein change, e.g. `KRAS p.Gly12Val`; HGVS.c
  or the consequence if there is no protein change), `gene`, `hgvsp`, `consequence`,
  `impact`, then sample, status, class and reasons, then position, alleles, VAF and the rest;
- coloured columns: `status` (TP green, FN orange, FP purple), `reasons_<query>`,
  `reasons` (one coloured pill per reason), `found_<query>` and other yes/no columns,
  `impact`, metrics (red to green), q-values; VAF columns as ticks on 0 to 1;
- links: `gene` to GeneCards and Ensembl, `locus` (`chrom:pos`) to the Ensembl genome
  browser (`--vep_genome`), `dbsnp` to dbSNP, `cosmic` to COSMIC; rows link to the truth
  variant (`variant_id` = `sample:chrom:pos:ref>alt`) and to the gene table;
- the less important columns (HGVS.c, transcript, biotype, evidence counts, ...) are in the
  detail view of a row (`+`);
- the variant tables can be downloaded as Excel.

To send the reports, use **`report/variantcomp_report.zip`**: the start page `main.html`,
the MultiQC report and the table browser in the layout of the output folder (`multiqc/`,
`report/tables/`); unzip it and open `main.html`. The MultiQC
report alone is a single file and can be sent as it is (its links to the tables then do
not work).

<details markdown="1">
<summary>Output files</summary>

- `main.html` (top of the output folder): the start page (above).
- `report/tables/`: the datavzrd table browser (`index.html` opens its first page, the
  MultiQC report).
- `report/variantcomp_report.zip`: start page, MultiQC report and table browser in one zip
  (above).
- `report/consistency_checks.tsv`: do the counts of the steps add up? (columns below)
- `report/variantcomp.truth_variants.tsv`: every compared truth variant per sample:
  `sample`, `truth`, `chrom`, `pos`, `ref`, `alt`, `found_by` (int: query technologies that
  found it), `queries` (int: query technologies compared), the annotation of the truth call
  (`type`, `length`, `truth_vaf`, `truth_dp`, `gene`, `hgvsp`, `hgvsc`, `consequence`,
  `consequence_group`, `impact`, `variant_class`, `in_repeat`, `transcript`, `dbsnp`,
  `cosmic`, `clin_sig`, `biotype`, if present) and per query technology
  `found_<query>` (yes/no) and `reasons_<query>` (all reasons of a miss, main reason
  first, `.` if found).
- `report/variantcomp.annotation_counts.tsv`: counts per `sample`, `technology`, `set`
  (`truth_variants`: TP+FN of the truth technology, each variant once; `query_calls`: TP+FP
  of a query technology; `missed_FN`; `extra_FP`), `annotation` (`type`,
  `consequence_group`, `impact`, `variant_class`, `in_repeat`), `value` and `count` (int).
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
  Limitation: `tp` counts the TP rows of the variant table, which are in the truth
  representation; `precision` = tp / (tp + fp) therefore equals the precision of the
  metrics only if every truth record matched exactly one query call (TP_baseline =
  TP_call, true for all tested data). With vcfeval matching one truth record to several
  query records (e.g. a complex variant against SNVs), the strata and gene precision is
  approximate; the metrics tables are exact.
- `report/variantcomp.genes.tsv`: per comparison and gene: `tp`, `fn`, `fp`,
  `discordance_rate`, `samples`, `odds_ratio`, `p_value`, `q_value` (discordant vs.
  concordant in the gene vs. all other genes); sorted by p-value.
- `report/variantcomp.classes.tsv`: per comparison and status (FN, FP): number of variants
  per discordance class.
- `report/plots/`: `recall_precision.png` (bars: all samples, dots: per sample),
  `recall_by_vaf.png`, `discordance_classes.png`, `discordance_by_annotation.png`,
  `enrichment.png` (log2 odds ratios).

</details>

`consistency_checks.tsv` (`consistency_checks.py check`, at the end of the run):

| Column | Type | Meaning |
| --- | --- | --- |
| `category` | text | `accounting` (the checks of every `*.checks.tsv`), `benchmark` (per pair: TP_baseline + FN ≥ compared truth records and TP_call + FP ≥ compared query records, equal unless `rtg vcfeval --decompose` splits complex records; FN and FP rows of the variant table = FN and FP of the metrics), `aggregation` (summed counts = sum over the samples), `discordance` (classes add up to FN / FP; every FN/FP has a known class), `regions` (last region step = common region; callable ≤ common region), `known` (passed + failed = required) |
| `scope` | text | the sample, VCF or comparison checked |
| `check` | text | what is checked |
| `expected`, `observed` | text | the values compared |
| `status` | text | `ok` or `FAIL` |

**A failed check stops the run** (`CHECKS_GATE`, exit status 1) after the MultiQC report,
the table browser and the zip were written, so the failure can be inspected
there. A failed check means the results may be wrong: please report it as a bug.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

The MultiQC report is the main report of the run; its sections are described under
[Report](#report). `assets/multiqc_config.yml` sets the order of the sections, the title and
fixed output names (`multiqc_report.html`, `multiqc_data/`, `multiqc_plots/`; the table
browser and the zip link to them).

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
