# core-unit-bioinformatics/workflow-nxf-somatic-variantcomp: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Regions](#regions) - Reference check of all target BEDs and the common comparison region
- [Coverage](#coverage) - Depth of every alignment and per-sample callable regions
- [Preprocessed VCFs](#preprocessed-vcfs) - Reference check, normalization, harmonization and uniform filtering of every VCF
- [Benchmark](#benchmark) - Recall, precision and F1 of every technology against the truth technology
- [Discordance](#discordance) - Why variants differ: classification, read evidence, genes and stratifications
- [Known variants](#known-variants) - Required variants (PASS/FAIL per technology) and benchmarks against external truth sets
- [Report](#report) - Enrichment statistics, gene summary and plots
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

Tools that produce each output folder (citations in [`CITATIONS.md`](../CITATIONS.md)):

| Folder | Tools |
| --- | --- |
| `regions/` | `check_reference.py`, `bedtools sort`, `bedtools multiinter`, `gawk`, `bedtools merge`, `bedtools intersect`, `bedtools subtract`, `region_summary.py` |
| `coverage/` | `samtools idxstats`, `check_reference.py`, `mosdepth`, `mosdepth_callable.py`, `bedtools sort`, `bedtools multiinter`, `bedtools intersect`, `region_summary.py` |
| `vcf_preprocessed/` | `bcftools view`, `bcftools query`, `check_reference.py`, `bcftools annotate`, `bcftools reheader`, `bcftools norm`, `harmonize_vcf.py` |
| `benchmark/` | `rtg vcfeval`, `bcftools isec`, `benchmark_metrics.py`, `aggregate_metrics.py` |
| `discordance/` | `bcftools query`, `bcftools mpileup` + `bcftools call`, `bcftools norm`, `classify_discordance.py`, `annotate_variants.py`, `bedtools sort`, `bedtools map`, `bedtools intersect` |
| `known_variants/` | `known_variants.py`, `bcftools norm`, `check_reference.py`, `bcftools isec`, `rtg vcfeval`, `benchmark_metrics.py`, `known_summary.py` |
| `report/` | `variantcomp_report.py` (pandas, SciPy, Matplotlib) |
| `multiqc/` | MultiQC |

### Regions

<details markdown="1">
<summary>Output files</summary>

- `regions/`
  - `common_targets.bed`: the comparison region, i.e. the intersection of the target BEDs of
    all targeted technologies (WGS counts as genome-wide), intersected with `--regions_bed`
    and minus `--exclude_bed` if given. Sorted in reference order and merged.
  - `common_targets.summary.tsv`: intervals and bp per technology, how many of them lie in
    the common region, and the size of the common region.
  - `targets/<technology>.bed`: target BED of each technology after the reference check
    (contigs renamed to the reference, unknown contigs dropped, sorted, merged).
  - `reference_check/<name>.bed.reference_check.tsv`: result of the reference check per BED
    (`status` = `ok`, `fixed` or `error`; renamed contigs, dropped intervals).

</details>

Every BED is checked against the reference `.fai`:

- **Fixed automatically:** contig names that differ only by the `chr` prefix (`1`
  ↔ `chr1`) or the mitochondrial name (`MT` ↔ `chrM`) are renamed. Intervals on contigs
  absent from the reference are dropped and reported.
- **Error:** intervals that end beyond their contig indicate a different genome build. The
  pipeline stops with a message naming the file.

### Coverage

<details markdown="1">
<summary>Output files</summary>

- `coverage/`
  - `<sample>.callable.bed`: bases of the common region that are callable in all alignments
    of the sample. Only for samples with alignments for every technology.
  - `callable_summary.tsv`: per sample, the number of alignments, the mode, the threshold,
    the size of the common region and of the callable region.
  - `<sample>/<sample>_<technology>_<tumor|normal>.*`: mosdepth outputs per alignment
    (`*.mosdepth.summary.txt`, `*.mosdepth.global.dist.txt`, `*.mosdepth.region.dist.txt`,
    `*.regions.bed.gz` = mean depth per common-region interval,
    `*.quantized.bed.gz` = `LOW_COVERAGE`/`CALLABLE` bins).
  - `reference_check/<...>.bam.reference_check.tsv`: reference check of each alignment header.

</details>

The mosdepth results are also shown in the MultiQC report.

### Preprocessed VCFs

<details markdown="1">
<summary>Output files</summary>

- `vcf_preprocessed/<sample>/`
  - `<sample>_<technology>.normalized.vcf.gz(.tbi)`: all records on reference contigs after
    normalization and harmonization, before filtering and region restriction.
  - `<sample>_<technology>.compare.vcf.gz(.tbi)`: the calls that enter the comparison: the
    uniform filters are applied and only records in the comparison region are kept.
  - `<sample>_<technology>.normalized.harmonize_summary.tsv`: how often each harmonization
    rule was used per sample (`T_source_AD`, `T_source_AF_DP`, `T_source_none`, depth from
    `DP` or `sum(AD)`, ...), records read, written and dropped (symbolic ALT); see
    [usage: harmonization](usage.md#harmonization-of-vaf-depth-and-alt-reads).
- `vcf_preprocessed/reference_check/`
  - `<sample>_<technology>.vcf.reference_check.tsv`: contig check per VCF (`status` = `ok`,
    `fixed` or `error`; renamed contigs, header length fixes, records on unknown contigs,
    sample names).
  - `<sample>_<technology>.ref_alleles.tsv`: REF-allele check, from the `bcftools norm`
    summary (records, REF mismatches removed by `--check-ref x`, mismatch rate, threshold).
  - `<sample>_<technology>.norm.log`: the `bcftools norm` log.

</details>

Every VCF goes through the same steps:

1. Reference check (contig names and lengths, REF alleles, sample names).
2. Tumor and normal columns selected in this order; recompressed with BGZF.
3. Contigs renamed to the reference (`chr` prefix, `MT`/`chrM`).
4. Records on non-reference contigs and pure gVCF reference blocks removed.
5. Contig header lines rewritten from the reference `.fai`.
6. `bcftools norm`: multiallelic records split, indels left-aligned, MNVs atomized into
   SNVs, records with a REF mismatch dropped, duplicates removed.
7. Harmonization: caller-independent INFO tags computed from FORMAT/AD, else FORMAT/AF +
   FORMAT/DP (`VC_T_VAF`, `VC_T_DP`, `VC_T_ALT`, `VC_T_SRC`, the same `VC_N_*` for the
   normal, `VC_TYPE`, `VC_LEN`); the rules are described in
   [usage](usage.md#harmonization-of-vaf-depth-and-alt-reads). Every record is a call of the
   tumor (GT tumor `0/1`, normal `0/0`), and remaining symbolic alleles are dropped.
8. Uniform filters (`--filter_pass_only`, `--min_vaf`, `--min_depth`, `--min_alt_reads`) and
   restriction to the comparison region (records overlapping it are kept).

### Benchmark

<details markdown="1">
<summary>Output files</summary>

- `benchmark/`
  - `metrics.aggregated.tsv`: **main result**. Per truth/query technology and stratum, the
    summed counts over all samples and the resulting `micro_recall`, `micro_precision` and
    `micro_f1`, and the mean and standard deviation of the per-sample values
    (`macro_*_mean`, `macro_*_sd`). `n_samples` is the number of samples.
  - `metrics.per_pair.tsv`: the metrics of every sample and technology pair.
- `benchmark/<sample>/`
  - `<sample>_<truth>_vs_<query>.metrics.tsv`: metrics of one pair.
  - `<sample>_<truth>_vs_<query>.variants.tsv`: every compared variant with status `TP`,
    `FN` (in truth, missed by the query) or `FP` (called by the query only), type, length
    and VAF, depth and alt reads in both technologies.
  - `vcfeval/`: rtg vcfeval output per pair (`*.tp.vcf.gz`, `*.tp-baseline.vcf.gz`,
    `*.fn.vcf.gz`, `*.fp.vcf.gz`, `*.summary.txt`, ROC tables over tumor VAF thresholds).
- `benchmark/power.aggregated.tsv`, `benchmark/power.per_pair.tsv`: power-adjusted recall
  (see [Discordance](#discordance)); `NA` without query alignments.

</details>

The truth technology's calls are the baseline and every other technology is compared with
it. With `--comparison_mode all_pairs`, all other technologies are also compared pairwise;
the alphabetically first technology is the baseline.

Matching uses `rtg vcfeval` (haplotype-aware, `--squash-ploidy`, `--decompose`, all records
of the compare VCFs), so different representations of the same variant still match.

| Stratum type | Strata | Meaning |
| --- | --- | --- |
| `overall` | `all` | All compared variants. |
| `variant_type` | `SNV`, `INDEL` | SNV includes atomized MNVs; INDEL includes insertions, deletions and complex variants. |
| `truth_vaf` | VAF bins (`--vaf_bins`) | Recall by the VAF of the truth call. Shows the detection limit of the query technology. |
| `query_vaf` | VAF bins | Precision by the VAF of the query call. Low-VAF false positives show up here. |
| `match` | `exact`, `position` | Set comparison of the two VCFs (`bcftools isec`) by `CHROM/POS/REF/ALT` (counted per allele) and by `CHROM/POS` only (counted per position: several alleles at one position count once). Differences to `overall` show how much discordance is representation or a different allele at the same position. |

`recall = TP_baseline / (TP_baseline + FN)`, `precision = TP_call / (TP_call + FP)`, and
`F1` is their harmonic mean. Undefined values are `NA`.

### Discordance

<details markdown="1">
<summary>Output files</summary>

- `discordance/variants.annotated.tsv`: **every compared variant of every comparison** (TP,
  FN, FP). The columns are the per-variant columns of the benchmark, plus:
  - `class` of each FN/FP (see below; `concordant` for TP);
  - `other_in_normalized`, `other_filter`, `other_vaf_normalized`: the call of the other
    technology before filtering, if any;
  - `other_tumor_dp`, `other_tumor_alt`: depth and alt reads in the other technology's tumor
    alignment;
  - `normal_dp`, `normal_alt`, `normal_vaf`: summed over the normal alignments of both
    technologies;
  - `query_power`, `detectable`: detection power of the query technology for truth variants;
  - `gene` (with `--gene_bed`) and one `yes`/`no` column per stratification
    (`--stratification`).
- `discordance/<sample>/<sample>_<truth>_vs_<query>.discordance.tsv`: the same, per
  comparison, without annotation.
- `discordance/<sample>/evidence/*.evidence.vcf.gz`: `bcftools mpileup` allele counts of
  every alignment at all compared positions of the sample.
- `discordance/reference_check/`: reference check of the gene and stratification BEDs.

</details>

A variant found in only one technology is explained from the point of view of the
technology that misses it (the query for an FN, the truth for an FP):

| Class | Meaning |
| --- | --- |
| `present_but_filtered` | The other technology called it, but the call failed FILTER or the uniform thresholds. |
| `germline_in_normal` | A normal alignment shows the alt allele (VAF ≥ `--germline_vaf`, ≥ `--evidence_min_alt_reads` reads): likely germline. |
| `low_coverage` | The other technology's tumor depth is below `--min_coverage`. |
| `subthreshold_evidence` | The other tumor alignment has alt reads, but the variant was not called. |
| `no_evidence` | The other tumor alignment has no alt read (MQ ≥ 20, BQ ≥ 13). |
| `absent` | No call and no alignment to look further. |

Power-adjusted recall: for each truth variant, the detection power of the query technology
is `P(X ≥ --power_min_alt_reads)` with `X ~ Binomial(query depth, truth VAF)`. Only variants
with power ≥ `--detection_power` count. This separates "the technology could not have
seen it" (depth/VAF) from "it should have seen it".

Read evidence needs the alignments to use the reference contig names.

There is no separate `outside_callable` class: both technologies of a sample are compared
in the same region (common targets, intersected with the callable region of the sample), so
a call of the other technology outside that region is in its normalized VCF but not
compared, and is reported as `present_but_filtered` (`other_filter` shows its FILTER;
`other_in_normalized` = `yes`).

### Known variants

<details markdown="1">
<summary>Output files</summary>

- `known_variants/`
  - `required_summary.tsv`: per list, sample and technology, the number of required
    variants, how many passed and `all_passed` (`yes`/`no`).
  - `required_variants.tsv`: every required variant per technology, with `status`
    (`PASS`, `FAIL_filtered`, `FAIL_outside_region`, `FAIL_not_called`), FILTER, VAF,
    depth, alt reads, read evidence and calls at the same position with other alleles.
  - `required/<list>/`: the same per list and sample.
  - `benchmark/known_benchmark.aggregated.tsv`, `benchmark/known_benchmark.per_pair.tsv`:
    metrics of every technology against each external truth set (same format as the main
    benchmark; `truth` is the list name).
  - `benchmark/<sample>/`: per-pair metrics, variant tables and vcfeval output.
  - `lists/`: the normalized known-variant VCFs per sample.
  - `reference_check/`: reference check of the list regions (`*.reference_check.tsv`) and
    of the REF alleles of the listed variants (`*.ref_alleles.tsv`; any mismatch is an
    error).

</details>

### Report

<details markdown="1">
<summary>Output files</summary>

- `report/`
  - `variantcomp.enrichment.tsv`: are discordant variants different from concordant ones?
    For every comparison and for FN, FP and FN+FP separately, the counts with and without a
    feature, Fisher's exact test (odds ratio, p) and the Benjamini-Hochberg q-value over all
    tests. Features: `indel`, `vaf_lt_<--low_vaf>` (truth VAF for FN/TP, query VAF for FP)
    and every stratification.
  - `variantcomp.genes.tsv`: per comparison and gene, TP/FN/FP counts, discordance rate,
    number of samples and Fisher's exact test of discordant vs. concordant variants in the
    gene vs. all other genes (BH q-value). Sorted by p-value.
  - `variantcomp.classes.tsv`: discordance classes per comparison and status.
  - `plots/`: `recall_precision.png` (bars: all samples, dots: per sample),
    `recall_by_vaf.png`, `discordance_classes.png`, `enrichment.png` (log2 odds ratios).

</details>

The MultiQC report starts with these sections: recall and precision per comparison and per
sample, recall by truth VAF and the discordance classes.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Validated samplesheet with resolved paths: `samplesheet.validated.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
