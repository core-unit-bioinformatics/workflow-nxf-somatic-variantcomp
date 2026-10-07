# core-unit-bioinformatics/workflow-nxf-somatic-variantcomp: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Regions](#regions) - Reference check of all target BEDs and the common comparison region
- [Coverage](#coverage) - Depth of every alignment and per-sample callable regions
- [Preprocessed VCFs](#preprocessed-vcfs) - Reference check, normalization, harmonization and uniform filtering of every VCF
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

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
- `vcf_preprocessed/reference_check/<sample>_<technology>.vcf.reference_check.tsv`: result of
  the reference check per VCF (`status` = `ok`, `fixed` or `error`; renamed contigs, header
  fixes, records on unknown contigs, REF mismatches, sample names).

</details>

Every VCF goes through the same steps:

1. Reference check (contig names and lengths, REF alleles, sample names).
2. Tumor and normal columns selected in this order; recompressed with BGZF.
3. Contigs renamed to the reference (`chr` prefix, `MT`/`chrM`).
4. Records on non-reference contigs and pure gVCF reference blocks removed.
5. Contig header lines rewritten from the reference `.fai`.
6. `bcftools norm`: multiallelic records split, indels left-aligned, MNVs atomized into
   SNVs, records with a REF mismatch dropped, duplicates removed.
7. Harmonization: caller-independent INFO tags computed from FORMAT/AD
   (`VC_T_VAF`, `VC_T_DP`, `VC_T_ALT`, `VC_N_VAF`, `VC_N_DP`, `VC_N_ALT`, `VC_TYPE`,
   `VC_LEN`). Every record is a call of the tumor (GT tumor `0/1`, normal `0/0`), and
   remaining symbolic alleles are dropped.
8. Uniform filters (`--filter_pass_only`, `--min_vaf`, `--min_depth`, `--min_alt_reads`) and
   restriction to the comparison region (records overlapping it are kept).

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
  - Validated samplesheet with resolved paths: `samplesheet_validated.tsv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
