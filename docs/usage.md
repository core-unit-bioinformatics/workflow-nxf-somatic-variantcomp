# core-unit-bioinformatics/workflow-nxf-somatic-variantcomp: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The pipeline compares somatic small-variant calls of the **same tumor/normal pairs**
sequenced with **different technologies** (e.g. WGS, WES, targeted panels). One
technology is the ground truth for all samples (`--truth_technology`). The design is
described in [`pipeline-plan.md`](pipeline-plan.md).

## Samplesheet input

The samplesheet is a **comma-separated** file (`.csv`) with a header row and one row per
sample and technology:

```bash
--input '[path to samplesheet.csv]' --truth_technology wgs
```

```csv title="samplesheet.csv"
sample,technology,library_type,tumor_id,normal_id,caller,vcf,tumor_bam,normal_bam,target_bed
P01,wgs,wgs,P01_T,P01_N,dragen-4.4.4,vcf/P01_wgs.hard-filtered.vcf.gz,bam/P01_T_wgs.bam,bam/P01_N_wgs.bam,n/a
P01,wes,capture,P01_T,P01_N,dragen-4.4.4,vcf/P01_wes.hard-filtered.vcf.gz,bam/P01_T_wes.bam,bam/P01_N_wes.bam,bed/exome_v2.bed
P01,panel,amplicon,P01_T,P01_N,dragen-4.4.4,vcf/P01_panel.hard-filtered.vcf.gz,n/a,n/a,bed/panel_v1.bed
```

| Column | Required | Description |
| --- | --- | --- |
| `sample` | yes | Patient/case identifier, identical for all technologies of the same tumor/normal pair. Allowed characters: `A-Z a-z 0-9 _ . -`. |
| `technology` | yes | Technology label, e.g. `wgs`, `wes`, `panel-oca`. Allowed characters: `a-z 0-9 -`. |
| `library_type` | yes | `wgs`, `capture` or `amplicon`. Controls how coverage is computed (duplicates are kept for amplicon data). |
| `tumor_id` | yes | Name of the tumor sample column in the VCF. |
| `normal_id` | yes | Name of the normal sample column in the VCF; empty or `n/a` for tumor-only calls. |
| `caller` | no | Variant caller and version, e.g. `dragen-4.4.4` (default `unknown`). Reported only. |
| `vcf` | yes | **Unfiltered** somatic VCF (all records, FILTER column set by the caller), e.g. DRAGEN `*.hard-filtered.vcf.gz`. Plain, gzip or bgzip compressed, or BCF. |
| `tumor_bam`, `normal_bam` | no | Tumor and normal alignment (BAM or CRAM). Enable coverage-based callable regions and read-evidence checks. |
| `target_bed` | for `capture`/`amplicon` | Target regions of the technology. Empty or `n/a` = genome-wide (WGS). |

**Index files are not listed.** The pipeline looks for them next to the data file:

| File | Index files that are used (first found) | If none is found |
| --- | --- | --- |
| VCF | `<vcf>.tbi`, `<vcf>.csi` | not needed: the first preprocessing step (`bcftools view`) reads the VCF without index and writes a new bgzipped VCF with index. Plain gzip VCFs can't be indexed anyway. |
| BAM | `<bam>.bai`, `<bam>.csi`, `<name>.bai` (`x.bam` → `x.bai`) | created with `samtools index` |
| CRAM | `<cram>.crai`, `<name>.crai` | created with `samtools index` |

The samplesheet is read as in the nf-core template: nf-schema checks every column
(`assets/schema_input.json`) when the pipeline starts. The rules across rows are then
checked by `bin/check_sheets.py` (module `CHECK_SAMPLESHEET`, first step of the workflow;
all problems are reported at once, nothing else starts before it passed):

- Empty optional fields are left empty or written as `n/a`.
- Relative paths are resolved against the directory of the samplesheet, so a project
  folder with samplesheet, VCFs and BEDs can be moved as a whole. Absolute paths and URLs
  are used as given. Missing files stop the pipeline with the missing path.
- At least two technologies are needed, and every sample needs a row for the truth
  technology and at least one other technology.
- Each `sample`/`technology` combination may occur only once.
- All rows of a technology must use the same `target_bed` and `library_type`.
- `capture` and `amplicon` technologies need a `target_bed`.
- A `normal_bam` needs a `tumor_bam` and a `normal_id`.
- Columns of earlier versions (`vcf_index`, `tumor_bai`, `normal_bai`) are rejected with a
  message; tab-separated sheets are rejected.

The validated samplesheet with resolved paths is written to
`<outdir>/pipeline_info/samplesheet.validated.csv`.

## VCF requirements

Every VCF must meet these requirements; the pipeline checks them and stops with an error
naming the file and record otherwise.

| Requirement | Details | Checked by |
| --- | --- | --- |
| Format | VCF 4.1–4.3 (or BCF), plain, gzip or bgzip compressed; read by HTSlib | `bcftools view` |
| Sample columns | the tumor column (`tumor_id`) and, for paired calls, the normal column (`normal_id`); other columns are removed | `check_reference.py`, `bcftools view -s` |
| Reference | contigs and lengths of the configured reference (safe renaming, see below); REF alleles match the FASTA (≤ `--max_ref_mismatch`) | `check_reference.py`, `bcftools norm --check-ref x` |
| Allele depths | `FORMAT/AD` (Number=R: ref, alt) for the tumor, or `FORMAT/AF` + `FORMAT/DP`; every tumor record needs a VAF (see [Harmonization](#harmonization-of-vaf-depth-and-alt-reads)) | `harmonize_vcf.py` |
| FILTER | the caller's FILTER column; `PASS` (or `.`) = passed; the uniform filters are applied on top | `bcftools view` |
| Unfiltered | all records written by the caller, including filtered ones (e.g. DRAGEN `*.hard-filtered.vcf.gz`, which despite its name contains all calls with FILTER set) | – |
| Somatic | every record is a somatic call of the tumor; the GT of the input is not used | `harmonize_vcf.py` |

gVCF reference blocks (`<NON_REF>`, `<*>` as only ALT) are removed; symbolic structural
variant alleles (`<DEL>`, ...) are dropped and counted. Multiallelic records are split,
MNVs atomized and indels left-aligned before matching.

### Tested inputs

| Input | Caller and version | VCF version | Data | Test |
| --- | --- | --- | --- | --- |
| SEQC2 somatic VCFs (HCC1395 / HCC1395BL, sites IL_1 and FD_1, WES and WGS, chr22 slice) | GATK MuTect2 3.8-0-ge9d806836 | VCFv4.2, FORMAT/AD, FORMAT/AF (Number=1), FORMAT/DP | `tests/seqc2/` (T1 on the Pi; T2 full data on the HPC) | `-profile test_seqc2` |
| SEQC2 high-confidence somatic call set v1.2.1 as external truth (known variants, mode `benchmark`) | – (sites-only VCF) | VCFv4.2 | `tests/seqc2/` | `-profile test_seqc2` |
| Synthetic DRAGEN-like VCFs with edge cases (duplicates, MNVs, multiallelics, gVCF blocks, `1`/`MT` contig names, wrong chrM header length, normal-first sample order) | synthetic, modelled on DRAGEN somatic output | VCFv4.2, FORMAT/AD, FORMAT/AF (Number=A), FORMAT/DP | test-data repo `synthetic/` | `-profile test` |
| Known-variant lists | – | CSV / TSV lists, sites-only VCF | test-data repo `synthetic/known/` | `-profile test` |

**Not tested yet:** real DRAGEN v4.4.4 somatic VCFs (the target use case). They are
expected to work: DRAGEN writes `FORMAT/AD`, `FORMAT/AF` and `FORMAT/DP` per sample, tumor
and normal columns named after the sample IDs, and all calls with FILTER set. Other callers
that write `FORMAT/AD` (e.g. GATK Mutect2 4.x, VarDict, LoFreq with AD) are expected to work
but are untested. Callers without `FORMAT/AD` or `FORMAT/AF` + `FORMAT/DP` (e.g. Strelka2,
which writes tier counts `AU/CU/GU/TU` and `TAR/TIR`) are **not supported**; the pipeline
stops with an error.

## Harmonization of VAF, depth and alt reads

Callers store read counts differently, so the pipeline computes the values used for
filtering, VAF strata, detection power and the tables in one way for every VCF
(`bin/harmonize_vcf.py`, after splitting and left-alignment). The results are written as
INFO tags: `VC_T_VAF`, `VC_T_DP`, `VC_T_ALT`, `VC_T_SRC` for the tumor and `VC_N_*` for the
normal, plus `VC_TYPE` (SNV, MNV, INS, DEL, COMPLEX) and `VC_LEN` (length difference).

For each sample, the first rule that applies is used:

| Rule | Condition | VAF | Depth | Alt reads | `VC_*_SRC` |
| --- | --- | --- | --- | --- | --- |
| 1 | `FORMAT/AD` has a ref and an alt value | `AD[1] / sum(AD)` | `FORMAT/DP` if present, else `sum(AD)` | `AD[1]` | `AD` |
| 2 | no usable AD, but `FORMAT/AF` and `FORMAT/DP` | `AF` | `DP` | `round(AF × DP)` | `AF_DP` |
| 3 | neither | missing | `DP` if present | missing | – |

If AD is `0,0` (no reads, e.g. filtered GATK MuTect2 3.8 calls), the VAF is undefined: it
is left missing (the caller's AF is not used), alt reads are 0, and the record is counted
as `zero_depth`, with a warning. These records don't count against `--max-missing-vaf`.

Why these choices:

- VAF uses `sum(AD)` as denominator, so numerator and denominator count the same reads.
  Callers often count only informative reads in AD (DRAGEN, Mutect2), so `DP` can be larger
  than `sum(AD)`; using `DP` would underestimate the VAF.
- Depth prefers `FORMAT/DP`, because it is the coverage the caller reports and is what the
  `--min_depth` filter is meant to act on.
- The caller's own `FORMAT/AF` is not used when AD is available: callers estimate it
  differently (e.g. Mutect2 reports a model-based estimate), which would make VAFs of
  different technologies incomparable.

Checks:

- Every record must have exactly one ALT (the input is split before); `AD` must have two
  values after splitting; AD and DP must not be negative; AF must be within [0, 1].
- `FORMAT/AF = nan` is treated as missing (GATK MuTect2 3.8 writes it for a normal without
  alt reads); it is counted (`*_af_nan` in the summary). AD is used where present.
- Every tumor record needs a VAF (rule 1 or 2). The allowed fraction of tumor records
  without VAF is set by `ext.args = '--max-missing-vaf 0'` of `HARMONIZE_VCF` in
  `conf/modules.config`. A normal without reads at a site is allowed (counted only).
- Rule 2 is reported as a warning with the number of records.

The counts per rule and sample are written to
`vcf_preprocessed/<sample>/<sample>_<technology>.normalized.harmonize_summary.tsv`.

The tumor GT is set to `0/1` and the normal GT to `0/0`: every record of a somatic VCF is a
tumor call, and callers fill GT differently (splitting a multiallelic record can leave
`0/0` for an allele).

## Reference genome

All inputs must be on the reference given with `--fasta` (plain or bgzipped; a `.fai` can be
given with `--fasta_fai` and is created otherwise). The pipeline checks every input against
it:

| Check | Action |
| --- | --- |
| Contig names differ only by the `chr` prefix or `MT`/`chrM` | renamed automatically |
| Contig on the reference but with a different length in a VCF header (mitochondrial only, e.g. `chrM` 16570 vs 16571) | header rewritten from the `.fai` |
| Contig with a different length (nuclear), BED interval beyond a contig end | **error**: different genome build |
| More than `--max_ref_mismatch` (default 1%) of VCF REF alleles differ from the FASTA | **error**: different build or reference |
| Records/intervals on contigs absent from the reference (decoys, alt contigs, ...) | dropped and reported |

There is no liftover. Convert inputs to one build before running the pipeline.

## Regions

The comparison is restricted to the intersection of the target BEDs of all targeted
technologies. Technologies without a BED (WGS) are genome-wide and don't restrict it.

| Parameter | Description |
| --- | --- |
| `--regions_bed` | Restrict the comparison further, e.g. to high-confidence regions of a truth set. |
| `--exclude_bed` | Remove regions, e.g. a blacklist or problematic regions. |

## Coverage and callable regions

If tumor alignments (BAM/CRAM) are given for **every** technology of a sample, the
comparison for that sample is restricted to the bases that are callable in all of them:

- `mosdepth` computes the depth of every alignment in the common region.
- In `--callable_mode base` (default), a base is callable if its depth is at least
  `--min_coverage` (default 10).
- In `--callable_mode region`, a common-region interval is callable if its mean depth is at
  least `--min_coverage`. This is the approach of the original OCA+/WES comparison; mean
  depth can hide low-coverage parts of an interval.
- With `--require_normal_coverage` (default `true`), the normal alignments must be callable,
  too.
- Reads flagged as duplicates are excluded, except for `library_type` `amplicon`, where all
  reads of an amplicon look like duplicates.
- Alignment headers are checked against the reference like VCFs and BEDs.

Samples without alignments for every technology keep the common region. The comparison
still works for them, but variants in poorly covered parts of a target can't be
distinguished from real discordance.

## Discordance analysis

| Parameter | Default | Description |
| --- | --- | --- |
| `--gene_bed` | – | BED4 with gene names in column 4, e.g. derived from GENCODE. |
| `--stratification` | – | Comma-separated sheet with columns `name` and `bed` (relative paths resolved against the sheet), e.g. GIAB stratifications for homopolymers, tandem repeats, segmental duplications, low mappability. |
| `--germline_vaf` | `0.1` | Normal VAF from which a discordant variant is classified as `germline_in_normal`. |
| `--evidence_min_alt_reads` | `2` | Alt reads in a normal needed for `germline_in_normal`. |
| `--power_min_alt_reads` | `3` | Alt reads a caller needs (power model). |
| `--detection_power` | `0.95` | Power from which a truth variant counts as detectable. |
| `--low_vaf` | `0.1` | VAF threshold of the low-VAF feature in the enrichment tests. |

```csv title="stratification.csv"
name,bed
homopolymer,giab/GRCh38_AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz
tandem_repeat,giab/GRCh38_AllTandemRepeats.bed.gz
segmental_duplication,giab/GRCh38_segdups.bed.gz
```

## Known variants

`--known_variants` takes a comma-separated sheet of variant lists that are checked in
addition to the comparison between technologies:

```csv title="known_variants.csv"
name,sample,mode,variants,regions
clinical,P01,required,known/P01_clinical.csv,n/a
seqc2_truth,all,benchmark,truth/sSNV.vcf.gz;truth/sINDEL.vcf.gz,truth/High-Confidence_Regions.bed
```

| Column | Description |
| --- | --- |
| `name` | Name of the list (unique). |
| `sample` | Samplesheet sample, or `all` for every sample. |
| `mode` | `required` or `benchmark` (see below). |
| `variants` | One or more VCF, CSV or TSV files separated by `;`. CSV/TSV columns: `chrom`, `pos` (1-based), `ref`, `alt`, optional `gene`, `label`, `vaf`. |
| `regions` | Optional BED (`n/a` if none), e.g. the confident regions of a truth set. |

`--stratification` and `--known_variants` are plain file parameters. Each sheet is checked
by `bin/check_sheets.py` in the subworkflow that uses it (`DISCORDANCE`,
`PREPARE_KNOWN_VARIANTS`): required columns, unique names, allowed modes, samples that
exist in the samplesheet, existing files with the right extensions. All problems are
reported at once; the checked sheets are written to `pipeline_info/*.validated.csv`.

Every list is checked against the reference (contig names are renamed if needed; a REF
mismatch is an error) and normalized like the VCFs.

**`required`**: a small list of variants (e.g. fewer than 10 manually checked, clinically
relevant variants) that **every technology must find**. For every technology the pipeline
reports per variant:

| Status | Meaning |
| --- | --- |
| `PASS` | Called and passing the uniform filters in the comparison region. |
| `FAIL_filtered` | Called, but failed FILTER or the uniform thresholds. |
| `FAIL_outside_region` | Not in the comparison region (common targets / callable region). |
| `FAIL_not_called` | In the region, but not called. Read evidence (if alignments are given) shows whether reads support it. |

Every FAIL is logged as a warning and shown in the MultiQC report. With
`--fail_on_missing_required`, the pipeline stops with an error.

**`benchmark`**: an external truth set, e.g. the SEQC2 reference call set or a GIAB
benchmark (such as HG008, if you sequence it yourself). Every technology, including the truth
technology, is compared with it by `rtg vcfeval` within `regions` ∩ the comparison region,
with the same metrics as the main benchmark. If the list contains a VAF (VCF `INFO/TVAF`,
`VAF` or `AF`, or the TSV column `vaf`), recall is also stratified by it.

## Filters

The input VCFs are **unfiltered** (all calls with the caller's FILTER column). The pipeline
applies one set of filters identically to every technology. VAF, depth and alt read counts
are computed from FORMAT/AD the same way for every caller.

| Parameter | Default | Description |
| --- | --- | --- |
| `--filter_pass_only` | `true` | Only calls with FILTER `PASS` (or `.`) are compared. |
| `--min_vaf` | `0` (off) | Minimum tumor VAF. |
| `--min_depth` | `0` (off) | Minimum tumor depth. |
| `--min_alt_reads` | `0` (off) | Minimum tumor alt read count. |

Filtered calls are kept in the normalized VCFs, so a missed variant can later be classified
as "present but filtered".

## Comparison options

| Parameter | Default | Description |
| --- | --- | --- |
| `--truth_technology` | required | Technology used as ground truth for all samples. |
| `--comparison_mode` | `truth` | `truth`: compare every technology with the truth technology. `all_pairs`: additionally compare all technologies pairwise. |
| `--vaf_bins` | `0,0.05,0.1,0.2,0.5,1` | VAF bin edges for the VAF-stratified metrics. |

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run core-unit-bioinformatics/workflow-nxf-somatic-variantcomp --input ./samplesheet.csv --truth_technology wgs --outdir ./results -profile singularity
```

This will launch the pipeline with the `singularity` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run core-unit-bioinformatics/workflow-nxf-somatic-variantcomp -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
truth_technology: 'wgs'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull core-unit-bioinformatics/workflow-nxf-somatic-variantcomp
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [core-unit-bioinformatics/workflow-nxf-somatic-variantcomp releases page](https://github.com/core-unit-bioinformatics/workflow-nxf-somatic-variantcomp/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
