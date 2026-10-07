# core-unit-bioinformatics/workflow-nxf-somatic-variantcomp: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The pipeline compares somatic small-variant calls of the **same tumor/normal pairs**
sequenced with **different technologies** (e.g. WGS, WES, targeted panels). One
technology is the ground truth for all samples (`--truth_technology`). The design is
described in [`pipeline-plan.md`](pipeline-plan.md).

## Samplesheet input

The samplesheet is a **tab-separated** file (`.tsv`) with a header row and one row per
sample and technology:

```bash
--input '[path to samplesheet.tsv]' --truth_technology wgs
```

```tsv title="samplesheet.tsv"
sample	technology	library_type	tumor_id	normal_id	caller	vcf	vcf_index	tumor_bam	tumor_bai	normal_bam	normal_bai	target_bed
P01	wgs	wgs	P01_T	P01_N	dragen-4.4.4	vcf/P01_wgs.hard-filtered.vcf.gz	n/a	n/a	n/a	n/a	n/a	n/a
P01	wes	capture	P01_T	P01_N	dragen-4.4.4	vcf/P01_wes.hard-filtered.vcf.gz	n/a	n/a	n/a	n/a	n/a	bed/exome_v2.bed
P01	panel	amplicon	P01_T	P01_N	dragen-4.4.4	vcf/P01_panel.hard-filtered.vcf.gz	n/a	n/a	n/a	n/a	n/a	bed/panel_v1.bed
```

| Column | Required | Description |
| --- | --- | --- |
| `sample` | yes | Patient/case identifier, identical for all technologies of the same tumor/normal pair. Allowed characters: `A-Z a-z 0-9 _ . -`. |
| `technology` | yes | Technology label, e.g. `wgs`, `wes`, `panel-oca`. Allowed characters: `a-z 0-9 -`. |
| `library_type` | yes | `wgs`, `capture` or `amplicon`. Controls how coverage is computed later on (duplicates are kept for amplicon data). |
| `tumor_id` | yes | Name of the tumor sample column in the VCF. |
| `normal_id` | yes | Name of the normal sample column in the VCF, `n/a` for tumor-only calls. |
| `caller` | no | Variant caller and version, e.g. `dragen-4.4.4` (default `unknown`). Reported only. |
| `vcf` | yes | **Unfiltered** somatic VCF (all records, FILTER column set by the caller), e.g. DRAGEN `*.hard-filtered.vcf.gz`. |
| `vcf_index` | no | `.tbi`/`.csi` index of the VCF. |
| `tumor_bam`, `tumor_bai` | no | Tumor alignment (BAM/CRAM) and index. Enables coverage-based callable regions and read-evidence checks. |
| `normal_bam`, `normal_bai` | no | Normal alignment and index. |
| `target_bed` | for `capture`/`amplicon` | Target regions of the technology. `n/a` = genome-wide (WGS). |

Rules:

- Empty optional fields are written as `n/a`.
- Relative paths are resolved against the directory of the samplesheet, so a
  project folder with samplesheet, VCFs and BEDs can be moved as a whole. Absolute
  paths and URLs are used as given.
- At least two technologies are needed, and every sample needs a row for the truth
  technology and at least one other technology.
- Each `sample`/`technology` combination may occur only once.
- All rows of a technology must use the same `target_bed` and `library_type`.
- All checks run before any process starts. All problems are reported at once.

The validated samplesheet with resolved paths is written to
`<outdir>/pipeline_info/samplesheet_validated.tsv`.

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
| `--stratification` | – | Tab-separated sheet with columns `name` and `bed` (relative paths resolved against the sheet), e.g. GIAB stratifications for homopolymers, tandem repeats, segmental duplications, low mappability. |
| `--germline_vaf` | `0.1` | Normal VAF from which a discordant variant is classified as `germline_in_normal`. |
| `--evidence_min_alt_reads` | `2` | Alt reads in a normal needed for `germline_in_normal`. |
| `--power_min_alt_reads` | `3` | Alt reads a caller needs (power model). |
| `--detection_power` | `0.95` | Power from which a truth variant counts as detectable. |
| `--low_vaf` | `0.1` | VAF threshold of the low-VAF feature in the enrichment tests. |

```tsv title="stratification.tsv"
name	bed
homopolymer	giab/GRCh38_AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz
tandem_repeat	giab/GRCh38_AllTandemRepeats.bed.gz
segmental_duplication	giab/GRCh38_segdups.bed.gz
```

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
nextflow run core-unit-bioinformatics/workflow-nxf-somatic-variantcomp --input ./samplesheet.tsv --truth_technology wgs --outdir ./results -profile singularity
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
input: './samplesheet.tsv'
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
