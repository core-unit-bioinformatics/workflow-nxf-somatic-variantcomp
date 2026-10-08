# core-unit-bioinformatics/workflow-nxf-somatic-variantcomp: Running on an HPC

The pipeline runs offline on an HPC with Singularity/Apptainer. Everything it needs comes
in through parameters; nothing is downloaded at runtime.

## 1. Build an offline bundle (machine with internet access)

```bash
git clone <repo-url> workflow-nxf-somatic-variantcomp
cd workflow-nxf-somatic-variantcomp
conda env create -f envs/variantcomp-dev.yaml && conda activate variantcomp-dev
scripts/make_offline_bundle.py --tag v1.0.0 --outdir /transfer/variantcomp_bundle --images
```

The bundle contains the pipeline code at the tag, all Singularity images (named like
Nextflow's Singularity cache), the Nextflow plugins and a `README.md`. Without `--images`,
only `containers.txt` lists the images, for example to pull them on the HPC itself.

Once the repository is on GitHub, `nf-core pipelines download` works as well.

## 2. Install on the HPC

```bash
cd /path/to/pipelines
tar -xzf workflow-nxf-somatic-variantcomp-v1.0.0.tar.gz
export NXF_OFFLINE=true
export NXF_PLUGINS_DIR=/path/to/bundle/nxf_plugins
export NXF_SINGULARITY_CACHEDIR=/path/to/bundle/singularity
```

## 3. HPC configuration

Put the scheduler and resources in your own config file, passed with `-c`. Pass
parameters with `-params-file`. Example for PBS Pro:

```groovy title="hpc.config"
process {
    executor       = 'pbspro'
    queue          = '<queue>'
    resourceLimits = [ cpus: 8, memory: 40.GB, time: 72.h ]
    errorStrategy  = { task.exitStatus in (1..200) ? 'retry' : 'finish' }
    maxRetries     = 2
}
executor {
    queueSize = 200
}
singularity {
    enabled    = true
    autoMounts = true
    // fail fast instead of waiting for a download that cannot happen offline
    pullTimeout = '1m'
}
```

```yaml title="params.yaml"
input: /path/to/project/samplesheet.csv
fasta: /path/to/reference/hg19.fa
truth_technology: wgs
gene_bed: /path/to/annotation/genes.bed
stratification: /path/to/annotation/stratification.csv
known_variants: /path/to/project/known_variants.csv
outdir: /path/to/project/results/variantcomp
```

```bash
nextflow run /path/to/pipelines/workflow-nxf-somatic-variantcomp \
    -profile singularity \
    -c hpc.config \
    -params-file params.yaml \
    -resume
```

Limit the memory of the Nextflow head job if the login/head node is shared, e.g.
`export NXF_OPTS="-Xms500M -Xmx8G"`.

## 4. Full SEQC2 test (T2)

The full SEQC2 HCC1395/HCC1395BL WES + WGS test uses 12 pseudo-samples genome-wide:

```bash
# on a machine with internet access (about 3 GB without BAMs)
tests/seqc2/seqc2_testdata.py prepare --mode full --outdir /transfer/seqc2_full
# on the HPC
tests/seqc2/seqc2_testdata.py run --mode full --outdir /data/seqc2_full \
    --rundir /data/seqc2_full_run --profile test_full,singularity -c hpc.config
```

`--with-bams` also downloads the full BAMs (10-160 GB each), which enables coverage-based
callable regions and read evidence. Expected results: WGS against the SEQC2 truth set
(`known_variants/benchmark/known_benchmark.aggregated.tsv`) should be close to the
published values for MuTect2 on BWA alignments (Xiao et al. 2021).
