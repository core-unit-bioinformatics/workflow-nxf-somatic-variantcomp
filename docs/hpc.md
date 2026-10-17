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

Download the annotation resources on the same machine (VEP cache, GIAB genome
stratifications; `--gtf` adds a GENCODE GTF for VEP without a cache):

```bash
scripts/download_resources.py --outdir /transfer/variantcomp_resources --genome GRCh37
```

It writes `resources.config` (a `params { }` block with the paths) and a `README.md` with
the sources and checksums. The GRCh38 VEP cache v116 is about 27 GB. Copy the folder to
the HPC to the same path, or edit the paths in `resources.config`.

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
known_variants: /path/to/project/known_variants.csv
outdir: /path/to/project/results/variantcomp
```

```bash
nextflow run /path/to/pipelines/workflow-nxf-somatic-variantcomp \
    -profile singularity \
    -c hpc.config \
    -c /path/to/variantcomp_resources/resources.config \
    -params-file params.yaml \
    -resume
```

`resources.config` sets `vep_cache`, `vep_cache_version`, `vep_genome` and one
`giab_<name>` per GIAB file (section 1); the pipeline never downloads them. Without a cache, `vep_gtf`
(bgzipped, tabix-indexed GTF, `download_resources.py --no-vep --gtf`) gives consequences
and genes without the cache-only annotations. The repeat BED of
msisensor-pro (`annotation/repeats/`) can be reused in later runs with `repeats_bed`.

Limit the memory of the Nextflow head job if the login/head node is shared, e.g.
`export NXF_OPTS="-Xms500M -Xmx8G"`.

## 4. Full SEQC2 test (T2)

The full SEQC2 HCC1395/HCC1395BL WES + WGS test uses the first replicate of every center
genome-wide. The data are prepared with seqc2-toolkit (separate repo):

```bash
# on a machine with internet access (about 3.5 GB without BAMs)
seqc2 fetch datasets/variantcomp_full.yaml --outdir /transfer/seqc2/full
seqc2 samplesheet datasets/variantcomp_full.yaml --outdir /transfer/seqc2/full --profile variantcomp
# on the HPC (the folder copied to /data/seqc2/full)
nextflow run /path/to/pipelines/workflow-nxf-somatic-variantcomp \
    -profile test_full,singularity -c hpc.config \
    --pipelines_testdata_base_path /data/ --outdir /data/seqc2_full_run/results -resume
```

`data.bam: true` in the dataset config also downloads the full BAMs (10-160 GB each), which
enables coverage-based callable regions and read evidence (set both params to true). Expected results: WGS against the SEQC2 truth set
(`known_variants/benchmark/known_benchmark.aggregated.tsv`) should be close to the
published values for MuTect2 on BWA alignments (Xiao et al. 2021).
