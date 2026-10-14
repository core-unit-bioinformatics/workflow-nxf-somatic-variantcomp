# core-unit-bioinformatics/workflow-nxf-somatic-variantcomp: Testing

## Overview

| Tier | What | Data | Where | Command |
| --- | --- | --- | --- | --- |
| T0 | Unit tests of `bin/` and test scripts | built on the fly | anywhere | `pytest` |
| T0 | nf-test function/module/subworkflow tests | synthetic dataset | anywhere | `nf-test test --tag subworkflows_local --profile +conda` |
| T0 | Pipeline on the synthetic edge-case dataset | `synthetic/` | anywhere | `-profile test` |
| T1 | Pipeline on small SEQC2 data (chr22) | `seqc2/small/` | laptop / Raspberry Pi | `-profile test_seqc2` |
| T2 | Pipeline on full SEQC2 data | `seqc2/full/` (prepared on the HPC) | HPC | `-profile test_full` |

All test **data** live in the separate repository `workflow-nxf-somatic-variantcomp-testdata`.
By default, the test profiles expect a checkout of it next to this repository:

```
<workspace>/
├── workflow-nxf-somatic-variantcomp/
└── workflow-nxf-somatic-variantcomp-testdata/
```

Use `--pipelines_testdata_base_path <path-or-URL>/` (or `NFT_TESTDATA` for nf-test)
to point elsewhere. Every dataset there is created by code in this repository and can be
regenerated at any time.

## Development environment

```bash
conda env create -f envs/variantcomp-dev.yaml
conda activate variantcomp-dev
```

It provides Nextflow, nf-core tools, nf-test, Python (pysam, pandas, scipy, matplotlib,
pytest, ruff) and bcftools, samtools, bedtools, mosdepth and rtg-tools. On aarch64 (e.g. a
Raspberry Pi), use `-profile conda`: the biocontainer images are amd64 only.

## Synthetic edge-case dataset

`tests/synthetic/make_synthetic_data.py` writes a small random reference (chr1, chr2,
chrM), DRAGEN-like tumor/normal VCFs for two samples and three technologies (`wgs` =
truth, `wes`, `panel`), target BEDs, a samplesheet and `design.tsv`. `design.tsv` records
the intended outcome of every designed variant. The variants cover:

- an MNV vs. split SNVs
- unnormalized indels in a homopolymer and a dinucleotide repeat
- a multiallelic site
- `<NON_REF>` alleles and a gVCF reference block
- duplicated records
- filtered true variants
- variants outside the common target region

The VCF headers also include problems the pipeline must fix:

- `1`/`2`/`MT` contig names without contig header lines
- a wrong chrM length in the header
- the normal column before the tumor column

```bash
tests/synthetic/make_synthetic_data.py --outdir ../workflow-nxf-somatic-variantcomp-testdata/synthetic
nextflow run . -profile test,conda --outdir results_test
```

## SEQC2 HCC1395/HCC1395BL (WES + WGS)

The SEQC2 Somatic Mutation Working Group sequenced one tumor/normal cell-line pair
(HCC1395 breast cancer, HCC1395BL matched B-lymphoblastoid normal) from a single DNA
batch with WGS and WES (Agilent SureSelect v6 + UTR) at several centers, with
replicates. Fang et al. 2021, *Nat Biotechnol* 39:1151; Xiao et al. 2021,
*Nat Biotechnol* 39:1141; SRA SRP162370.

Each center/replicate pair that has both WGS and WES (e.g. `IL_1` = `WGS_IL_1` +
`WES_IL_1`) is used as one pseudo-sample. The VCFs are the consortium's unfiltered
MuTect2 (GATK 3.8) calls on BWA alignments, with records kept and FILTER set. They use
GRCh38 (`GRCh38.d1.vd1`). The SEQC2 truth set v1.2.1 and its high-confidence regions
are included for the known-variants feature.

`tests/seqc2/seqc2_testdata.py` (configured by `tests/seqc2/seqc2_testdata.yaml`) is
the single entry point:

| Subcommand | Does |
| --- | --- |
| `prepare` | Download (resumable, md5-verified where SEQC2 publishes checksums, cached in `<outdir>/.downloads`), subset, index, and write `samplesheet.csv`, `params.yaml`, `checksums.md5` and `README.md` with full provenance |
| `run` | Run the workflow on a prepared dataset: resolve `params.yaml`, create a run folder, call `nextflow run` with `-resume` |
| `all` | `prepare` + `run` |

| Mode | Content | Size | Default profile |
| --- | --- | --- | --- |
| `small` | chr22 only; pseudo-samples IL_1, FD_1, NV_1; BAM slices (WES targets ±100 bp in chr22:23-24 Mb and 50-50.8 Mb, reads subsampled to 50% with a fixed seed) for IL_1 | ~50 MB (+130 MB download cache) | `test_seqc2,conda` |
| `full` | Genome-wide; 12 pseudo-samples; VCFs only (`--with-bams` adds the full BAMs, 10-160 GB each) | ~3 GB without BAMs (reference 2.9 GB) | `test_full,singularity` |

```bash
# small (this is what the test-data repo contains)
tests/seqc2/seqc2_testdata.py prepare --mode small
tests/seqc2/seqc2_testdata.py run --mode small

# full, on the HPC
tests/seqc2/seqc2_testdata.py prepare --mode full --outdir /scratch/seqc2/full
tests/seqc2/seqc2_testdata.py run --mode full --outdir /scratch/seqc2/full \
    --rundir /scratch/seqc2/run_full --profile test_full,singularity -c my_hpc.config
```

Further `nextflow` arguments can be passed after `--`, e.g.
`... run --mode small -- -with-report report.html`.

Notes:

- Some SEQC2 VCFs are plain gzip, not BGZF. The script streams them with
  `bcftools view -t`, and the pipeline must also accept such input.
- The SEQC2 VCF headers list the 2,779 contigs of `GRCh38.d1.vd1`, while the small
  reference only contains chr22. This deliberately exercises the reference check.
- Downloads are done in Python (`urllib`) with resume and retries. Re-running `prepare`
  reuses cached and finished files.
- HTTPS needs a CA bundle; the script uses the first existing one of `--ca-bundle`,
  `$CURL_CA_BUNDLE`, `$SSL_CERT_FILE`, `$REQUESTS_CA_BUNDLE`, the Python package
  `certifi` and the usual system paths, and exports it as `CURL_CA_BUNDLE` /
  `SSL_CERT_FILE` for samtools/bcftools reading remote files. This fixes
  `curl: (77) error adding trust anchors from file: /home/conda/feedstock_root/...`,
  which occurs when a conda build of curl/htslib keeps its unrelocated build-time CA path.
- The GIAB v3.5 stratifications are subset to the dataset's contigs and written in the
  release layout to `annotation/giab/GRCh38/` (`--giab_dir`, default `--giab_files`).

## nf-test

```bash
nf-test test --profile +conda                         # everything
nf-test test tests/default.nf.test --profile +conda   # synthetic pipeline test
nf-test test tests/seqc2_small.nf.test --profile +conda
nf-test test --tag subworkflows_local --profile +conda
```

On the Raspberry Pi, each function test takes about 35 s (JVM start-up), and the
pipeline tests take about 4 minutes.
