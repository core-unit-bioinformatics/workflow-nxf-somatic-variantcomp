# core-unit-bioinformatics/workflow-nxf-somatic-variantcomp: Testing

## Overview

| Tier | What | Data | Where | Command |
| --- | --- | --- | --- | --- |
| T0 | Unit tests of `bin/` and test scripts | built on the fly | anywhere | `pytest` |
| T0 | nf-test function/module/subworkflow tests | synthetic dataset | anywhere | `nf-test test --tag subworkflows_local --profile +conda` |
| T0 | Pipeline on the synthetic edge-case dataset | `synthetic/` | anywhere | `-profile test` |
| T1 | Pipeline on small SEQC2 data (chr22) | `seqc2/small/` | laptop / Raspberry Pi | `-profile test_seqc2` |
| T2 | Pipeline on full SEQC2 data | `seqc2/full/` (prepared with seqc2-toolkit) | HPC | `-profile test_full` |

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

The data are prepared with **seqc2-toolkit** (`../seqc2-toolkit/`, a separate repo used
by several workflows): it selects the SEQC2 files, downloads only what is missing into a
shared md5-checked cache, builds the dataset folder (subsets, indices, BAM slices) and
writes the vcf-comp samplesheet, params and known-variant sheet with its `variantcomp`
profile. Install it into the dev env with `pip install -e ../seqc2-toolkit`.

| Dataset config (seqc2-toolkit) | Content | Size | Profile here |
| --- | --- | --- | --- |
| `datasets/variantcomp_small.yaml` | chr22; pairs IL_1, FD_1, NV_1; WGS + WES; bwa + MuTect2; BAM slices (WES targets ±100 bp in chr22:23-24 Mb and 50-50.8 Mb, 50% of the reads, fixed seed) for every pair; truth v1.2.1, GENCODE 46 GTF, GIAB v3.5 | ~150 MB | `test_seqc2` |
| `datasets/variantcomp_full.yaml` | Genome-wide; first replicate of every center with WGS and WES; VCFs only (`data.bam: true` adds the full BAMs, 10-160 GB each) | ~3.5 GB (reference 3.2 GB) | `test_full` |

```bash
cd ../seqc2-toolkit
# small: this is what the test-data repo contains (seqc2/small/)
seqc2 fetch datasets/variantcomp_small.yaml \
    --outdir ../workflow-nxf-somatic-variantcomp-testdata/seqc2/small
seqc2 samplesheet datasets/variantcomp_small.yaml --profile variantcomp \
    --outdir ../workflow-nxf-somatic-variantcomp-testdata/seqc2/small
# run the pipeline on it (or: nextflow run . -profile test_seqc2,conda)
seqc2 run datasets/variantcomp_small.yaml --profile variantcomp \
    --outdir ../workflow-nxf-somatic-variantcomp-testdata/seqc2/small \
    --pipeline ../workflow-nxf-somatic-variantcomp --nf-profile conda
```

The dataset folder holds `dataset.json`, `checksums.md5`, `README.md` (sources, settings)
and the toolkit's manifest; the workflow files are in `seqc2/small/variantcomp/`
(`samplesheet.csv`, `params.yaml`, `known_variants.csv`, `known/`). `seqc2 fetch` again
does nothing if all files are up to date; `seqc2 verify` checks all md5 sums.

Notes:

- Some SEQC2 VCFs are plain gzip, not BGZF; the toolkit recompresses them with
  `bcftools view` (the pipeline also accepts plain gzip input).
- The SEQC2 VCF headers list the 2,779 contigs of `GRCh38.d1.vd1`, while the small
  reference only contains chr22. This deliberately exercises the reference check.
- The toolkit downloads in Python with an explicit CA bundle (no `curl`), which avoids
  `curl: (77) error adding trust anchors` with conda builds of curl.

## nf-test

```bash
nf-test test --profile +conda                         # everything
nf-test test tests/default.nf.test --profile +conda   # synthetic pipeline test
nf-test test tests/seqc2_small.nf.test --profile +conda
nf-test test --tag subworkflows_local --profile +conda
```

On the Raspberry Pi, each nf-test takes at least about 35 s (JVM start-up of nf-test and
Nextflow). Two settings keep the rest short:

- `tests/nextflow.config` activates the conda environments with micromamba (in the dev
  env) instead of `conda info` + `source activate`, which cost about 7 s per task on the
  Pi (`NFT_MICROMAMBA=false` switches it off).
- The test profiles request small resources per label (1 CPU / 768 MB, 2 CPUs / 1 GB):
  the default labels ask for 6 GB, which the 3 GB limit caps to the whole machine, so
  only one task would run at a time.

Measured on the Pi (2026-10-14): the `PREPROCESS_VCF` test went from 847 s to 87 s.
