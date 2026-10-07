# core-unit-bioinformatics/workflow-nxf-somatic-variantcomp

> [!WARNING]
> Under active development (`v0.1.0-dev`). Not ready for production use. See
> [`docs/pipeline-plan.md`](docs/pipeline-plan.md) for the design and roadmap.

## Introduction

**workflow-nxf-somatic-variantcomp** compares somatic small-variant calls (SNVs and indels)
of the **same tumor/normal pairs** sequenced with **different technologies**
(e.g. WGS, WES, targeted panels). One technology is the ground truth for all samples.
The pipeline reports recall, precision and F1, and explains *why* variants differ
between technologies.

Unlike benchmarking pipelines that compare different callers on the same data
(e.g. nf-core/variantbenchmarking), here the input data differ. Meaningful numbers
therefore require:

1. Checking that every input is on the same reference genome, and harmonizing it
2. Identical normalization of all VCFs (split, left-align, atomize, deduplicate, tag
   harmonization) and identical filtering
3. Restriction to the intersection of all target regions and, if BAMs are given, to
   regions callable in all technologies
4. VAF-stratified, detection-power-aware metrics
5. Classification and annotation of every discordant variant (type, gene, repeat
   context, coverage, filter status, read evidence)

Planned workflow steps (see the plan for details):

1. Input and reference consistency checks
2. Region harmonization and N-way intersection of target BEDs
3. Callable regions from coverage (`mosdepth`, optional BAMs)
4. VCF normalization and harmonization (`bcftools`)
5. Comparison (`rtgtools vcfeval`, `bcftools isec`) and metrics
6. Discordance classification and annotation
7. Statistics and reports (`MultiQC`, `datavzrd`)
8. Optional known-variant checks (clinical must-find list, external truth sets)

## Usage

The samplesheet is a tab-separated file with one row per sample and technology.
Its full specification will follow in `docs/usage.md`.

```bash
nextflow run core-unit-bioinformatics/workflow-nxf-somatic-variantcomp \
   -profile <conda/singularity/.../institute> \
   --input samplesheet.tsv \
   --fasta genome.fa \
   --truth_technology wgs \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or the Nextflow `-params-file` option.
> Custom config files, including those provided by the `-c` Nextflow option, can be used
> to provide any configuration _**except for parameters**_; see the
> [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Development

- Development environment: `conda env create -f envs/variantcomp-dev.yaml`
- Branches: `main` (releases), `dev` (development), `feat-*` (features, branched from
  `dev`). The history is kept linear.
- Test data (SEQC2 HCC1395/HCC1395BL, WES + WGS) live in a separate repository and
  are prepared reproducibly with the scripts in `tests/seqc2/`.

## Credits

core-unit-bioinformatics/workflow-nxf-somatic-variantcomp was originally written by Kai Horny.

## Contributions and Support

If you would like to contribute to this pipeline, please see the
[contributing guidelines](docs/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the
[`CITATIONS.md`](CITATIONS.md) file. A machine-readable citation is in
[`CITATION.cff`](CITATION.cff).

This pipeline uses code and infrastructure developed and maintained by the
[nf-core](https://nf-co.re) community, reused here under the
[MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
