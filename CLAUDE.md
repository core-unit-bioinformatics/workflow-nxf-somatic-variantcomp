# CLAUDE.md: workflow-nxf-somatic-variantcomp

The global CUBI rules in `~/.claude/CLAUDE.md` always apply: naming, semver,
git process, ISO dates, sample sheets, file types. This file adds the
project-specific instructions. The full design, roadmap and agreed decisions
are in `docs/pipeline-plan.md`. Read the plan before starting any feature,
and keep it up to date when decisions change (§11).

## Role

Act as a bioinformatician and Nextflow / nf-core developer.

## Names

- Repo / manifest name: `core-unit-bioinformatics/workflow-nxf-somatic-variantcomp`.
- Internal nf-core template name: `variantcomp`. Nextflow identifiers can't
  contain dashes, so the main workflow is `VARIANTCOMP` and the local utils
  subworkflow is `utils_nfcore_variantcomp_pipeline`.
- Test-data repo: `../workflow-nxf-somatic-variantcomp-testdata/`, a separate
  git repo. It holds data only. All test *code* lives in this repo under
  `tests/`.

## What we build

An nf-core-template Nextflow DSL2 pipeline that compares somatic SNV/indel
VCFs of the **same tumor/normal pair** sequenced with **different
technologies** (WGS, WES, targeted panels; Illumina short reads). One
technology is the ground truth for all samples (`params.truth_technology`).

The main outputs are:

- recall / precision / F1;
- tables and statistics of discordant variants: type, indel length, VAF,
  gene, repeat/homopolymer context, coverage, filter status, read evidence.

This is **not** nf-core/variantbenchmarking (different callers, same data).
Here the *data* differ, so normalization, region restriction and
VAF/detection-power awareness are essential.

## User requirements (do not drop any)

1. ≥2 technologies and ≥1 sample. The real use case is 8 samples × 3
   technologies, called with **DRAGEN v4.4.4**, all **tumor/normal paired**.
   DRAGEN docs: help.dragen.illumina.com.
2. Input: **unfiltered** VCFs. BAMs are optional and enable coverage and
   evidence checks. FASTQ is out of scope for v1.
3. Ground truth: one technology for all samples, set in the config.
4. Recall / precision statistics as the main output.
5. Tables of variants that differ between technologies, with statistics:
   indel vs. SNV, gene clustering, repeat regions, etc.
6. Target BEDs are input. The comparison is restricted to the **intersection
   of all BEDs** (and, with BAMs, to regions callable in all technologies).
7. Design inspiration: `cio-abcd/variantinterpretation`. Use its nf-schema
   validation, VCF checks, split/left-align, datavzrd reports and MultiQC.
8. git from the beginning. Version every feature and step properly.
9. Dev/test tooling goes into the conda env `envs/variantcomp-dev.yaml`
   (name `variantcomp-dev`).
10. Test data: SEQC2 HCC1395/HCC1395BL, **WES + WGS only for now**. A
    reusable, reproducible script (`tests/seqc2/`) handles download,
    preparation, samplesheet/params generation and the pipeline run.
11. If something can't run on this machine, prepare everything (bundle,
    scripts, docs) so the user can transfer it and run it on the HPC.
12. **Reference consistency check + safe auto-fixes.** Check that every input
    (VCF, BED, BAM) is on the configured reference: contig names and
    lengths, build fingerprint, REF alleles.
    - Fix automatically: contig header lengths (e.g. chrM 16570 vs 16571),
      `chr`/no-`chr` and `MT`/`chrM` naming, missing contig lines.
    - Fail hard: build mismatch, or a REF mismatch rate above the threshold.
      No liftover.
13. **Last feature (`feat-known-variants`).** An extra list of known
    variants (TSV or VCF), in two modes:
    - `required`: a clinical checklist of fewer than 10 variants. Each one
      must be found in every technology; report clear PASS/FAIL per
      sample × technology × variant.
    - `benchmark`: an external truth set (SEQC2 v1.2.1, GIAB HG008).
      Recall/precision per technology within its confident regions.

## Design principles

- **Existing tools first.** Order of preference:
  1. an nf-core module (`nf-core modules install`);
  2. bcftools / bedtools / samtools / mosdepth / rtg-tools commands in a
     local module;
  3. a small custom script.

  Before writing custom code, check whether a bcftools command or plugin
  (`norm`, `isec`, `+fill-tags`, `+setGT`, `annotate`, `query`) already does
  the job.
- **Python, not R.** Custom scripts live in `bin/`:
  - `#!/usr/bin/env python3` shebang, executable, argparse;
  - PEP8 (ruff, configured in `pyproject.toml`);
  - no workflow-engine logic;
  - unit-tested with pytest (`tests/python/`).

  Plots use matplotlib (+ upsetplot), statistics use scipy/pandas.
- **Identical treatment of every technology.** One normalization chain and
  one filter rule set for all VCFs.
- **Normalization chain:** this order is mandatory; don't reorder or skip a
  step without discussing it.
  1. reference check + reheader contigs from the `.fai`
  2. contig renaming if needed
  3. select tumor/normal columns
  4. drop gVCF/reference records (`<*>`, `<NON_REF>`)
  5. REF check + `-m-any` split + left-align
  6. atomize MNPs
  7. deduplicate
  8. harmonize VAF/DP/AD tags and GT (caller profiles: `dragen`, `mutect2`)
  9. uniform filtering (keep the unfiltered normalized VCF as well)
  10. restrict to the common/callable regions
- **Matching.**
  - Primary: `rtgtools vcfeval --squash-ploidy` (haplotype-aware).
  - Also report exact-allele (`bcftools isec`) and position-only matches.
- **Regions.**
  - N-way intersection of all target BEDs, with the optional
    `regions_bed` / `exclude_bed`.
  - Callable regions: mosdepth on the tumor (+ normal) of *every*
    technology. Base mode (`--quantize`) is the default; region mode is
    optional.
  - Amplicon libraries: don't filter duplicate-flagged reads.
- **Metrics.** Report overall, by SNV/indel, by VAF bin, by stratum, and as
  power-adjusted recall. Aggregate micro (summed counts) and macro
  (mean ± sd) over samples.
- **Discordance classes.**
  - `absent`
  - `present_but_filtered`
  - `outside_callable`
  - `no_evidence`
  - `subthreshold_evidence`
  - `germline_in_normal`
  - `low_coverage`
- **No hardcoded paths, env vars or downloads** in the pipeline. All
  resources come in via params. It must run offline on an HPC (PBS Pro +
  Singularity).
- **Samplesheet:**
  - `.tsv` with a `sample` column, `n/a` for empty fields, validated with
    nf-schema;
  - columns: `sample`, `technology`, `library_type`, `tumor_id`,
    `normal_id`, `caller`, `vcf`, `vcf_index`, `tumor_bam`, `tumor_bai`,
    `normal_bam`, `normal_bai`, `target_bed`.
- **Configs** are `.yaml` (params files) or `.config` (Nextflow). No `.csv`
  for new files. `.toml` is only used for `pyproject.toml`.

## Git and versioning

- Branches:
  - `main` holds releases only and is fast-forwarded from `dev`;
  - `dev` is the integration branch;
  - `feat-<name>` branches come off `dev` and are rebased into `dev` by
    fast-forward (`git merge --ff-only`), keeping the history linear. No
    merge or squash commits. Delete the feature branch after the merge.
  - `TEMPLATE` is the nf-core template sync branch. Don't commit to it
    manually.
- Every feature includes code, tests, docs (`docs/usage.md`,
  `docs/output.md`) and a `CHANGELOG.md` entry.
- **A version bump is its own commit.** It touches only the version strings:
  `nextflow.config` `manifest.version`, `pyproject.toml`, `.nf-core.yml`, and
  the CHANGELOG header. Tag it `vX.Y.Z-dev` on `dev`. Version strings use the
  CUBI form `0.2.0-dev`. The roadmap is in `docs/pipeline-plan.md` §9.
- Commit locally only. Don't add remotes or push unless the user asks.
- Small, focused commits with imperative messages.

## Environment and testing

- This machine is a Raspberry Pi (aarch64, 4 cores, ~3 GB RAM, ~29 GB free
  disk) with no Docker. Singularity/biocontainer images are amd64, so local
  runs use `-profile conda` with tiny data only. Use low thread and memory
  settings.
- Activate tools with `conda run -n variantcomp-dev <cmd>` or
  `conda activate variantcomp-dev`.
- Test tiers:
  - **T0**: pytest for `bin/` + module/subworkflow nf-tests on synthetic
    edge-case data.
  - **T1**: `-profile test`, SEQC2 small slices (prepared by
    `tests/seqc2/`), stored in the test-data repo, runs on the Pi.
  - **T2**: `-profile test_full`, full SEQC2 WES + WGS on the HPC.
- Run before every commit that touches pipeline code:
  - `nf-core pipelines lint`
  - the relevant `nf-test test`
  - `pytest`
  - a `-stub` run when real data can't run locally
- Report test results faithfully. If something couldn't be run on the Pi,
  say so and give the exact HPC command instead.

## Working agreement

- Discuss changes to the plan or design with the user before implementing
  them. Record agreed decisions in `docs/pipeline-plan.md` §11.
- At each version milestone, summarize what was done, what was tested
  (where), and what is open.
