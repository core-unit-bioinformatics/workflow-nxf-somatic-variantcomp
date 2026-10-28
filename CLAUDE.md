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
10. Test data: SEQC2 HCC1395/HCC1395BL, **WES + WGS only for now**. They are
    prepared with **seqc2-toolkit** (`../seqc2-toolkit/`, separate repo, usable
    for any workflow): download, selection, preparation, samplesheet/params
    generation (profile `variantcomp`) and the pipeline run (user decision,
    2026-10-17).
11. If something can't run on this machine, prepare everything (bundle,
    scripts, docs) so the user can transfer it and run it on the HPC.
12. **Reference consistency check + safe auto-fixes.** Check that every input
    (VCF, BED, BAM, known-variant list) is on the configured reference: contig
    names and lengths, build fingerprint, REF alleles.
    - Fix automatically: `chr`/no-`chr` and `MT`/`chrM` naming, missing contig
      lines.
    - Never fix contig lengths: a length mismatch of a nuclear contig is an
      error. chrM is dropped from all inputs and reported (`--keep_mito`
      keeps it, with normal length checks) (user decision, 2026-10-13).
    - Fail hard: build mismatch, or a REF mismatch rate above the threshold.
      No liftover.
13. **Last feature (`feat-known-variants`).** An extra list of known
    variants (TSV or VCF), in two modes:
    - `required`: a clinical checklist of fewer than 10 variants. Each one
      must be found in every technology; report clear PASS/FAIL per
      sample × technology × variant.
    - `benchmark`: an external truth set (SEQC2 v1.2.1, GIAB HG008).
      Recall/precision per technology within its confident regions.
14. **Sheets are CSV** (samplesheet, stratification sheet, known-variant
    sheet). This overrides the CUBI TSV default for this project (user
    decision, 2026-10-08).
15. **No index columns** in the samplesheet. A VCF index (`.tbi`/`.csi`) or
    alignment index (`.bai`/`.csi`/`.crai`) next to the file is used;
    otherwise the pipeline creates it (`bcftools index`, `samtools index`).
16. Review and follow the general guidelines below for all code: nf-core
    modules for tools, logic in scripts, inline comments, strict checks,
    detailed documentation with tested inputs and tool citations.

## General nf-core workflow development guidelines

These guidelines are project-independent and can be reused for other
nf-core-style pipelines.

### Architecture: modules run tools, scripts hold logic

- **One tool per module.** Every process runs one tool (or one script). Prefer,
  in this order:
  1. an nf-core module (`nf-core modules install`);
  2. an nf-core module adapted with `nf-core modules patch` (the patch is
     tracked in the module folder);
  3. a local module that runs one tool command (bcftools, bedtools, samtools,
     mosdepth, rtg-tools, ...) without pipes into other tools;
  4. a local module that only calls a script from `bin/`.
- **Don't pipe tools into each other inside a module** (e.g. `bcftools mpileup |
  bcftools norm`). Split the chain into modules connected by channels, so each
  step can be configured, tested and replaced on its own.
- **Logic lives in scripts, not in module code.** Anything beyond calling a
  tool (decisions, parsing, statistics, merging tables, building report
  files) goes into a script in `bin/`. The module only stages files and calls
  the script. Module `script:` blocks contain no awk programs, loops over
  files or conditional logic beyond assembling arguments.
- **All tool and script options are configurable.** Modules pass
  `task.ext.args` (and `ext.args2`, ...); defaults are set in
  `conf/modules.config`, so users can change them with a config file
  without touching code.
- **Before writing custom code**, check whether an existing tool or plugin
  does the job (e.g. `bcftools norm --check-ref`, `bcftools +fixref`,
  `bcftools +fill-tags`, `bedtools multiinter`). If a custom script is still
  needed, document in its docstring which tools were considered and why
  they are not enough.

### Checks and errors

- **Minimise logic in Nextflow code.** Workflows and subworkflows connect
  channels; they don't validate inputs, parse files or decide about data.
  Validation (samplesheet cross-row checks, file existence, consistency of
  inputs, pass/fail decisions) is done by Python scripts in processes.
- **Be strict.** If a value is needed downstream, check that it exists and is
  valid, and stop with an error if not. An error is better than a wrong
  result. Never silently fall back to a default that changes results; if a
  fallback is intended (e.g. VAF from AF when AD is missing), log it and
  count it in a report.
- **Error messages must help debugging.** They name the file, the record or
  row, the expected and the observed value, and what the user can do
  (e.g. `ERROR: S1_wes.vcf.gz: sample 'TUMOR1' not found (samples in the VCF:
  TUMOR, NORMAL); check tumor_id in the samplesheet`). Scripts exit with a
  non-zero status on errors and print `ERROR:` / `WARNING:` / `INFO:`
  prefixed messages to stderr.
- Scripts log what they did: which fields were used, how many records were
  read, kept, dropped (with the reason) and written.

### Nextflow strict syntax

- **All Nextflow code follows the strict syntax** (Nextflow language v2, the default
  parser since 25.10): `nextflow run`, `nextflow inspect` and therefore
  `nf-core pipelines download` reject code that doesn't compile ("Script compilation
  failed"). Rules: no Groovy outside the allowed statements (no `import`, classes or
  top-level code in scripts; helper functions are `def` in a workflow/utils file); `def`
  for every local variable; lower-case `channel`; no `shell:` blocks or `addParams`;
  closures with explicit parameters (`{ meta, bed -> ... }`, not `it`); no deprecated
  operators. See https://www.nextflow.io/docs/latest/strict-syntax.html.
- **Check with `scripts/check_syntax.py` before every commit that touches `.nf` or
  `.config` files.** It runs `nextflow inspect -profile singularity,test,test_full` (the
  command of `nf-core pipelines download`) with the minimum Nextflow version of the
  manifest and with the current one, and `nextflow lint` of the whole repo. It must
  pass: no inspect failure, no lint error, no lint warning except the accepted ones.
  `nextflow lint` alone is not enough: it reported "no errors" for a stray closing brace
  that broke compilation (commit bcefbcb, fixed in d9973a7).
- **Accepted lint warning:** "The `when` section is deprecated". nf-core module lint
  requires `when: task.ext.when == null || task.ext.when` in every module, so it stays
  until the nf-core template changes. Every other warning is fixed.
- The reference linter is the newest Nextflow of the dev env (`envs/variantcomp-dev.yaml`);
  linters of other versions disagree on style warnings (e.g. single named emits).

### Code style and comments

- **Python, not R.** Scripts live in `bin/`: `#!/usr/bin/env python3`
  shebang, executable, argparse, PEP8: formatted with `black` and linted with
  `ruff` (both configured in `pyproject.toml`; run `black bin scripts tests`
  and `ruff check` before every commit), no workflow-engine logic, unit-tested with pytest
  (`tests/python/`). Plots use matplotlib, statistics scipy/pandas.
- **Inline comments everywhere.** Every block of code executed by the
  pipeline (Python scripts and the shell code of modules) has comments that
  describe what the next line(s) do exactly, not only the overall design.
  Example: `# keep only records whose only ALT is a symbolic allele
  (gVCF reference blocks) for removal`.
- Shell code in modules: one tool call per module, each option on its own
  line, with a comment for non-obvious options.

### Documentation

- `docs/usage.md` documents every input and parameter; `docs/output.md`
  every output file and every column of the main tables.
- **Document the logic, not only the interface.** If a step picks values
  from several possible sources (e.g. VAF from FORMAT/AD, else FORMAT/AF),
  the documentation describes each source, the order of preference, the
  checks and what happens when a value is missing.
- **Tested inputs.** The documentation lists exactly which inputs were tested
  (e.g. caller, caller version, VCF version, test dataset) and which are
  expected to work but untested, and the formal input requirements (VCF
  version, required FORMAT/INFO fields, sample columns, contig naming).
- **Tools and citations.** README and `docs/output.md` name, for every step,
  the tools used. `CITATIONS.md` and the MultiQC methods description cite
  every tool and dataset with its authors and publication (or repository and
  version if there is none). Credit the authors of all tools used.
- The README contains an nf-core-style overview figure (metro map) of the
  pipeline (`docs/images/`), kept up to date with the workflow.

## Project design principles

- **Identical treatment of every technology.** One normalization chain and
  one filter rule set for all VCFs.
- **Normalization chain:** this order is mandatory; don't reorder or skip a
  step without discussing it. Implemented in `subworkflows/local/preprocess_vcf`.
  1. reference check (contigs, lengths, REF alleles, sample names), done for
     all inputs at the start of the run (`subworkflows/local/check_inputs`)
  2. select tumor (+ normal) in this order, recompress to BGZF (inputs may be
     plain gzip)
  2b. `--filter_pass_only` (default): split into PASS (FILTER `PASS`/`.`) and
     non-PASS records; only PASS records are compared, the non-PASS records
     run through the same steps 3-7 (`subworkflows/local/normalize_vcf`) so
     that DISCORDANCE can report them as `filtered_non_pass`
  3. contig renaming
  4. keep reference contigs only; drop records whose only ALT is symbolic
     (gVCF blocks)
  5. rewrite contig header lines from the `.fai`
  6. `-m-any` split + left-align + atomize MNPs + drop REF mismatches +
     deduplicate (`bcftools norm`)
  7. harmonize VAF/DP/alt counts (`VC_*` INFO tags; `--vaf_source`, auto
     prefers FORMAT/AF, else AD), variant type,
     GT (tumor `0/1`, normal `0/0`); drop remaining symbolic alleles
  8. uniform filtering + restriction to the common/callable regions (the
     unfiltered normalized VCF is kept as well)
- **Matching.**
  - Primary: `rtgtools vcfeval --squash-ploidy` (haplotype-aware).
  - Also report exact-allele (`bcftools isec`) and position-only matches, each with its
    own published outputs and metrics (`benchmark/<tool>/`, column `method`), and the
    isec status of every variant. All methods count the same records in the same region.
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
- **Reasons of discordant variants** (D38). Exactly one **main reason** (`class`, first
  entry of `reasons`) from the view of the technology that misses the variant; they exclude
  each other:
  - called, not compared: `filtered_non_pass` (FILTER not PASS), `present_but_filtered`
    (PASS, but removed by the uniform filters/region);
  - not called: `absent` (no alignment), `low_coverage`, `subthreshold_evidence` (alt
    reads), `no_evidence` (no alt reads). Read-evidence reasons never for called variants.
  - **Further reasons** in addition: `germline_in_normal`, and the context `indel`,
    `low_vaf`, `in_repeat`, `giab_<name>`, stratifications.
  - `reasons` is what readers see; `consistency_checks.py` recomputes every rule.
  - There is no `outside_callable` reason: both technologies of a sample are compared in
    the same region (open question in `docs/roadmap.md`).
- **Read evidence** is allele-specific (`read_evidence.py`): indels by haplotype sequence,
  never by a pileup's own allele discovery (it reports other indel alleles than callers).
- **No hardcoded paths, env vars or downloads** in the pipeline. All
  resources come in via params. It must run offline on an HPC (PBS Pro +
  Singularity).
- **Samplesheet:**
  - `.csv` with a `sample` column, `n/a` (or empty) for empty optional
    fields; per-column formats validated with nf-schema, cross-row checks in
    a Python script;
  - columns: `sample`, `technology`, `library_type`, `tumor_id`,
    `normal_id`, `caller`, `vcf`, `tumor_bam`, `normal_bam`, `target_bed`.
- **Configs** are `.yaml` (params files) or `.config` (Nextflow).
  `.toml` is only used for `pyproject.toml`.

### Design essence (feedback 2026-10-13)

- **All input checks first.** `subworkflows/local/check_inputs` checks every
  sheet, VCF, alignment, BED and known-variant list against the reference
  before any other work starts; every other subworkflow waits for it. Each
  check can be switched off (`--check_vcf_contigs`, `--check_ref_alleles`,
  `--check_bam_contigs`, `--check_bed_contigs`, `--check_known_variants`);
  contig renaming is still computed, only the errors become warnings
  (`--check_ref_alleles false` skips the REF comparison).
- **No silent fixes of contig lengths;** chrM is dropped and reported.
- **pysam in Python scripts** for VCF/BAM/CRAM/FASTA; no processes that only
  extract headers or contig lists for a script. Readers (`read_bed`,
  `read_fai`, `open_vcf`, `open_alignment`) are functions in
  `bin/check_reference.py`; other scripts import them from there (no shared
  helper module).
- **Every script** has an INPUT/OUTPUT docstring (columns, types, 0- or
  1-based coordinates) and a docstring for every function; `main()` only
  parses arguments and calls functions.
- **Every threshold** of a script or tool is a documented pipeline param
  (`nextflow.config`, `nextflow_schema.json`, `docs/usage.md`), passed via
  `ext.args`.
- **Every step summary reaches MultiQC** (`bin/mqc_tables.py`, `MQC_TABLES`).
- **One config file per subworkflow** in `conf/modules/`, included from
  `conf/modules.config`.
- **Optional features are params, and require their inputs.**
  `--callable_regions` and `--read_evidence` need BAMs, `--vep` a cache or a
  GTF, `--require_normal_coverage` normal BAMs; missing inputs stop the run in
  the first step with a hint which param to switch off. Without these
  features, BAMs are ignored.
- **Annotation like variantinterpretation:** Ensembl VEP (`+split-vep`) for
  consequence, gene, impact; msisensor-pro for repeats; GIAB genome
  stratifications (one file param per stratification, `--giab_<name>`: a
  given file switches it on, none by default); `--stratification` BEDs as
  extra context. No gene BED.
- **Resources are prepared, never downloaded by the pipeline.**
  `scripts/download_resources.py` gets the VEP cache, GIAB stratifications and
  a GENCODE GTF beforehand and writes `resources.config`. Download code uses
  Python (`urllib`, resumable) with an explicit CA bundle, not `curl`.
- **Statistics from scipy** (`fisher_exact`, `false_discovery_control`),
  not own implementations.

### Report and end-of-run checks (feedback 2026-10-22 and 2026-10-26, D34/D37)

- **MultiQC is the main report** (D37, replaces the narrative report of D33). The variantcomp
  sections (`bin/report_content.py`, ids `vc_*`, order in `assets/multiqc_config.yml`) come
  first, ordered by importance: Summary (cards of recall/precision/F1 per matching method,
  run cards, key findings), Introduction (bullet lists), key results, then details, then
  the QC of the steps, methods and versions. Concise, not overloaded: every section starts
  with text that explains what is shown.
- **Plots for readers.** Bar plots are vertical with the samples on the x-axis (Plotly of
  the MultiQC report via `JS_LIB`, following the toolbox: hide, highlight, rename);
  recall/precision axes always span 0 to 1; sizes in Mbp; matching methods (vcfeval,
  isec exact, isec position) can be switched where metrics are shown.
- **Variant tables live in datavzrd** (nf-core `datavzrd`, `report/tables/`). MultiQC shows
  at most 15 variant rows (`MAX_ROWS`) and links to the table. datavzrd views start with a
  description, show gene + HGVS.p first, colour status columns, link to GeneCards, Ensembl,
  dbSNP and COSMIC, put minor columns in the detail view, and link the tables by
  `variant_id` (`sample:chrom:pos:ref>alt`) and `gene`. The MultiQC report is a datavzrd page.
- **Shareable.** `report/variantcomp_report.zip` (`report_bundle.py`) keeps the output layout
  (`multiqc/`, `report/tables/`), so the relative links between the reports work.
- **One source of explanations.** Class, reason and column texts live once in
  `report_content.py` (`CLASS_TEXT`, `CLASS_LABEL`, `REASON_TEXT`, `COLUMN_TEXT`) and are
  used by both reports; keep them in line with `docs/output.md`.
- **Every count is traceable.** Each VCF has an accounting table (records kept, removed
  with the reason, added per step; `variant_accounting.py`) and a table of every removed
  record; the comparison region has a step table (`region_summary.py steps`); coverage is
  summarized in the common and in the target region (`coverage_summary.py`, MultiQC
  section Coverage).
- **Consistency checks at the end.** `consistency_checks.py check` verifies that the counts
  of all steps add up. New steps that change counts add their checks there. MultiQC
  shows the summary, datavzrd the table; `CHECKS_GATE` stops the run if a check failed, after the
  reports were written (so the failure can be inspected). Never weaken a check to make a
  run pass; find the cause.
- **No versions topic after the versions are collated.** A process that waits for the
  collated versions, directly or through another report (MULTIQC, REPORT_BUNDLE, DATAVZRD
  (patched), CHECKS_GATE), must not emit into the `versions` topic: the topic would stay
  open and the run would hang.

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
  - **T1**: `-profile test_seqc2`, SEQC2 small slices (prepared by
    seqc2-toolkit, `datasets/variantcomp_small.yaml`), stored in the test-data
    repo, runs on the Pi.
  - **T2**: `-profile test_full`, full SEQC2 WES + WGS on the HPC.
- nf-test shares conda environments via `tests/nextflow.config`
  (`NFT_CONDA_CACHE`).
- Run before every commit that touches pipeline code:
  - `scripts/check_syntax.py` (strict syntax: compiles with `nextflow inspect`, lint clean)
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
