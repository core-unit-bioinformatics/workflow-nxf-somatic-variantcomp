# Plan: implement feedback `20261013_feedback_claude.md` (v0.14.0-dev → v0.15.0-dev)

## Context

The user reviewed v0.14.0-dev and listed 8 general and 17 specific changes. The
main pain points:

- **Late failures.** Reference checks can stop a run hours in, e.g. overnight.
- **Hidden logic.** Hard-coded thresholds and lost summaries make results hard to trace.
- **Annotation is too basic.** It uses a gene BED instead of VEP and repeat calls.
- **Readability.** Some code is hard to read: long `main()` functions, a big
  `modules.config`, too few comments, column formats not documented.

The feedback also sets the working rules:

- turn the items into a prioritized TODO list;
- work through it step by step, **one small git commit per step**;
- run the **long tests only at the end**;
- save the essence of the design choices in `CLAUDE.md`;
- write a separate **roadmap file** of further improvements.

## Decisions already made (user, 2026-10-13)

| Topic | Decision |
| --- | --- |
| PASS filter | Early: PASS/`.` records only, right after sample selection (step 2), param `filter_pass_only` (default true). The non-PASS records go into a side branch with the same normalization; DISCORDANCE uses it to classify a missed variant as `filtered_non_pass` with its FILTER value. This changes the mandatory normalization chain → update CLAUDE.md and record decision D23. |
| chrM | Always dropped from all inputs and reported (records/intervals); its length is not checked. `--keep_mito` (default false) keeps it, with normal length checks. Contig **lengths are never fixed**: any nuclear length mismatch is an error. |
| Old annotation | `--gene_bed` is removed (genes come from VEP). `--stratification` BEDs stay as an extra, optional context. |
| Cross-row sheet checks | Unchanged (D22): `CHECK_SAMPLESHEET` at workflow start. |

Interpretations taken directly from the feedback (no further question needed):

- **VAF source:** `--vaf_source auto|af|ad`. `auto` (the default) prefers FORMAT/AF (+ DP), else AD; this reverses D21.
- **Missing VAFs:** above `--max_missing_vaf` → a warning, and the VAF strata and detection power of that VCF become NA (no error).
- **Coverage/evidence toggles:** `--callable_regions` and `--read_evidence` (default true) need BAMs; without BAMs the run stops early with a hint to switch them off.
- **Check toggles:** each upfront check can be switched off; contig renaming is still computed, only the errors are skipped.
- **VEP** like variantinterpretation (`vep` default true, cache params, `--flag_pick_allele_gene`, `--terms SO`, ...). In addition, `--vep_gff` runs VEP offline without a cache (used by the tests on the Pi). If VEP is on and neither cache nor GFF is given, the run stops in the input check.

## TODO list (priority order; each bullet = one commit)

Work on one feature branch per phase off `dev` (`feat-<phase>`), fast-forward
merge into `dev`. Before every commit: `black --check`, `ruff`, `pytest`,
`nextflow lint`. nf-tests and pipeline runs only in phase H.

### Phase A — foundations (`feat-code-standards`)

- **A1** Add `black` to `envs/variantcomp-dev.yaml`; configure `[tool.black]`
  (line-length 100) and ruff (lint only, no format rules) in `pyproject.toml`. In CLAUDE.md:
  "Python is formatted with black, linted with ruff".
- **A2** Run `black` on `bin/`, `scripts/` and `tests/`, and remove the now-unneeded
  `# fmt: skip` markers (formatting only).
- **A3** Split `conf/modules.config` into `conf/modules/<subworkflow>.config`
  (input_checks, prepare_reference, prepare_regions, coverage, preprocess_vcf,
  benchmark, discordance, annotation, known_variants, report). `conf/modules.config` keeps
  only the default `publishDir` and the `includeConfig` lines. Pure move.
- **A4** ~~Shared helper module~~ dropped (user, 2026-10-13: no io helper script). `read_bed`,
  `read_fai` and the pysam readers become named functions in `check_reference.py`
  (B1); `known_variants.py` imports them from there.

### Phase B — all input checks first (`feat-input-checks`) — highest value

- **B1** Rewrite `check_reference.py` on pysam:
  - `vcf`: `pysam.VariantFile` header, contigs and samples; REF alleles compared with
    `pysam.FastaFile` after contig mapping; mismatch rate vs. `--max-ref-mismatch`.
  - `bam`: `pysam.AlignmentFile` header (contigs, lengths).
  - `bed`: `read_bed`.
  - No length fixing (mismatch → error); chrM dropped and reported (`--keep-mito`).
  - A short docstring for every function (e.g. `ref_alleles`).
  - INPUT/OUTPUT column description in the module docstring.
  - Tests updated.
- **B2** Module `CHECK_REFERENCE`: input `tuple val(meta), path(input), path(input_index)`
  (renamed from `facts`), plus `fasta`/`fai`.
- **B3** New subworkflow `CHECK_INPUTS` (`subworkflows/local/check_inputs/`), run right after
  `CHECK_SAMPLESHEET`. It checks every VCF, BAM, target / regions / exclude / stratification /
  known-list BED, and the known lists (`known_variants.py prepare` with REF check against
  the FASTA). All reports are collected; every other subworkflow waits for them. It emits the
  rename maps and contig lists used downstream.
  - Params (default true): `check_vcf_contigs`, `check_ref_alleles`, `check_bam_contigs`,
    `check_bed_contigs`, `check_known_variants`.
- **B4** Remove the now-duplicate steps:
  - `BCFTOOLS_VIEW_HEADER`, `BCFTOOLS_QUERY_CONTIGS`, `CHECK_REF_ALLELES(_KNOWN)` from
    PREPROCESS_VCF / PREPARE_KNOWN_VARIANTS;
  - `SAMTOOLS_IDXSTATS` from COVERAGE;
  - the per-subworkflow `CHECK_REFERENCE_*` calls.

  These subworkflows now take the rename/keep outputs of `CHECK_INPUTS`.
  `bcftools norm --check-ref x` stays (it removes mismatches; its log goes to MultiQC).
- **B5** Minimum region size: `--min_region_bp` (default 10000). `region_summary.py` stops
  with an error if the common region or a sample's callable region is smaller (test
  configs lower it).

### Phase C — options and parameters (`feat-options`)

- **C1** `--callable_regions` / `--read_evidence` (default true):
  - `check_sheets.py` requires BAMs when they are on;
  - COVERAGE resp. the mpileup branch is skipped when off;
  - `classify_discordance.py` runs without evidence (no evidence columns/classes).
- **C2** `--min_coverage_tumor` / `--min_coverage_normal` replace `--min_coverage`:
  - mosdepth quantize threshold by `meta.role`; classification uses the tumor threshold.
  - `--require_normal_coverage`: today it only adds normals to the callable intersection.
    Now it also stops in `check_sheets.py` if a paired row has no normal BAM. Docs answer
    the user's question.
- **C3** Early PASS filter in PREPROCESS_VCF:
  - after step 2: `bcftools view -f PASS,.` (compare branch) and the complement
    (non-PASS branch);
  - the non-PASS branch goes through the same steps 3–7 (shared subworkflow
    `NORMALIZE_VCF` called twice) and is passed to DISCORDANCE;
  - new class `filtered_non_pass` (+ FILTER value);
  - update CLAUDE.md chain + D23.
- **C4** `harmonize_vcf.py`:
  - `--vaf-source auto|af|ad` (auto prefers AF);
  - `--max-missing-vaf` as pipeline params `vaf_source`, `max_missing_vaf`;
  - no VAF → warning, header line `##variantcomp_vaf=unavailable`, read by
    `benchmark_metrics.py` / `classify_discordance.py` → VAF strata and power NA;
  - docstring lists the supported VCF versions (4.1–4.3) and tested callers/versions;
  - `main()` split into `read_vcf`, `check_samples`, `harmonize_record`, `write_summary`.
- **C5** Parameter audit: every hard-coded threshold becomes a documented param in
  `nextflow.config` + `nextflow_schema.json`. Candidates:
  - mpileup MQ/BQ/max depth (`evidence_min_mapq`, `evidence_min_baseq`, `evidence_max_depth`);
  - amplicon duplicate flag;
  - known-list REF mismatch;
  - vcfeval flags;
  - report low-VAF;
  - `germline_*`.

  Done via `ext.args` in `conf/modules/*.config`.

### Phase D — annotation subworkflow (`feat-vep-annotation`)

- **D1** Install the nf-core modules `ensemblvep/vep` and `msisensorpro/scan`. Add a local
  module `BCFTOOLS_SPLITVEP` (one call to `bcftools +split-vep`).
- **D2** New subworkflow `ANNOTATE_VARIANTS`, replacing part 3 of DISCORDANCE:
  1. sites VCF of all compared variants (`annotate_variants.py to-vcf`, tabix);
  2. VEP;
  3. `+split-vep` (Consequence, SYMBOL, IMPACT, VARIANT_CLASS; PICK);
  4. repeats: `msisensorpro scan` on the FASTA → BED (`--repeats_bed` reuses a
     precomputed one);
  5. `bedtools intersect -c` for repeats and stratification BEDs;
  6. `annotate_variants.py merge`.

  Remove `--gene_bed`, `BEDTOOLS_MAP_GENES`.
- **D3** Params:
  - `vep` (true), `vep_cache`, `vep_cache_version` (115), `vep_species`, `vep_genome`;
  - `vep_gff` (offline without cache);
  - VEP flags (`vep_args`);
  - `repeats_bed`, `msisensor_args`.

  Schema docs.
- **D4** Report: discordance stratified by
  - consequence group (synonymous / non-synonymous / non-coding, mapped from SO terms in
    `variantcomp_report.py`);
  - repeat yes/no and repeat unit;
  - VEP variant class;
  - impact.

  Gene table from SYMBOL.
- **D5** Test data:
  - GENCODE chr22 GFF3 slice for SEQC2 (`tests/seqc2/` prepare script);
  - a tiny synthetic GFF;
  - test configs use `vep_gff`.

### Phase E — summaries, known variants, statistics (`feat-reporting`)

- **E1** All step summaries go into MultiQC. New `bin/mqc_tables.py` + module `MQC_TABLES`:
  - harmonize summaries;
  - `bcftools norm` logs;
  - reference checks;
  - region / callable summaries;
  - sheet checks;
  - required-variant summary.

  The native MultiQC VEP module reads the VEP stats. All are collected into `ch_multiqc_files`.
- **E2** `known_variants.py`:
  - remove `MAX_REQUIRED` (and its warning);
  - split `prepare` / `check` into sub-functions (`read_lists`, `map_variants`,
    `write_sites_vcf`, `load_technology`, `status_of`, ...);
  - highlight the comparison logic with comments.
- **E3** Recall of required variants per sample × technology and over all samples
  (`known_summary.py` + MultiQC + report).
- **E4** `variantcomp_report.py`: use `scipy.stats.false_discovery_control` instead of
  `bh()`; keep `fisher_exact`. Update the dependency pins (scipy ≥ 1.11) in the env and the
  report container (rebuild the Seqera container if needed).

### Phase F — code clarity (`feat-code-clarity`)

- **F1** BED_INTERSECT_ALL contains `bedtools sort`. The docs/code explain why not
  `bedtools intersect`: `-b A B C` gives overlap with *any* file (union); an N-way
  intersection needs N−1 chained calls, which a fixed DSL2 module chain can't express.
  `multiinter` does it in one call.
- **F2** Refactor `mosdepth_callable.py`:
  - functions `open_mosdepth_bed`, `read_rename_map`, `is_callable`, `write_callable`;
  - clear variable names;
  - comments.
- **F3** BENCHMARK subworkflow: inline comments per process.
- **F4** Every `bin/` script gets an INPUT/OUTPUT docstring with all columns and types
  (contig name, 0-based start, ...) and a short docstring for every function.

### Phase G — documentation (`feat-docs`)

- **G1** `docs/output.md`: every published file, its columns, column types and meaning.
- **G2** Metro map (`scripts/make_metro_map.py`), as the user described:
  - **Central green line:** VCF → PREPROCESS_VCF → BENCHMARK → DISCORDANCE →
    ANNOTATE_VARIANTS → REPORT.
  - **Input lines:** FASTA, BED, BAM, known lists and the VEP cache each start at their
    input, run through their own subworkflow (PREPARE_REFERENCE, PREPARE_REGIONS → COVERAGE,
    ...) and **merge** into the steps that consume them, so each line is one path.
  - Boxes per subworkflow stay. The README table is regenerated.
- **G3** Record the design essence in CLAUDE.md (see below). Add decisions D23–D27 to
  `docs/pipeline-plan.md` §11. Update `docs/usage.md` (new params, BAM requirements), the
  CHANGELOG and CITATIONS (Ensembl VEP: McLaren et al. 2016; MSIsensor-pro: Jia et al. 2020).
- **G4** New file **`docs/roadmap.md`** (content below), linked from the README and the plan §9.

### Phase H — long tests and release

- **H1** Update the nf-tests:
  - CHECK_INPUTS, PREPROCESS_VCF (PASS side branch), COVERAGE (tumor/normal thresholds),
    ANNOTATE_VARIANTS;
  - test configs: synthetic without BAMs → `callable_regions=false`, `read_evidence=false`.
- **H2** Run all nf-tests and the `test` / `test_seqc2` pipelines on the Pi. Check the
  designed synthetic counts and the SEQC2 metrics against v0.14.0-dev; update the snapshots
  only for explained differences.
- **H3** `nf-core pipelines lint`. Version bump `0.15.0-dev` as its own commit, tag. Give a
  milestone summary to the user (done / tested where / open; T2 command for the HPC).

## Essence for CLAUDE.md (G3)

- All input checks run first (CHECK_INPUTS); each can be switched off.
- No silent fixes of contig lengths; chrM is dropped and reported.
- Python: pysam for VCF/BAM/FASTA, readers as named functions (no shared helper module), black + ruff, an INPUT/OUTPUT
  docstring per script and a docstring per function.
- Every threshold of a script or tool is a documented pipeline param.
- Every step summary reaches MultiQC.
- One config file per subworkflow in `conf/modules/`.
- Optional features (callable regions, read evidence, VEP) are switched by params; when
  they are on, their inputs are required.
- Annotation is done by VEP (+ msisensor-pro repeats), the same way as in
  variantinterpretation.
- Statistics come from scipy/statsmodels, not own code.

## Files touched (representative)

- `bin/check_reference.py`, `bin/harmonize_vcf.py`, `bin/known_variants.py`,
  `bin/mosdepth_callable.py`, `bin/variantcomp_report.py`, `bin/annotate_variants.py`,
  new `bin/mqc_tables.py`
- `subworkflows/local/{check_inputs,annotate_variants}/` (new), `preprocess_vcf`,
  `coverage`, `discordance`, `prepare_regions`, `bed_intersect_all`, `benchmark`,
  `prepare_known_variants`
- `workflows/variantcomp.nf`, `conf/modules/*.config`, `nextflow.config`,
  `nextflow_schema.json`, `envs/variantcomp-dev.yaml`, `pyproject.toml`
- `docs/{usage,output,pipeline-plan,roadmap}.md`, `README.md`, `CLAUDE.md`,
  `CITATIONS.md`, `scripts/make_metro_map.py`
- Test-data repo: GFF slices, sheet/config updates

## Verification

- Per commit: `black --check bin scripts tests`, `ruff check`, `pytest tests/python`,
  `nextflow lint`.
- At the end (H2):
  - nf-tests for all local modules/subworkflows;
  - `nf-test test tests/default.nf.test tests/seqc2_small.nf.test`;
  - the synthetic counts (S1 wes 9/1/9/1, panel 8/2/8/1; S2 wes 2/1/2/1, panel 3/0/3/1)
    must stay the same;
  - SEQC2 metrics are compared with v0.14.0-dev.
- Failure tests: a wrong-build BED, a chrM-length VCF and a missing BAM with
  `callable_regions=true` must fail **in CHECK_INPUTS within the first minutes**;
  switching the toggle off must let the run pass.
- `nf-core pipelines lint`: 0 failed.

## Content of `docs/roadmap.md` (G4)

1. **Validation:**
   - real DRAGEN v4.4.4 VCFs (8 × 3 cohort);
   - full SEQC2 WES+WGS on the HPC (T2);
   - GIAB HG008 tumor/normal truth set.
2. **Interactive reports:** datavzrd tables of discordant variants (design input from
   variantinterpretation, not done yet).
3. **all_pairs comparison:** concordance matrix (Jaccard) and UpSet plot of all
   technologies.
4. **Statistics:**
   - Wilson/bootstrap confidence intervals for recall/precision;
   - per-gene discordance normalized by the callable bp of the gene;
   - purity/ploidy-aware detection power.
5. **Read-level FP features:** strand bias, MQ/BQ, read position, soft clips at
   discordant sites; FP clustering.
6. **Caller support:**
   - harmonization profiles for Strelka2 (tier counts), VarScan2, Mutect2 4.x;
   - DRAGEN-specific FILTER/INFO mapping.
7. **VEP extensions:**
   - plugins (REVEL, SpliceAI, AlphaMissense);
   - ClinVar/COSMIC lookups of discordant variants;
   - HGVS input for required known variants.
8. **Discordance class `outside_callable`:** define it per technology or drop it (open
   since v0.13).
9. **Scale:**
   - sharding by chromosome for large WGS cohorts;
   - cached reference artefacts (msisensor repeats, VEP cache index);
   - resource labels tuned on the HPC.
10. **Infrastructure:**
    - CI (GitHub/GitLab) running lint, pytest and stub runs;
    - conda lock files;
    - pinned Singularity images in the offline bundle;
    - pre-commit hooks (black, ruff, prettier).
11. **Inputs:**
    - UMI/amplicon-specific settings;
    - FASTQ entry point (out of scope for v1);
    - tumor-only mode.
12. **Output schema:** nf-schema output definitions and an RO-Crate for provenance.
