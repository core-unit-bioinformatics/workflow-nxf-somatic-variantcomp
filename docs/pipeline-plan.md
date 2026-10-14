# Pipeline plan: cross-technology somatic VCF comparison

- Status: **approved 2026-10-07**, development started. The decisions in §11
  and the addition in §12 take precedence over earlier sections.
- Repository: `core-unit-bioinformatics/workflow-nxf-somatic-variantcomp`
  (internal nf-core template name `variantcomp`).

## 1. Goal

We compare somatic small-variant calls (SNVs and indels) of the **same
tumor/normal pair**, sequenced with **different technologies** (WGS, WES,
targeted panels; Illumina short reads). One technology, set in the config, is
the ground truth. The main outputs are:

- recall, precision and F1 per sample and technology pair, plus aggregated
  values;
- tables of the variants that differ between technologies, with statistics
  that explain *why* they differ: variant type, indel length, VAF, gene,
  repeat/homopolymer context, coverage, filter status, and read evidence.

### How this differs from nf-core/variantbenchmarking

variantbenchmarking compares *different callers on the same data*. Here the
*input data differ* (assay, depth, limit of detection, target regions, and
often the caller as well). Getting meaningful numbers therefore needs:

1. strict, identical normalization of every VCF;
2. restriction to regions that **all** technologies cover, both by design
   (BED) and in practice (coverage);
3. metrics split by VAF and adjusted for detection power. A ~1000x panel and a
   ~40x WGS have very different limits of detection;
4. an explanation for every discordant variant: missed vs. filtered vs. not
   covered vs. no read evidence.

## 2. User requirements (from the request on 2026-10-07)

| # | Requirement | Where it is covered |
|---|---|---|
| R1 | ≥2 technologies, ≥1 sample (use case: 8 samples × 3 technologies) | Samplesheet design, §4 |
| R2 | Input: unfiltered VCFs; FASTQ/BAM optional | §4. BAM is optional and enables coverage and evidence checks. FASTQ is out of scope for v1 (use nf-core/sarek upstream) |
| R3 | Ground truth defined in the config | `params.truth_technology`, §4 |
| R4 | Recall / precision statistics | BENCHMARK, §5.5 |
| R5 | Tables of discordant variants + statistics (indels, genes, repeats) | DISCORDANCE + STATS_REPORT, §5.6–5.7 |
| R6 | Target BEDs as input; comparison restricted to the intersection of all BEDs | REGIONS, §5.1 |
| R7 | Design inspired by cio-abcd/variantinterpretation | §3 |
| R8 | git from the start, every feature/step versioned | §9 |
| R9 | Conda env for dev/test tooling | §8.1 |
| R10 | Test with public data: same sample sequenced with WGS and targeted assays, somatic tumor/normal | §8.2 |
| R11 | If it can't run on the Pi, prepare a bundle to transfer to another system | §8.4 |
| R12 | Prefer existing tools (bcftools, …) and nf-core modules | §6 |
| R13 | Python over R for scripts | §6 |

## 3. Design inspiration: cio-abcd/variantinterpretation

Patterns we take over:

- nf-core DSL2 template, nf-schema samplesheet validation, one container per
  process, nf-test.
- VCF checks before anything else: index, reference/contig consistency, sample
  names.
- Split multiallelic sites (needed for downstream tools that only handle
  biallelic records) and left-align against the FASTA.
- Filtering marks variants instead of silently dropping them, wherever
  possible.
- Reports: TSV export plus interactive **datavzrd** HTML tables configured from
  YAML templates in `assets/`, and MultiQC for the summary.
- The user's previous project ran variantinterpretation v1.1.0 offline on a
  PBS Pro HPC with Singularity. The HPC profile and the offline design follow
  that setup.

## 4. Input

### 4.1 Samplesheet (`.tsv`, CUBI rules: `sample` column, `n/a` for empty fields)

| Column | Required | Description |
|---|---|---|
| `sample` | yes | Patient / case ID. Identical across technologies |
| `technology` | yes | Free label, e.g. `wgs`, `wes`, `panel_oca` |
| `library_type` | yes | `wgs` \| `capture` \| `amplicon`. Controls duplicate handling, see §5.2 |
| `tumor_id` / `normal_id` | yes / opt. | Sample names inside the VCF. `normal_id` is `n/a` for tumor-only calls |
| `caller` | yes | Key into `assets/caller_profiles.yaml` (where VAF/DP/AD live, germline tags, …) |
| `vcf` (+ `vcf_index`) | yes | **Unfiltered** somatic VCF |
| `tumor_bam`/`tumor_bai`, `normal_bam`/`normal_bai` | no | Enables coverage-based callable regions and read-evidence checks |
| `target_bed` | no | `n/a` = genome-wide (WGS). Must be identical for all rows of a technology |

Validation at startup:

- each sample has the truth technology plus at least one other;
- each technology has one BED;
- each (sample, technology) pair is unique;
- the IDs exist in the VCF headers;
- all VCFs and BEDs use the configured genome build. There is no liftover in
  v1.

An external truth set (e.g. the SEQC2 reference calls) can be added as a
pseudo-technology: a VCF with its confident-region BED and no BAM.

### 4.2 Main parameters

| Parameter | Default | Purpose |
|---|---|---|
| `input` | – | Samplesheet |
| `fasta`, `fasta_fai` | – | Reference. One build for all inputs |
| `truth_technology` | – | Ground truth technology |
| `comparison_mode` | `truth` | `truth` (truth vs. each) or `all_pairs` (additionally all-vs-all concordance + UpSet) |
| `filter_pass_only` | `true` | Uniform FILTER handling |
| `min_vaf`, `min_depth`, `min_alt_reads` | `0.0`, `0`, `0` | Uniform hard filters, applied identically to all technologies |
| `min_coverage_tumor`, `min_coverage_normal` | `10` | Callable thresholds (from the old R script); replace `min_coverage` (v0.15.0) |
| `callable_regions`, `read_evidence` | `true` | Use the BAMs for callable regions / read evidence; need BAMs (v0.15.0) |
| `callable_mode` | `base` | `base` (per-base, mosdepth `--quantize`) or `region` (mean per target region, as in the old script) |
| `require_normal_coverage` | `true` | The normal must be callable too |
| `regions_bed`, `exclude_bed` | – | Optional extra restriction / blacklist |
| `vaf_bins` | `0,0.05,0.1,0.2,0.5,1` | VAF strata |
| `detection_power` | `0.95` | Threshold for power-adjusted recall |
| `stratification_beds` | – | TSV `name<TAB>path`, e.g. GIAB stratifications (homopolymers, STRs, segdups, low mappability, GC extremes) |
| `vep` | `true` | Ensembl VEP annotation (D25); `gene_bed` was removed in v0.15.0 |
| `vep_cache` or `vep_gtf`, `vep_cache_version`, `vep_genome`, `vep_species` | – | Annotation resources (offline) |
| `annotate_repeats`, `repeats_bed` | `true`, – | msisensor-pro repeat regions (D25) |
| `check_vcf_contigs`, `check_ref_alleles`, `check_bam_contigs`, `check_bed_contigs`, `check_known_variants`, `keep_mito` | `true` (`keep_mito` `false`) | Upfront input checks (D24, D26) |

## 5. Workflow

```
INPUT_CHECK (nf-schema + cross-row checks)
│
├─ 5.1 REGIONS ───────── per technology BED: sort → merge → rename contigs → keep contigs in .fai
│                         intersect ALL target BEDs  (∩ regions_bed, − exclude_bed)
│                         → common_targets.bed  (WGS-only comparison: genome or regions_bed)
│
├─ 5.2 COVERAGE (opt.) ─ mosdepth per BAM on common_targets.bed
│                         base mode:   --quantize 0:<min_cov>: → CALLABLE intervals
│                         region mode: --by common_targets.bed, mean ≥ min_cov
│                         intersect CALLABLE of all technologies (tumor [+ normal]) per sample
│                         → callable_<sample>.bed + coverage QC (histogram/scatter, Python)
│                         amplicon libraries: don't exclude duplicate-flagged reads
│
├─ 5.3 VCF_PREPROCESS ── per VCF (bcftools wherever possible):
│     1. reheader contigs from .fai           (bcftools reheader --fai; fixes chrM length bug)
│     2. rename contigs (optional)            (bcftools annotate --rename-chrs)
│     3. select tumor [+ normal] columns      (bcftools view -s)
│     4. drop gVCF/reference records          (exclude ALT="<*>"/"<NON_REF>", trim alt alleles)
│     5. REF check + split + left-align       (bcftools norm -f ref --check-ref x -m-any)
│     6. atomize MNPs/complex                 (bcftools norm --atomize)
│     7. deduplicate                          (bcftools norm -d exact)
│     8. harmonize tags                       (bcftools +fill-tags FORMAT/VAF, bcftools annotate
│                                              renames from caller profile; GT → 0/1 via +setGT;
│                                              small Python helper only where bcftools can't)
│     9. uniform filtering (mark, then view)  → *.norm.vcf.gz (all) + *.filtered.vcf.gz
│    10. restrict to callable_<sample>.bed or common_targets.bed (bcftools view -T)
│
├─ 5.4 REFERENCE ─────── rtgtools format (SDF from FASTA), once
│
├─ 5.5 BENCHMARK ─────── per sample × (truth, query):
│     rtgtools vcfeval --squash-ploidy --decompose     → TP / FP / FN VCFs (haplotype-aware)
│     bcftools isec (exact allele match)                → representation cross-check
│     position-only match (bcftools isec -c all)
│     metrics: recall, precision, F1 — overall, SNV vs indel, VAF bins, strata
│     power-adjusted recall: truth variants where the query had ≥ detection_power to see them,
│       from query depth (bcftools mpileup / mosdepth) + binomial model (Python/scipy)
│     comparison_mode=all_pairs: concordance matrix (Jaccard) + multi-way bcftools isec → UpSet
│     aggregation over samples: micro (summed counts) and macro (mean ± sd)
│
├─ 5.6 DISCORDANCE ───── for every FP/FN record:
│     type (SNV/MNV/ins/del), indel length, VAF/DP in both technologies
│     status in the other technology:
│       absent | present_but_filtered (from *.norm.vcf.gz) | outside_callable
│     read evidence (if BAMs): bcftools mpileup -a AD,DP on the other technology's tumor & normal
│       → no_evidence | subthreshold_evidence | germline_in_normal | low_coverage
│     annotation: bcftools annotate with gene BED (or VEP), stratification BEDs
│       (homopolymer/STR/segdup/low-mappability), technology target-edge distance
│     export: bcftools query → TSV
│
└─ 5.7 STATS_REPORT ──── Python (pandas/scipy/matplotlib):
      enrichment discordant vs concordant (Fisher's exact, BH-corrected):
        indel fraction, per stratum, low-VAF, per gene (normalized by callable bp of the gene)
      plots: VAF-stratified recall, discordance class bars, VAF scatter truth vs query,
        coverage QC, UpSet (upsetplot)
      MultiQC (custom-content summary tables) + datavzrd (interactive discordant tables)
```

### 5.8 Output layout

```
results/
├── regions/           common_targets.bed, callable_<sample>.bed, region_summary.tsv
├── coverage/          mosdepth outputs, coverage QC plots
├── vcf_preprocessed/  <sample>_<tech>.{norm,filtered,callable}.vcf.gz(.tbi)
├── benchmark/<sample>/<truth>_vs_<query>/   vcfeval TP/FP/FN VCFs, isec outputs
├── stats/             metrics_per_pair.tsv, metrics_aggregated.tsv, metrics_by_vaf.tsv,
│                      metrics_by_stratum.tsv, concordance_matrix.tsv, upset_<sample>.png
├── discordance/       discordant_variants.tsv, enrichment_tests.tsv, gene_summary.tsv
├── report/            datavzrd/, multiqc_report.html
└── pipeline_info/
```

## 6. Tooling principles

1. **Use existing tools first.** Choose in this order:
   1. an nf-core module;
   2. bcftools / bedtools / samtools / mosdepth / rtg-tools commands inside a
      local module;
   3. a small Python script in `bin/`.

   No R. Python scripts get a `#!/usr/bin/env python3` shebang, are
   executable, use argparse, and contain no workflow-engine logic.
2. **nf-core modules to install:**
   - bcftools: `index`, `reheader`, `annotate`, `view`, `norm`, `isec`,
     `query`, `mpileup`, `pluginsetgt`/`plugin*` where available, `stats`
   - bedtools: `sort`, `merge`, `intersect`, `subtract`
   - `mosdepth`
   - rtgtools: `format`, `vcfeval`
   - `samtools/faidx`, `tabix/bgziptabix`, `ensemblvep/vep`, `datavzrd`,
     `multiqc`
3. **Local modules / `bin/` scripts (Python):**
   - `compute_metrics.py`
   - `power_adjusted_recall.py`
   - `classify_discordance.py`
   - `enrichment_stats.py`
   - `plot_report.py`
   - `check_samplesheet` logic beyond nf-schema, if needed
4. Container and conda specs come from the nf-core modules. Local modules get
   biocontainers or Seqera/Wave containers with pinned versions.

## 7. Scientific caveats (to put in the docs and the report)

- **Different callers per technology mix two effects.** Discordance then
  comes from both the sequencing and the caller. The report shows the caller
  for each pair.
- **Tumor-only calling produces germline false positives.** The
  `germline_in_normal` evidence class and an optional population-AF filter
  (VEP/gnomAD) make this visible.
- **Separate DNA aliquots and subclonality.** Low-VAF discordance is not
  automatically a technical error.
- **Amplicon libraries.** Duplicate marking destroys their coverage, so
  coverage and pileup steps must not filter duplicates for
  `library_type=amplicon`. Primer trimming is the user's upstream job, and the
  docs say so.
- **Validation panels are biased.** Panels designed around known calls (e.g.
  SEQC2 AmpliSeq) overestimate recall compared with clinical panels.

## 8. Development environment, testing, and HPC transfer

### 8.1 Conda dev environment (`envs/variantcomp-dev.yaml`)

Pinned versions: nextflow (26.04.x, the current local version), nf-core tools
4.1.x, nf-test, python 3.12, pysam, pandas, scipy, matplotlib, upsetplot,
pytest, ruff, bcftools/samtools/htslib 1.24, bedtools 2.31.1, mosdepth 0.3.x,
rtg-tools 3.13.

Checked on 2026-10-07 via the anaconda.org API:

- bcftools, samtools, htslib, bedtools, mosdepth and pysam have
  `linux-aarch64` builds;
- rtg-tools, vembrane, ensembl-vep and multiqc are `noarch`;
- datavzrd has not been checked yet.

So the pipeline should run on the Raspberry Pi (aarch64, 4 cores, 3 GB RAM)
with `-profile conda` on **tiny** data. Docker isn't installed, and
Singularity images are amd64, so container runs go to the HPC.

### 8.2 Test data: literature research

| Dataset | Content | Access | Fit |
|---|---|---|---|
| **SEQC2 HCC1395 / HCC1395BL** (breast cancer cell line + matched B-lymphoblastoid normal). Fang et al. 2021, *Nat Biotech* 39:1151, doi:10.1038/s41587-021-00993-6. Xiao et al. 2021, *Nat Biotech* 39:1141, doi:10.1038/s41587-021-00994-5 | WGS and WES (Agilent SureSelect v6+UTR, BED provided) from 6 centers with replicates; WES also on Ion Torrent; **AmpliSeq validation panel** (Illumina MiSeq 2×150, ~2000×, 1,368 amplicons, T + N in triplicate); per-run VCFs (MuTect2, Strelka, SomaticSniper × BWA/Bowtie2/NovoAlign); BAMs; truth set v1.2.1 + High-Confidence_Regions_v1.2.bed; GRCh38 (`GRCh38.d1.vd1`) | Open: SRA SRP162370, `ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/` (checked 2026-10-07) | **Primary choice**: same DNA batch, WGS + WES + panel, tumor/normal, plus an independent truth set |
| GIAB HG008 (pancreatic adenocarcinoma cell line + matched normal tissue). McDaniel et al., bioRxiv 2024 / 2025 | WGS on many platforms (Illumina, Element, Ultima, PacBio, ONT, …). Somatic SV/CNV draft benchmarks; small-variant benchmark in progress | Open (GIAB FTP) | Secondary: cross-platform WGS. No WES/panel data found |
| COLO829 / COLO829BL. Craig et al. 2016, *Sci Rep* 6:24607 | Multi-center WGS, somatic truth set | dbGaP phs001839 / EGA (controlled) | Not suitable (controlled access) |

**Notes on SEQC2:**

- The AmpliSeq panel BED is not on the FTP. We derive it from AmpliSeq normal
  coverage (mosdepth, ≥ min depth), or take it from the SEQC2 website if
  available.
- There are no AmpliSeq VCFs, so calling has to happen upstream (§8.3).
- The AmpliSeq panel was designed around candidate calls, so it is biased
  (§7).
- The pipeline expects patients, but SEQC2 has a single cell-line pair.
  Center/replicate pairs are used as **pseudo-samples**: WGS_X_n and WES_X_n
  from the same center X, and AmpliSeq replicate n. This exercises the
  multi-sample aggregation.

### 8.3 Test tiers

| Tier | Data | Runs on | Purpose |
|---|---|---|---|
| T0 unit | Hand-made mini VCFs on a ~10 kb synthetic reference. Edge cases: multiallelic, MNP vs. split SNVs, unnormalized indels, gVCF blocks, chrM length mismatch, `chr`/no-`chr`, missing AF field, tumor-only | Pi | pytest for `bin/` scripts, nf-test for modules/subworkflows |
| T1 `-profile test` | SEQC2 slices: a few Mb on chr22 (GRCh38) or windows around AmpliSeq amplicons. VCF subsets (WGS/WES MuTect2 + Strelka), BAM slices fetched remotely by region (`samtools view` on the FTP URLs with index), AmpliSeq slices, truth subset, matching FASTA. Target size ≤ ~10 MB | Pi (conda) + CI | End-to-end pipeline test, nf-test snapshot |
| T2 `-profile test_full` | Full SEQC2, one center (e.g. IL or NV): WGS (truth: SEQC2 v1.2.1, and WGS as technology truth), WES, AmpliSeq. Upstream calling with **nf-core/sarek** (`--step variant_calling`, Mutect2 [+ Strelka2]) from the provided BAMs, restricted to the union of targets, so all three technologies use the same caller | HPC | Realistic validation. Sanity check: WGS vs. SEQC2 truth close to the published numbers |

Test-data preparation is codified, with no manual steps:

- `scripts/prepare_seqc2_testdata.py` (Python, calls samtools/bcftools/bedtools)
- `scripts/envs/testdata.yaml`
- `docs/testing.md`

T1 files are either committed to `tests/data/` (if ≤ ~10 MB) or published as
a separate test-data repo. That is decision D6.

### 8.4 HPC transfer bundle

Everything needed to run offline elsewhere:

1. A pipeline release tag.
2. `nf-core pipelines download` → pipeline code + Singularity images (and
   sarek for T2), as an archive.
3. `conf/` profile template for PBS Pro / offline Singularity (modelled on the
   old `run.config`), without hardcoded paths.
4. A params-file template (`.yaml`).
5. Test-data download/preparation script plus md5 checks.
6. `docs/hpc-run.md`: step-by-step commands, expected outputs, and how to send
   results back.

## 9. Git and versioning

- Local git from the first commit, with branches `main` and `dev`. Every
  feature goes on `feat-<name>` off `dev` and is rebased-and-merged into `dev`
  to keep the history linear. Merge or squash commits are not allowed.
- Each feature ships with nf-test tests, docs, and a CHANGELOG entry.
- Version bumps are **separate commits**: `manifest.version` plus the
  CHANGELOG header only. They are tagged `vX.Y.Z-dev` on `dev`. `1.0.0` comes
  from a fast-forward of `dev` into `main`.
- Remotes (GitHub / GitLab@HHU) are set up only on request (decision D5).

| Version | Branch | Content |
|---|---|---|
| `0.1.0-dev` | `feat-template` | nf-core template, dev conda env, CUBI metadata (LICENSE, CITATION.cff, pyproject.toml), README, CLAUDE.md + this plan |
| `0.2.0-dev` | `feat-input` | Samplesheet schema + cross-row validation, channel structure, T0 data |
| `0.3.0-dev` | `feat-testdata` | SEQC2 prep script, T1 test profile |
| `0.4.0-dev` | `feat-regions` | BED harmonization + N-way intersection |
| `0.5.0-dev` | `feat-vcf-preprocess` | Header fix, sample selection, normalization, tag harmonization, uniform filters |
| `0.6.0-dev` | `feat-coverage` | mosdepth callable regions (base/region), coverage QC |
| `0.7.0-dev` | `feat-benchmark` | vcfeval + isec, metrics (overall/VAF/strata), power-adjusted recall, aggregation, all_pairs mode |
| `0.8.0-dev` | `feat-discordance` | Discordance classification, annotation, BAM evidence |
| `0.9.0-dev` | `feat-report` | Enrichment stats, plots, MultiQC + datavzrd |
| `0.10.0-dev` | `feat-known-variants` | Required clinical variants + external truth-set benchmark (§12) |
| `0.11.0-dev` | `feat-hpc-bundle` | HPC profile template, download bundle, T2 instructions |
| `1.0.0` | `dev` → `main` | After a successful T2 run on the HPC and review by the user |

Versions `0.12.0-dev` to `0.14.0-dev` refactored the modules and sheet reading (see
CHANGELOG); `0.15.0-dev` implements the feedback of 2026-10-13
([`plan-feedback-20261013.md`](plan-feedback-20261013.md)). Further work is collected in
[`roadmap.md`](roadmap.md).

Status (2026-10-08): `0.1.0-dev` to `0.11.0-dev` implemented and tested on the Raspberry Pi
(pytest, nf-test, synthetic and SEQC2 small data with `-profile conda`). Not yet done:
the T2 run on the HPC with Singularity, and datavzrd reports (MultiQC custom content is used
for now).

## 10. Workspace layout (agreed)

```
vcf-comp-dev/                                  # workspace (no git)
├── CLAUDE.md                                  # imports the repo CLAUDE.md
├── workflow-nxf-somatic-variantcomp/          # pipeline repo (git): code, tests/, docs/, envs/
└── workflow-nxf-somatic-variantcomp-testdata/ # test-data repo (git): data only
```

## 11. Decisions (agreed with the user on 2026-10-07)

| ID | Question | Decision |
|---|---|---|
| D1 | Name / prefix | Repo and pipeline `workflow-nxf-somatic-variantcomp`, prefix `core-unit-bioinformatics` (manifest `core-unit-bioinformatics/workflow-nxf-somatic-variantcomp`). Internal nf-core template name `variantcomp` (Nextflow identifiers can't contain dashes) |
| D2 | Genome build (length fixes reversed by D24) | GRCh38 for tests (SEQC2); GRCh37 via params. **New requirement:** check that all inputs (VCFs, BEDs, BAMs) are on the configured reference (contig names + lengths, REF alleles, build fingerprint via chr1 length) and auto-fix what can be fixed safely: contig header lengths (e.g. chrM 16570 vs 16571), `chr`/no-`chr` and `MT`/`chrM` naming, missing contig lines. Hard error on a build mismatch or a REF mismatch rate above a threshold (no liftover) |
| D3 | Callers | User data: **DRAGEN v4.4.4**, all samples **tumor/normal paired** (docs: help.dragen.illumina.com). Test data (SEQC2): MuTect2 VCFs. `caller_profiles.yaml` starts with `dragen` and `mutect2` |
| D4 | Truth | One technology for all samples (`truth_technology`) |
| D5 | Remotes | Commit locally only for now |
| D6 | Test data | Data in a separate repo/folder (`workflow-nxf-somatic-variantcomp-testdata/`). The test code (download/prepare/run script) lives in the workflow repo |
| D7 | Annotation (reversed by D25) | `gene_bed` default, VEP optional |
| D8 | Layout | Legacy scripts deleted. Workspace holds two subfolders: the pipeline repo and the test-data repo |
| D9 | SEQC2 scope | WES + WGS only for now (no AmpliSeq) |
| D10 | Call definition (implementation, 2026-10-07) | Every record of a somatic VCF is a call of the tumor sample: after splitting, GT is set to tumor `0/1` and normal `0/0` (callers fill GT differently, and split multiallelic records can otherwise be lost). VAF/depth/alt counts are computed from FORMAT/AD for every caller (`VC_*` INFO tags), so no caller profiles are needed. |
| D12 | Exact/position matching (implementation, 2026-10-08; reversed by D17) | The `match` strata are computed as set comparisons of the compared VCFs in `benchmark_metrics.py` instead of with `bcftools isec`: it is one pass over two small VCFs and needs no extra task per pair. |
| D14 | Sheet format (user, 2026-10-08) | Samplesheet, stratification and known-variant sheets are CSV (overrides the CUBI TSV default for this project). No index columns: indices next to the files are used, missing ones are created with `bcftools index` / `samtools index`. |
| D15 | Module architecture (user, 2026-10-08) | One tool per module, preferably nf-core modules (patched if needed); no tool pipes inside modules; logic only in `bin/` scripts called by modules, with arguments set in `conf/modules.config`. See CLAUDE.md, "General nf-core workflow development guidelines". |
| D16 | Checks (user, 2026-10-08) | Validation and pass/fail decisions in Python scripts with meaningful errors, not in Nextflow code; be strict, error rather than wrong results. |
| D17 | Custom scripts (user, 2026-10-08) | Replace custom code where tools can do the job (e.g. contig/length extraction with `bcftools view -h` / `samtools view -H`, REF check with `bcftools +fixref`, exact matching with `bcftools isec`); Python only for the decisions. This reverses D12. |
| D18 | Documentation (user, 2026-10-08) | Inline comments in all executed code; documented field-selection logic; list of tested callers/VCFs and input requirements; tools and citations for every step; nf-core-style overview figure. |
| D19 | REF check tool (implementation of D17, 2026-10-08; replaced by D26) | REF alleles are checked with `bcftools norm --check-ref x` (its log is emitted by the patched nf-core module and evaluated by `check_reference.py ref-alleles`), not with `bcftools +fixref`: `+fixref` only handles biallelic SNVs. |
| D20 | Position-only match (implementation, 2026-10-08) | The `match`/`position` stratum counts positions: `bcftools isec -c all` still lists an extra allele at a shared position separately, its flags are combined per position. |
| D21 | Missing VAF (implementation of D16, 2026-10-08; revised by D27) | Every tumor record needs a VAF (FORMAT/AD, else FORMAT/AF + FORMAT/DP); otherwise the run stops. The allowed fraction is configurable (`HARMONIZE_VCF` `--max-missing-vaf`, default 0). Records with AD 0,0 (no reads) have an undefined VAF (`zero_depth`, counted and warned, not an error); `FORMAT/AF = nan` counts as missing AF. Both occur in the SEQC2 MuTect2 3.8 VCFs. |
| D22 | Sheet reading (user, 2026-10-13) | The samplesheet is read as in the nf-core template (`samplesheetToList` in `PIPELINE_INITIALISATION`, columns checked by nf-schema); the checks across rows run in `CHECK_SAMPLESHEET` (`check_sheets.py`) at the start of `VARIANTCOMP`. `--stratification` and `--known_variants` are plain file parameters, checked by `check_sheets.py` in the subworkflow that uses them (`DISCORDANCE`, `PREPARE_KNOWN_VARIANTS`). |
| D23 | Early PASS filter (user, 2026-10-13) | With `--filter_pass_only` (default) the records are split into PASS (`PASS`/`.`) and non-PASS right after sample selection. Only PASS records are compared; the non-PASS records run through the same normalization (`NORMALIZE_VCF`) so that DISCORDANCE classifies a missed variant called with another FILTER as `filtered_non_pass` (with its FILTER value). `present_but_filtered` now means a PASS call removed by the uniform filters or the region. Changes the normalization chain (CLAUDE.md). |
| D24 | Contig lengths and chrM (user, 2026-10-13) | Contig lengths are never fixed: a length mismatch of a nuclear contig is an error. chrM is dropped from all inputs and reported (records/intervals), its length is not checked; `--keep_mito` keeps it with normal checks. |
| D25 | Annotation (user, 2026-10-13) | Like variantinterpretation: Ensembl VEP (nf-core module, offline cache or `--vep_gtf`) + `bcftools +split-vep` for consequence, gene, impact and variant class; msisensor-pro `scan` for repeat regions (`--repeats_bed` reuses them). `--gene_bed` is removed; `--stratification` BEDs stay. If VEP is on, a cache or a GTF is required. |
| D26 | All input checks first (user, 2026-10-13) | `CHECK_INPUTS` runs right after `CHECK_SAMPLESHEET` and checks every VCF, alignment, BED, sheet and known-variant list with `check_reference.py` (pysam: headers, contigs, REF alleles against the FASTA); all other steps wait for it. Each check has a toggle (warnings instead of errors). Header-extraction processes (`bcftools view -h`, `bcftools query`, `samtools idxstats`) are removed; this revises D17/D19. Optional features need their inputs (`--callable_regions`, `--read_evidence` need BAMs; `--require_normal_coverage` normal BAMs), checked by `check_sheets.py`. Common and callable regions must be at least `--min_region_bp` long. Every step summary reaches MultiQC. |
| D27 | VAF source (user, 2026-10-13) | `--vaf_source` (`auto`, `af`, `ad`); `auto` (default) prefers FORMAT/AF (+ DP), else AD. Above `--max_missing_vaf` (default 0.05) missing VAFs, the VCF gets a warning and the header line `##variantcomp_vaf=unavailable`; its VAF strata and detection power become NA instead of stopping the run. |
| D28 | GIAB stratifications (user, 2026-10-14) | The GIAB genome stratifications v3.5 are part of the annotation: `--giab_stratifications` (default true) with `--giab_dir` (local folder, release layout) and `--giab_files` (`name=path`, GRCh38 defaults). They are checked in `CHECK_INPUTS` like every BED and give `giab_<name>` columns. If on without `--giab_dir`, the run stops at the start. |
| D29 | Offline resources (user, 2026-10-14) | No process downloads data. `scripts/download_resources.py` downloads the VEP cache, the GIAB stratifications and optionally a GENCODE GTF (checksums verified) and writes `resources.config` for `-c`. Downloads in Python with an explicit CA bundle (fixes curl exit 77 with unrelocated conda CA paths). |
| – | Open (2026-10-08) | The discordance class `outside_callable` can't occur, because both technologies of a sample share the comparison region; such calls are `present_but_filtered`. To discuss: drop the class or define it differently (e.g. call outside the *query* technology's own callable region). |
| D13 | Callers (implementation, 2026-10-08) | No caller profiles are needed (D10): the `caller` samplesheet column is informational. |
| D11 | Normalization order (implementation, 2026-10-07) | Records on non-reference contigs are dropped before the header is rewritten from the `.fai`. Records with a symbolic ALT only are dropped before `norm`, and remaining symbolic alleles of split records are dropped during harmonization. See CLAUDE.md for the full order. |

## 12. Additional feature (last on the roadmap): known/required variants

An optional extra input `--known_variants` (sheet: `sample`, `variants` [VCF
or TSV `chrom pos ref alt [gene label]`], `regions` [optional BED], `mode`).
It works in two modes:

- **`required`**: a small clinical checklist (<10 manually curated,
  clinically relevant variants). Every variant must be found in **every**
  technology. The output is a clear PASS/FAIL table per sample × technology ×
  variant. It shows the status (called PASS / called but filtered / not
  called), VAF, depth, read evidence (if BAMs are given), and whether the
  position is inside the target/callable regions. Any FAIL is reported
  prominently in the report and summary, and optionally as a non-zero
  pipeline status (`--fail_on_missing_required`).
- **`benchmark`**: an external truth set (e.g. SEQC2 v1.2.1 with
  High-Confidence_Regions_v1.2.bed, or the GIAB HG008 benchmark when the user
  sequences HG008). It is normalized like every other VCF, and each
  technology gets recall/precision against it, restricted to its confident
  regions ∩ common targets. This runs in addition to the comparison against
  the truth technology.

Planned as `feat-known-variants` → `0.10.0-dev`, before the HPC bundle
(`0.11.0-dev`).
