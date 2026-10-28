# Roadmap: further improvements

This file collects work that is planned or proposed after v0.16.0-dev. The implemented
versions and agreed decisions are in [`pipeline-plan.md`](pipeline-plan.md) (§9, §11).
Items are grouped by topic, roughly in priority order; nothing here is agreed yet.

## 1. Validation

- Run the real use case: DRAGEN v4.4.4 VCFs, 8 samples × 3 technologies, tumor/normal
  paired. Add DRAGEN to the tested inputs in `docs/usage.md`.
- Full SEQC2 WES + WGS on the HPC (T2, `-profile test_full`, Singularity) and comparison
  with the published SEQC2 results.
- GIAB HG008 tumor/normal truth set as a second external benchmark (`--known_variants`,
  mode `benchmark`).

## 2. Report and statistics precision

- Exact precision in the strata and gene tables: they count the TP rows of the variant
  table, which are in the truth representation; with vcfeval matching one truth record to
  several query records, their precision differs from the metrics. Add the query-side
  matches (vcfeval `tp.vcf`) as rows of their own (documented limitation in
  `docs/output.md`).
- Read evidence of indels near read ends: a read counts only if it contains the whole
  haplotype window (`--evidence_indel_flank` bases around the indel), so the alt reads of
  long indels are a lower bound. Local realignment of the reads against the REF and ALT
  haplotypes would count partial overlaps too.
- (Done in 0.16.0-dev: datavzrd tables of all variant tables, MultiQC as the main report.)

## 3. `all_pairs` comparison

- Concordance matrix of all technologies (Jaccard index) and an UpSet plot of the shared
  and private calls per sample.

## 4. Statistics

- Confidence intervals for recall and precision (Wilson score or bootstrap over variants).
- Discordance per gene normalized by the callable bases of the gene.
- Detection power that takes tumor purity and ploidy into account.

## 5. Read-level features of false positives

- Strand bias, mapping and base quality, read position and soft clips at the discordant
  positions; clustering of FP calls (e.g. near indels or in low-complexity sequence).

## 6. Caller support

- Harmonization profiles for callers without `FORMAT/AF` or `FORMAT/AD`: Strelka2 (tier
  counts `AU/CU/GU/TU`, `TAR/TIR`), VarScan2, and tests with GATK Mutect2 4.x.
- DRAGEN-specific mapping of FILTER and INFO values (e.g. for `filtered_non_pass`).

## 7. VEP extensions

- VEP plugins: REVEL, SpliceAI, AlphaMissense.
- ClinVar / COSMIC: IDs and significance come from the VEP cache since 0.16.0-dev; a
  local ClinVar/COSMIC VCF as a VEP custom annotation would also work with a GTF.
- HGVS notation as input for the required known variants.

## 8. Discordance reason `outside_callable`

- Open since v0.13: the class can't occur, because both technologies of a sample are
  compared in the same region. Decide whether to drop it for good or to define it per
  technology (a call outside the *query* technology's own callable region).

## 9. Scale

- Sharding by chromosome for large WGS cohorts (normalization, vcfeval, read evidence).
- Cached reference artefacts: msisensor-pro repeats (`--repeats_bed`), VEP cache index.
- Resource labels tuned with the trace reports of the HPC runs.

## 10. Infrastructure

- CI (GitHub / GitLab@HHU): `nf-core pipelines lint`, black, ruff, pytest and stub runs.
- Conda lock files for the dev environment and the local modules.
- Pinned Singularity images in the offline bundle (`scripts/make_offline_bundle.py`).
- pre-commit hooks: black, ruff, prettier.

## 11. Inputs

- UMI- and amplicon-specific settings (e.g. consensus reads, primer trimming checks).
- FASTQ entry point (out of scope for v1).
- Tumor-only mode (no normal sample, no germline reason).

## 12. Output schema and provenance

- nf-schema output definitions for every published file.
- RO-Crate of the run for provenance (`nextflow` RO-Crate plugin).
