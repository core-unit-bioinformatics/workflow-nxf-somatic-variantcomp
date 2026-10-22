# Plan: report redesign and consistency checks (2026-10-22)

User request of 2026-10-22; decisions D33 and D34 in `docs/pipeline-plan.md` §11.

## Decisions (user, 2026-10-22)

| Topic | Decision |
| --- | --- |
| Report tool | A **narrative HTML report** (one offline file, written by a Python script) is the entry point: Summary, Introduction, Results, Methods. Browsable, linked tables (all variants, truth variants, each discordance class, regions, metrics) are a **datavzrd** report linked from it. **MultiQC** stays as the QC appendix (linked from Methods). |
| Consistency checks | Checks of counts and sums at the end of the pipeline. A failed check **stops the run**, after the report (with the failed check marked) has been written. |

## Report structure (ordered by importance)

1. **Header:** pipeline version, date, samples, technologies, truth technology; the
   result of the consistency checks (all passed / which failed).
2. **Summary:** key numbers; **recall/precision scatter** (one dot per sample, colour =
   query technology); truth variants found per technology; discordance classes per
   comparison; links to the tables.
3. **Introduction:** question, data (samplesheet table: samples, technologies, library
   types, callers), truth technology, what recall/precision/F1 and power-adjusted recall
   mean here.
4. **Results**
   1. *Comparison regions:* step diagram (genome → each target BED → intersection →
      `--regions_bed` → `--exclude_bed` → callable per sample), bar plot of bp per step,
      table (bp, fraction of the genome, of each target, of the common region).
   2. *Input variant sets:* per VCF: records, SNV/indel, PASS/non-PASS.
   3. *Filtering of the variant sets:* per VCF the records removed at each step and why
      (non-PASS, mitochondrial/unknown contig, symbolic/gVCF, REF mismatch, duplicate,
      uniform filters VAF/depth/alt reads, outside the comparison region); stacked bar
      plot and table.
   4. *Comparison metrics:* overall, SNV/indel, VAF bins, power-adjusted; per pair.
   5. *Truth variants:* every compared truth variant, found / missed by each query
      technology (top rows here, all rows in datavzrd).
   6. *Discordant variants:* per class an explanation, the counts and the top variants;
      all in datavzrd.
   7. *Annotation:* consequence, impact and repeat context per sample and technology,
      for truth variants and for discordant variants; enrichment tests; genes.
   8. *Known variants:* required checklist (PASS/FAIL), external truth sets.
5. **Methods:** every step with tool and version, the parameters used, normalization
   (`bcftools norm`) and harmonization counts, coverage (mosdepth; common region and
   target regions), the consistency checks with expected/observed values, links to
   MultiQC and datavzrd.

Text for every section explains what was done and what each value and column means, for
a reader who knows molecular biology and bioinformatics but not the pipeline.

## New data and steps

| Step | What | Where |
| --- | --- | --- |
| R1 | **Variant accounting** per VCF: counts per step and per removal reason, by SNV/indel | `bin/variant_accounting.py`, module `VARIANT_ACCOUNTING` |
| R2 | **Region accounting**: bp after every region step | `region_summary.py steps`, PREPARE_REGIONS |
| R3 | **Coverage** in the common region and in each technology's targets (mean depth, % bases ≥ min coverage / 30X); MultiQC General Statistics | mosdepth `--thresholds`; second mosdepth run on the target BED; `bin/coverage_summary.py` |
| R4 | **Report tables**: truth variants × query technologies, per-class tables, annotation per sample/technology/truth/discordant, P/R scatter | `bin/variantcomp_report.py` |
| R5 | **Consistency checks** + gate after the report | `bin/consistency_checks.py`, modules `CONSISTENCY_CHECKS`, `CHECKS_GATE` |
| R6 | **Narrative HTML report** | `bin/report_html.py`, module `REPORT_HTML` |
| R7 | **datavzrd** tables (config written by a script) | `bin/datavzrd_config.py`, nf-core module `datavzrd` |
| R8 | Docs, CLAUDE.md design, CITATIONS (datavzrd: Wiegand et al. 2025), CHANGELOG | |
| R9 | Tests: pytest per script, nf-tests, snapshots (at the end) | |

## Consistency checks (R5)

- Per VCF: input records = PASS + non-PASS (+ removed before the split); every later step:
  in = out + removed (per reason).
- Per pair: TP_baseline + FN = truth records compared; TP_call + FP = query records
  compared (from the compared VCFs).
- Aggregated metrics (micro) = sum of the per-pair counts.
- Discordance: FN + FP rows of `variants.annotated.tsv` = per-pair FN + FP; every FN/FP
  row has a known class; the class counts add up to FN + FP.
- Regions: common ≤ each target; callable ≤ common; every step ≤ the step before.
- Known variants: passed + failed = required.

Output `report/consistency_checks.tsv`: check, scope, expected, observed, status.
