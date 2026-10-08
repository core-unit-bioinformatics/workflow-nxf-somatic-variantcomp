# core-unit-bioinformatics/workflow-nxf-somatic-variantcomp

> [!WARNING]
> Development version. Not yet validated on a full dataset (see the roadmap in
> [`docs/pipeline-plan.md`](docs/pipeline-plan.md)).

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

## Pipeline summary

![Metro map of workflow-nxf-somatic-variantcomp](docs/images/variantcomp_metro_map.png)

The map runs from top to bottom. Each grey box is a subworkflow, each row a step that
every sample goes through (numbers as in the table below). Each coloured line is a use
case and stops at the steps that run in it:

- **green — always:** the VCFs of all technologies are checked, normalized, filtered
  identically and compared with the truth technology;
- **blue — targeted technologies:** target BEDs (WES, panels) define the common
  comparison region; `--regions_bed` / `--exclude_bed` refine it;
- **orange — with alignments (BAM/CRAM):** callable regions from coverage restrict the
  comparison, and read counts explain missed or extra calls;
- **pink — with `--gene_bed` / `--stratification`:** genes and genomic context
  (repeats, homopolymers, ...) of every variant and enrichment statistics;
- **yellow — `--known_variants`, mode `required`:** a clinical checklist, PASS/FAIL per
  sample, technology and variant;
- **red — `--known_variants`, mode `benchmark`:** recall and precision of every technology
  against an external truth set.

Every step runs one tool (mostly nf-core modules) or one script of `bin/`; the options of
every step are set in `conf/modules.config`.

<!-- steps:start (generated: scripts/make_metro_map.py --markdown) -->
| Step | What happens | Tools | Runs for |
| --- | --- | --- | --- |
| **PIPELINE_INITIALISATION, PREPARE_ALIGNMENTS** | _check the inputs_ | | |
| 1.1 | **Check the samplesheet, stratification and known-variant sheets.** Checks the per-column formats (nf-schema) and the rows across the sheet: one truth technology per sample, one target BED and library type per technology, existing files. All problems are reported at once with line numbers; VCF and alignment indices are looked up next to the files. | `nf-schema, check_sheets.py` | all use cases |
| 1.2 | **Index alignments without an index.** Alignments without a `.bai`/`.csi`/`.crai` next to them are indexed once; they are needed for coverage and read evidence. | `samtools index` | With alignments (BAM/CRAM) |
| **PREPARE_REFERENCE** | _index the reference genome_ | | |
| 2.1 | **Unpack and index the reference FASTA.** A gzipped FASTA is unpacked and indexed; the `.fai` (contig names and lengths) is the reference every input is checked against. | `gunzip, samtools faidx` | Always |
| 2.2 | **Genome BED: one interval per reference contig.** One interval per contig, used as comparison region when no technology has a target BED (WGS only) and to drop records on other contigs. | `gawk` | Always |
| **PREPARE_REGIONS** | _build the common comparison region_ | | |
| 3.1 | **Check target BEDs against the reference, rename contigs.** Contig names that differ only by the `chr` prefix or `MT`/`chrM` are renamed, intervals on unknown contigs are dropped and reported; intervals beyond a contig end stop the run (other genome build). | `check_reference.py` | Targeted technologies |
| 3.2 | **Sort the BEDs in reference order.** All BEDs are sorted in the contig order of the reference, which the following bedtools steps need. | `bedtools sort` | Targeted technologies |
| 3.3 | **Keep bases targeted by every technology.** Only bases covered by the target BED of every targeted technology are kept (WGS counts as genome-wide), so every technology could have called a variant there. | `bedtools multiinter, gawk, merge` | Targeted technologies |
| 3.4 | **Restrict to --regions_bed, remove --exclude_bed (optional).** Optional: the region is restricted to `--regions_bed` and `--exclude_bed` is removed (e.g. difficult regions). | `bedtools intersect, subtract` | Targeted technologies |
| 3.5 | **Common region (genome if no target BEDs) and its size.** The final common comparison region (`regions/common_targets.bed`), merged; its size and the share of each technology's targets are summarized. | `bedtools merge, region_summary.py` | Targeted technologies, Always |
| **COVERAGE** | _find the bases callable in all technologies of a sample_ | | |
| 4.1 | **Check alignment contigs against the reference.** Contig names and lengths of every alignment are compared with the reference (same rules as for BEDs); the common region is translated to the alignment's contig names. | `samtools idxstats, check_reference.py` | With alignments (BAM/CRAM) |
| 4.2 | **Depth per base and per common-region interval.** mosdepth classifies every base as below or above `--min_coverage` and computes the mean depth per common-region interval (duplicates are kept for amplicon libraries). | `mosdepth` | With alignments (BAM/CRAM) |
| 4.3 | **Callable bases per alignment (depth >= --min_coverage).** The callable bases of each alignment (base mode; region mode uses the mean depth per interval), renamed to reference contigs and sorted. | `mosdepth_callable.py, bedtools sort` | With alignments (BAM/CRAM) |
| 4.4 | **Bases callable in all alignments of a sample (in common region).** Per sample: bases callable in the tumor (and normal) alignment of every technology, inside the common region. Variants are only compared there, so missing coverage in one technology does not count as a miss. | `bedtools multiinter, gawk, intersect` | With alignments (BAM/CRAM) |
| **PREPROCESS_VCF** | _normalize every VCF in the same way_ | | |
| 5.1 | **Check VCF contigs, contig lengths and sample names.** Contig names and lengths of the VCF header and records are checked against the reference, and the tumor/normal sample columns must exist. Safe fixes are prepared, a different build stops the run. | `bcftools view -h, query, check_reference.py` | Always |
| 5.2 | **Keep tumor (+ normal) columns, recompress with BGZF.** Only the tumor (and normal) columns are kept, in this order, and the VCF is written as bgzipped VCF with index (inputs may be plain gzip). | `bcftools view` | Always |
| 5.3 | **Rename contigs to the reference (chr prefix, MT/chrM).** Contigs are renamed to the reference names (`1` -> `chr1`, `MT` -> `chrM`). | `bcftools annotate` | Always |
| 5.4 | **Remove records on other contigs and gVCF blocks.** Records on contigs absent from the reference and gVCF reference blocks (only ALT `<NON_REF>`/`<*>`) are removed. | `bcftools view` | Always |
| 5.5 | **Rewrite the contig header from the reference.** Contig header lines are rewritten from the reference index (fixes wrong lengths, e.g. chrM, and missing contig lines). | `bcftools reheader` | Always |
| 5.6 | **Split multiallelics, left-align, split MNVs, deduplicate.** Multiallelic records are split, indels left-aligned, MNVs split into SNVs, records whose REF differs from the FASTA removed and counted, duplicates removed. Same representation for every caller. | `bcftools norm` | Always |
| 5.7 | **Stop if too many REF alleles differ (other genome build).** The fraction of REF mismatches from step 5.6 is checked: above `--max_ref_mismatch` the VCF is on another reference and the run stops. | `check_reference.py` | Always |
| 5.8 | **Same VAF, depth and alt reads for every caller (from AD).** VAF, depth and alt reads are computed the same way for every caller (from FORMAT/AD, else AF + DP; documented in docs/usage.md) and stored as `VC_*` INFO tags; GT is set to tumor 0/1, normal 0/0. Every tumor record needs a VAF. | `harmonize_vcf.py` | Always |
| 5.9 | **Uniform filters; keep the common (or callable) region.** The same filters for every technology (`--filter_pass_only`, `--min_vaf`, `--min_depth`, `--min_alt_reads`) and restriction to the comparison region: the common region, or with alignments the sample's callable region. The unfiltered VCF is kept for the discordance step. | `bcftools view` | Always, With alignments (BAM/CRAM) |
| **BENCHMARK** | _compare every technology with the truth technology_ | | |
| 6.1 | **Match calls haplotype-aware: TP, FN, FP.** Each technology's calls are matched with the truth technology's calls of the same sample, haplotype-aware, so different representations of one variant still match (TP, FN = missed, FP = only in the query). | `rtg vcfeval` | Always |
| 6.2 | **Exact-allele and same-position matches.** Record-by-record comparison by CHROM/POS/REF/ALT and by position only, to see how much discordance is representation or another allele. | `bcftools isec` | Always |
| 6.3 | **Recall, precision, F1 per sample (overall, type, VAF).** Recall, precision and F1 per sample and technology: overall, SNV vs. indel, by truth VAF (detection limit) and query VAF (low-VAF false positives), exact and position matches. | `benchmark_metrics.py` | Always |
| 6.4 | **Combine samples: summed (micro) and mean (macro).** Over all samples: from the summed counts (micro) and as mean +- sd of the per-sample values (macro). Main result: `benchmark/metrics.aggregated.tsv`. | `aggregate_metrics.py` | Always |
| **DISCORDANCE** | _explain why variants differ_ | | |
| 7.1 | **Positions of all compared variants of a sample.** All positions of compared variants of a sample are collected, so every alignment is looked at the same sites. | `bcftools query` | With alignments (BAM/CRAM) |
| 7.2 | **Read counts of every alignment at these positions.** Allele counts of every tumor and normal alignment at these positions (MQ >= 20, BQ >= 13), normalized like the calls. | `bcftools mpileup + call, bcftools norm` | With alignments (BAM/CRAM) |
| 7.3 | **Reason for every FN/FP; detection power of the query.** Each FN/FP gets a reason from the other technology's view: called but filtered, germline (alt reads in the normal), low coverage, alt reads below the caller's threshold, no alt reads, or absent (no alignment). For truth variants, the binomial detection power of the query gives the power-adjusted recall. | `classify_discordance.py` | Always, With alignments (BAM/CRAM) |
| 7.4 | **Check and sort gene / stratification BEDs.** Gene and stratification BEDs (repeats, homopolymers, ...) are checked against the reference and sorted. | `check_reference.py, bedtools sort` | With --gene_bed / --stratification |
| 7.5 | **Add genes and genomic context to every variant.** Every compared variant gets its gene(s) and a yes/no per stratification; result: `discordance/variants.annotated.tsv`. | `bedtools map, intersect -c, annotate_variants.py` | With --gene_bed / --stratification |
| **PREPARE_KNOWN_VARIANTS, KNOWN_VARIANTS** | _check known variants_ | | |
| 8.1 | **Read the lists (CSV/TSV/VCF), check contigs and positions.** Known variants from CSV, TSV or VCF lists are mapped to reference contigs, positions and alleles checked, and written as a sorted VCF. | `known_variants.py prepare` | --known_variants, mode required, --known_variants, mode benchmark |
| 8.2 | **Normalize like the calls; stop on any REF mismatch.** The list is normalized like the calls (split, left-align, atomize); any REF allele that differs from the FASTA stops the run. | `bcftools norm, check_reference.py` | --known_variants, mode required, --known_variants, mode benchmark |
| 8.3 | **PASS/FAIL of each variant per sample and technology.** Required mode: each listed variant must be found in every technology of the sample: PASS, or FAIL as filtered, outside the comparison region or not called (with read evidence if alignments are given). | `known_variants.py check` | With alignments (BAM/CRAM), --known_variants, mode required |
| 8.4 | **Checklist summary over all samples.** Required mode: one table over all lists, samples and technologies, a warning per failed variant, and optionally a stop (`--fail_on_missing_required`). | `known_summary.py` | --known_variants, mode required |
| 8.5 | **Recall/precision of every technology vs the truth set.** Benchmark mode: every technology is compared with the external truth set inside its confident regions, with the same metrics as the main benchmark. | `rtg vcfeval, bcftools isec, benchmark_metrics.py` | --known_variants, mode benchmark |
| **REPORT** | _summarize_ | | |
| 9.1 | **Enrichment statistics, gene table and plots.** Are discordant variants enriched for indels, low VAF, stratifications or genes compared with concordant ones (Fisher's exact test, Benjamini-Hochberg)? Plus gene table, class counts and plots. | `variantcomp_report.py` | Always, With --gene_bed / --stratification |
| 9.2 | **Interactive report of all results.** One report with metrics, discordance classes, recall by VAF, coverage (mosdepth), required variants and software versions. | `MultiQC` | all use cases |
<!-- steps:end -->

All tools and datasets, with their authors and publications, are listed in
[`CITATIONS.md`](CITATIONS.md). Map and table are generated by `scripts/make_metro_map.py`
(`--outdir docs/images`, `--markdown`).

Outputs are described in [`docs/output.md`](docs/output.md), parameters in
[`docs/usage.md`](docs/usage.md), HPC/offline use in [`docs/hpc.md`](docs/hpc.md) and testing
in [`docs/testing.md`](docs/testing.md).

## Usage

The samplesheet is a comma-separated file with one row per sample and technology,
described in [`docs/usage.md`](docs/usage.md).

```bash
nextflow run core-unit-bioinformatics/workflow-nxf-somatic-variantcomp \
   -profile <conda/singularity/.../institute> \
   --input samplesheet.csv \
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
- Test data (synthetic edge cases; SEQC2 HCC1395/HCC1395BL, WES + WGS) live in a separate
  repository and are prepared reproducibly with the scripts in `tests/`.
- Offline bundle for an HPC: `scripts/make_offline_bundle.py` (see `docs/hpc.md`).

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
