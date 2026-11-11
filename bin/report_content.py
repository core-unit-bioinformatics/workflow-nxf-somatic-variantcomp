#!/usr/bin/env python3
"""The content of the variantcomp report: MultiQC sections and the datavzrd table browser.

The MultiQC report is the main report of a run (decision D37). This script writes its
variantcomp sections as MultiQC custom content, ordered by importance by
assets/multiqc_config.yml: Summary (key numbers as cards, key findings), Introduction
(this run and the comparison methods as lists), recall and precision per sample, truth
variants found, known variants, why variants differ, the discordant variants in detail,
metrics per sample, compared region, filtering of the input VCFs, coverage and the
consistency checks. The QC sections of the tools (bcftools norm, mosdepth, VEP, input
checks) follow below them. Tables with variants in rows are not put into MultiQC (only
short top lists of at most MAX_ROWS rows); the sections link to the datavzrd tables.

Outputs (one call writes both)
  --multiqc   folder with MultiQC custom-content files:
              vc_<section>_mqc.html  HTML sections (cards, lists, small tables, plots). The
                                     plots are drawn with the Plotly library of the MultiQC
                                     report (vertical bars, samples on the x-axis; MultiQC's
                                     own bar plots are horizontal) and follow the toolbox
                                     (hide, highlight and rename samples)
              vc_<section>_mqc.json  native MultiQC sections (recall-precision scatter plot
                                     with one dataset per matching method; metrics table)
  --datavzrd  folder with the tables and the configuration (datavzrd.yaml) of the datavzrd
              report: linked tables of every variant, truth variant, discordant variant,
              gene, known variant, removed record, region step, coverage row and check, with
              a description at the top, coloured status columns, links to dbSNP, COSMIC,
              GeneCards and Ensembl, a detail view for the less important columns, and a page
              that shows the MultiQC report (multiqc_page.js)

The explanations of the discordance classes, filtering reasons and columns (CLASS_TEXT,
REASON_TEXT, COLUMN_TEXT) are used by both outputs and kept in line with docs/output.md.

INPUT (all tables are TSV with header; optional inputs may be left out)
  --per-pair, --metrics       benchmark/metrics/metrics.per_pair.tsv, metrics.aggregated.tsv
                              (stratum_type overall = RTG vcfeval; match/exact and
                              match/position = bcftools isec)
  --power                     benchmark/metrics/power.aggregated.tsv (optional)
  --variants                  discordance/variants.annotated.tsv (every compared variant)
  --truth-variants            report/<prefix>.truth_variants.tsv (variantcomp_report.py)
  --classes, --enrichment, --genes, --annotation-counts   report tables
  --accounting, --removed     vcf_preprocessed/accounting/*.accounting.tsv (and the
                              *.checks.tsv of the same folder, which are skipped) /
                              *.removed.tsv
  --region-steps, --targets   regions/common_targets.steps.tsv / .summary.tsv
  --callable, --coverage      coverage/callable_summary.tsv, coverage/coverage.coverage.tsv
  --required-summary, --required-variants, --known-benchmark   known_variants/ (optional)
  --checks                    report/consistency_checks.tsv
  --igv-links                 report/igv/igv_links.tsv (--igv_report; optional): variant_id,
                              sample, part, igv_file; adds the column igv (link to the IGV
                              view of the variant) to every datavzrd table with a variant_id
                              and a link to the IGV overview to chapter 3
  --igv-link                  path of the IGV overview (report/igv/index.html) relative to the
                              MultiQC report
  --samplesheet               pipeline_info/samplesheet.validated.csv
  --params                    parameter table (group, parameter, value)
  --pipeline-version, --run-name, --start   run information
  --genome                    genome build of the links to Ensembl (GRCh38 or GRCh37)
  --fdr                       q-value threshold of the enrichment tests (marks rows)
  --filter-pass-only, --min-vaf, --min-depth, --min-alt-reads, --callable-regions
                              filter settings of the run (shown in the Introduction)
  --tables-link               path of the datavzrd report relative to the MultiQC report
  --multiqc-link              path of the MultiQC report relative to a datavzrd view page
OUTPUT
  --multiqc   folder of the MultiQC custom-content files (see above)
  --datavzrd  <folder>/datavzrd.yaml and the tables it names (paths are file names:
              datavzrd runs next to the tables)

Coordinates are 1-based (VCF) everywhere. Sizes are shown in Mbp (kbp, bp for small ones).

Why a script: the report combines the tables of all steps with explanatory text and
plots; no tool writes such a report from these tables. MultiQC renders the sections and
datavzrd the tables; MultiQC's own plot types were checked first: its bar plots are
horizontal only (MultiQC 1.35), so vertical bars with the samples on the x-axis are drawn
with the Plotly library that MultiQC ships in every report.
"""

import argparse
import html
import json
import os
import re
import sys
from datetime import datetime

import pandas as pd

# ----------------------------------------------------------------------------- texts

# main reasons (classify_discordance.py): exactly one per discordant variant; they
# describe the call and the reads of the technology that misses the variant
CLASS_ORDER = [
    "filtered_non_pass",
    "present_but_filtered",
    "low_coverage",
    "subthreshold_evidence",
    "no_evidence",
    "absent",
]
CLASS_TEXT = {
    "filtered_non_pass": (
        "The other technology called the variant, but its caller set a FILTER other than "
        "PASS. Only PASS calls are compared (--filter_pass_only), so the call is missing from "
        "the comparison. The FILTER value is in the column other_filter. This points to "
        "differences in caller filtering (e.g. strand bias, low VAF, germline risk) rather "
        "than in the data."
    ),
    "present_but_filtered": (
        "The other technology called the variant with FILTER PASS, but the call was removed "
        "by the uniform filters of the pipeline (--min_vaf, --min_depth, --min_alt_reads) or "
        "lies outside the comparison region. The variant is in the data of both technologies; "
        "the difference comes from a threshold of this comparison."
    ),
    "low_coverage": (
        "Not called. The tumor alignment of the technology that misses the variant has a depth "
        "below --min_coverage_tumor at the variant: it could not be detected reliably. Typical "
        "for positions at the edge of or outside a capture target."
    ),
    "subthreshold_evidence": (
        "Not called, although the tumor alignment of the technology that misses the variant "
        "has enough depth and at least one read with the alternative allele. Typical for "
        "low-VAF variants near the detection limit of the technology."
    ),
    "no_evidence": (
        "Not called; the tumor alignment of the technology that misses the variant has enough "
        "depth but no read with the alternative allele (MQ >= --evidence_min_mapq, SNV bases "
        "with BQ >= --evidence_min_baseq; indels by their haplotype sequence). For a false "
        "negative this suggests a variant below the detection limit or a sampling difference; "
        "for a false positive an artefact of the technology that called it."
    ),
    "absent": (
        "Not called and no read evidence available to look further (no alignment given, or "
        "--read_evidence false)."
    ),
}
# further reasons: in addition to the main reason (germline from the reads of the normal,
# the others from the annotation; every stratification is a reason of its own name)
FURTHER_ORDER = ["germline_in_normal", "indel", "low_vaf", "in_repeat"]
FURTHER_TEXT = {
    "germline_in_normal": (
        "The normal alignments of the sample show the alternative allele (VAF >= "
        "--germline_vaf with at least --evidence_min_alt_reads reads): probably a germline "
        "variant (or a contaminated normal), which one caller rejected and the other not."
    ),
    "indel": "The variant is an insertion, deletion or complex variant (harder to call and to align).",
    "low_vaf": (
        "The VAF of the call that exists (truth call of an FN, query call of an FP) is below "
        "--low_vaf: near the detection limit."
    ),
    "in_repeat": "The variant lies in a repeat (microsatellite) region of the reference (msisensor-pro).",
}
STRATUM_TEXT = (
    "The variant overlaps the stratification of this name (GIAB genome stratifications are "
    "named giab_<name>, e.g. giab_homopolymer; others come from --stratification)."
)
# short labels of the classes for legends and lists
CLASS_LABEL = {
    "filtered_non_pass": "called, FILTER not PASS",
    "present_but_filtered": "PASS, removed by uniform filters/region",
    "germline_in_normal": "germline (alt reads in the normal)",
    "indel": "indel",
    "low_vaf": "low VAF",
    "in_repeat": "in a repeat",
    "low_coverage": "low coverage",
    "subthreshold_evidence": "alt reads, not called",
    "no_evidence": "no alt reads",
    "absent": "absent (no alignment)",
}
CLASS_COLOR = {
    "filtered_non_pass": "#4e79a7",
    "present_but_filtered": "#76b7b2",
    "germline_in_normal": "#b07aa1",
    "low_coverage": "#9c755f",
    "subthreshold_evidence": "#f28e2b",
    "no_evidence": "#e15759",
    "absent": "#bab0ac",
    "concordant": "#59a14f",
    "indel": "#edc948",
    "low_vaf": "#ff9da7",
    "in_repeat": "#86bcb6",
}
# colour of a stratification reason (giab_* and --stratification)
STRATUM_COLOR = "#a0cbe8"
STATUS_COLOR = {"TP": "#59a14f", "FN": "#f28e2b", "FP": "#b07aa1"}
REASON_TEXT = {
    "filter_not_pass": "FILTER not PASS ('.' counts as PASS); only PASS calls are compared",
    "mitochondrial": "on the mitochondrial contig, which is not compared (--keep_mito false)",
    "contig_not_compared": "on a contig not in the compared reference (decoy, alt, unplaced)",
    "gvcf_block": "gVCF reference block (only ALT <NON_REF> or <*>), not a variant",
    "ref_mismatch": "REF allele differs from the reference FASTA (bcftools norm --check-ref x)",
    "duplicate": "exact duplicate after normalization (bcftools norm --rm-dup exact)",
    "split": "added: a multiallelic record split into one record per ALT, MNVs atomized",
    "symbolic_allele": "symbolic ALT allele (<DEL>, <*>, ...) left after splitting",
    "overlap_allele": "overlapping allele '*' left after splitting: the position is covered by a variant of another record (e.g. a spanning deletion), not a variant of its own",
    "outside_region": "outside the comparison region (common targets / callable bases)",
    "vaf_below_min": "tumor VAF below --min_vaf (or no VAF)",
    "depth_below_min": "tumor depth below --min_depth",
    "alt_reads_below_min": "tumor alt reads below --min_alt_reads",
}
COLUMN_TEXT = {
    "variant": "gene and protein change (HGVS.p; HGVS.c or the consequence if there is none)",
    "sample": "sample (one tumor/normal pair)",
    "truth": "truth (baseline) technology",
    "query": "query technology compared with the truth",
    "status": "TP: found by both; FN: truth variant missed by the query (in the truth only); FP: query call without truth variant (in the query only)",
    "class": "main reason: the call and reads of the technology that misses the variant (exactly one; concordant for TP)",
    "reasons": "all reasons of a variant: the main reason first (concordant for a TP), then germline_in_normal and the context (indel, low_vaf, in_repeat, stratifications)",
    "other_tumor_other": "reads with another allele at the variant in that alignment (e.g. another indel in a repeat)",
    "locus": "position as chrom:pos (1-based), linked to the Ensembl genome browser",
    "chrom": "contig",
    "pos": "position (1-based)",
    "ref": "reference allele",
    "alt": "alternative allele",
    "type": "SNV, MNV, INS, DEL or COMPLEX",
    "length": "len(ALT) - len(REF)",
    "truth_vaf": "tumor VAF of the truth call",
    "query_vaf": "tumor VAF of the query call",
    "truth_dp": "tumor depth of the truth call",
    "query_dp": "tumor depth of the query call",
    "query_filter": "FILTER of the query call",
    "other_filter": "FILTER of the call in the technology that misses the variant",
    "other_tumor_dp": "tumor depth in the alignment of the technology that misses the variant",
    "other_tumor_alt": "alt reads in that alignment",
    "normal_vaf": "VAF in the normal alignments (germline evidence)",
    "query_power": "probability that the query detects the variant (binomial model)",
    "detectable": "query_power >= --detection_power",
    "gene": "gene symbol (Ensembl VEP), linked to GeneCards and Ensembl",
    "hgvsp": "protein change (HGVS.p, VEP --hgvs), of the transcript of the most severe consequence",
    "hgvsc": "coding change (HGVS.c, VEP --hgvs)",
    "transcript": "transcript of the most severe consequence (VEP Feature)",
    "dbsnp": "dbSNP IDs of known variants at this allele (VEP cache, --check_existing), linked to dbSNP",
    "cosmic": "COSMIC IDs (VEP cache), linked to COSMIC",
    "clin_sig": "ClinVar significance (VEP cache)",
    "consequence": "most severe consequence (VEP, Sequence Ontology term)",
    "consequence_group": "synonymous, non_synonymous or other (from the consequence)",
    "impact": "VEP impact: HIGH, MODERATE, LOW, MODIFIER",
    "biotype": "biotype of the transcript (VEP)",
    "in_repeat": "inside a repeat (microsatellite) region of the reference",
    "found_by": "number of query technologies that found the truth variant",
    "variant_id": "sample:chrom:pos:ref>alt (link key between the tables)",
    "igv": "IGV view of the variant (igv-reports): reads of the tumor and normal alignments and the normalized calls of every technology of the sample (--igv_report)",
    "tp_baseline": "truth variants found by the query",
    "fn": "truth variants missed by the query",
    "tp_call": "query calls matching a truth variant",
    "fp": "query calls without a truth variant",
    "recall": "tp_baseline / (tp_baseline + fn): fraction of the truth variants found",
    "precision": "tp_call / (tp_call + fp): fraction of the query calls that are true",
    "f1": "harmonic mean of recall and precision",
}
STEP_TEXT = {
    "genome": "compared genome",
    "targets_intersection": "intersection of all target BEDs",
    "regions_bed": "within --regions_bed",
    "exclude_bed": "without --exclude_bed",
}
# matching methods: key, stratum_type, stratum, label, explanation
METHODS = [
    (
        "vcfeval",
        "overall",
        "all",
        "RTG vcfeval",
        "haplotype-aware matching (rtg vcfeval --squash-ploidy): different representations "
        "of the same variant (e.g. an MNV and two SNVs, shifted indels) match. Main result.",
    ),
    (
        "exact",
        "match",
        "exact",
        "isec: same allele",
        "bcftools isec -c none: a match needs the same CHROM, POS, REF and ALT after "
        "normalization.",
    ),
    (
        "position",
        "match",
        "position",
        "isec: same position",
        "bcftools isec -c all: a match needs the same CHROM and POS only (any allele); "
        "counted once per position.",
    ),
]
# colours of the outcomes of the input records (filtering plot): compared green, every
# reason its own colour
REASON_COLOR = {
    "compared": "#59a14f",
    "filter_not_pass": "#edc948",
    "outside_region": "#4e79a7",
    "gvcf_block": "#b07aa1",
    "mitochondrial": "#9c755f",
    "contig_not_compared": "#bab0ac",
    "ref_mismatch": "#e15759",
    "duplicate": "#76b7b2",
    "symbolic_allele": "#ff9da7",
    "overlap_allele": "#f28e2b",
    "vaf_below_min": "#8cd17d",
    "depth_below_min": "#86bcb6",
    "alt_reads_below_min": "#d4a6c8",
}
PALETTE = ["#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#76b7b2", "#edc948", "#b07aa1", "#9c755f"]
# at most this many rows of a variant (or gene) list in MultiQC; the rest is in datavzrd
MAX_ROWS = 15
# VAF bins of the characterization of discordant variants
VAF_BINS = [0, 0.05, 0.1, 0.2, 0.5, 1.0001]
VAF_LABELS = ["<5 %", "5-10 %", "10-20 %", "20-50 %", ">=50 %"]

# ----------------------------------------------------------------------------- helpers


def read(path, sep="\t"):
    """A table as a data frame ('NA' as missing); an empty frame if path is None."""
    if not path:
        return pd.DataFrame()
    return pd.read_csv(path, sep=sep, na_values=["NA"], keep_default_na=False, low_memory=False)


def read_many(paths, required=None):
    """Concatenation of several tables (empty frame if none); with required, only tables
    that have all these columns (e.g. accounting tables among the checks tables)."""
    frames = [read(p) for p in paths or []]
    frames = [f for f in frames if not f.empty and set(required or []) <= set(f.columns)]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def esc(value):
    """HTML-escaped text."""
    return html.escape(str(value))


def fmt_int(value):
    """Integer with thousands separators ('NA' if missing)."""
    return "NA" if pd.isna(value) else f"{int(value):,}"


def fmt_num(value, digits=3):
    """A metric value with the given decimals ('NA' if missing)."""
    return "NA" if value is None or pd.isna(value) else f"{float(value):.{digits}f}"


def fmt_pct(part, whole, digits=1):
    """part / whole as a percentage string ('NA' if whole is 0)."""
    return f"{100 * part / whole:.{digits}f} %" if whole else "NA"


def fmt_bp(value):
    """A size in Mbp (>= 1 Mbp), kbp (>= 1 kbp) or bp."""
    value = int(value)
    if value >= 1e6:
        return f"{value / 1e6:,.2f} Mbp"
    return f"{value / 1e3:.1f} kbp" if value >= 1e3 else f"{value} bp"


def step_label(name):
    """Readable label of a region step (names like 2_targets_intersection)."""
    return STEP_TEXT.get(re.sub(r"^\d+_", "", str(name)), str(name))


def value_of(row, column):
    """A cell as text; '' for a missing column, NA or '.'."""
    value = row.get(column)
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return ""
    value = str(value)
    return "" if value in (".", "NA") else value


def variant_label(row):
    """Readable name of a variant: gene + protein change (HGVS.p), else gene + HGVS.c, else
    gene + consequence; chrom:pos ref>alt without a gene."""
    gene = value_of(row, "gene")
    change = (
        value_of(row, "hgvsp")
        or value_of(row, "hgvsc")
        or value_of(row, "consequence").split("&")[0]
    )
    if not gene:
        return f"{row['chrom']}:{row['pos']} {row['ref']}>{row['alt']}"
    return f"{gene} {change}".strip()


def add_variant_columns(df):
    """variant_id (sample:chrom:pos:ref>alt; the link key between the tables), variant (gene
    + protein change) and locus (chrom:pos for the genome browser links)."""
    if df.empty or not {"chrom", "pos", "ref", "alt"} <= set(df.columns):
        return df
    df = df.copy()
    site = df["chrom"].astype(str) + ":" + df["pos"].astype(str)
    if "sample" in df.columns:
        ids = df["sample"].astype(str) + ":" + site + ":" + df["ref"].astype(str)
        df["variant_id"] = ids + ">" + df["alt"].astype(str)
    df["variant"] = [variant_label(r) for r in df.to_dict("records")]
    df["locus"] = site
    return df


def order_columns(df, front, back=()):
    """The columns of df: front first (those present), back last, the rest in between."""
    # dict.fromkeys: every column once, in the order of its first mention
    first = [c for c in dict.fromkeys(front) if c in df.columns]
    last = [c for c in back if c in df.columns and c not in first]
    middle = [c for c in df.columns if c not in first and c not in last]
    return df[first + middle + last]


def overall_rows(df, stratum_type="overall", stratum="all"):
    """The rows of one stratum (default: vcfeval overall) of a metrics table."""
    if df.empty:
        return df
    return df[(df["stratum_type"] == stratum_type) & (df["stratum"].astype(str) == stratum)]


def param(params, name, default="NA"):
    """Value of a pipeline parameter from the parameter table."""
    if params.empty or "parameter" not in params.columns:
        return default
    rows = params[params["parameter"] == name]
    return str(rows["value"].iloc[0]) if not rows.empty else default


def comparisons(per_pair):
    """The (truth, query) pairs of the run, sorted."""
    if per_pair.empty:
        return []
    return sorted({(t, q) for t, q in per_pair[["truth", "query"]].itertuples(index=False)})


def html_table(df, link=None, link_text=None, total=None, classes="table table-sm table-hover"):
    """A small HTML table (at most MAX_ROWS rows) with a link to the full table."""
    if df.empty:
        return "<p class='text-muted'>No rows.</p>"
    head = "".join(f"<th>{esc(c)}</th>" for c in df.columns)
    body = []
    for row in df.head(MAX_ROWS).itertuples(index=False):
        cells = []
        for value in row:
            text = (
                fmt_num(value, 4)
                if isinstance(value, float) and not float(value).is_integer()
                else (fmt_int(value) if isinstance(value, (int, float)) else str(value))
            )
            # pre-rendered HTML (badges) is passed through, everything else escaped
            cells.append(f"<td>{text if text.startswith('<span') else esc(text)}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    out = (
        f"<div class='table-responsive'><table class='{classes}'><thead><tr>{head}</tr></thead>"
        f"<tbody>{''.join(body)}</tbody></table></div>"
    )
    n = total if total is not None else len(df)
    if link:
        more = f"Showing {min(n, MAX_ROWS)} of {fmt_int(n)} rows. " if n > MAX_ROWS else ""
        out += (
            f"<p class='vc-note'>{more}<a href='{esc(link)}' target='_blank'>"
            f"{esc(link_text or 'All rows in the table browser')} &rarr;</a></p>"
        )
    return out


def venn2_svg(truth, query, n_truth_only, n_both, n_query_only, n_both_query=None):
    """SVG Venn diagram of the variants of two technologies of one sample: the truth
    technology only (FN), both (TP), the query technology only (FP). n_both_query: TP in
    the query representation, shown if it differs (one truth record can match several
    query records with vcfeval)."""
    left, right = STATUS_COLOR["FN"], STATUS_COLOR["FP"]
    both = f"{n_both}" if n_both_query in (None, n_both) else f"{n_both} ({n_both_query})"
    tip_both = f"in both: {n_both} truth variants found by {query}" + (
        "" if n_both_query in (None, n_both) else f" ({n_both_query} {query} calls)"
    )
    return (
        "<svg viewBox='0 0 300 190' width='280' role='img' style='font-family: inherit'>"
        f"<circle cx='115' cy='95' r='75' fill='{left}' fill-opacity='0.3' stroke='{left}' "
        f"stroke-width='2'><title>{esc(truth)}</title></circle>"
        f"<circle cx='185' cy='95' r='75' fill='{right}' fill-opacity='0.3' stroke='{right}' "
        f"stroke-width='2'><title>{esc(query)}</title></circle>"
        f"<text x='72' y='101' text-anchor='middle' font-size='18' font-weight='700' "
        f"fill='currentColor'>{n_truth_only}<title>only in {esc(truth)} (FN): {n_truth_only}"
        "</title></text>"
        f"<text x='150' y='101' text-anchor='middle' font-size='18' font-weight='700' "
        f"fill='currentColor'>{both}<title>{esc(tip_both)}</title></text>"
        f"<text x='228' y='101' text-anchor='middle' font-size='18' font-weight='700' "
        f"fill='currentColor'>{n_query_only}<title>only in {esc(query)} (FP): {n_query_only}"
        "</title></text>"
        f"<text x='72' y='185' text-anchor='middle' font-size='12' fill='{left}'>"
        f"{esc(truth)} only (FN)</text>"
        f"<text x='150' y='14' text-anchor='middle' font-size='12' fill='currentColor'>"
        "truth variants in both (TP)</text>"
        f"<text x='228' y='185' text-anchor='middle' font-size='12' fill='{right}'>"
        f"{esc(query)} only (FP)</text>"
        "</svg>"
    )


def badge(text, color):
    """A coloured label (HTML)."""
    return f"<span class='vc-badge' style='background:{color}'>{esc(text)}</span>"


# ----------------------------------------------------------------------------- MultiQC

# anchors of the MultiQC modules that follow a chapter (named and ordered in
# assets/multiqc_config.yml); styled as sections of that chapter
TOOL_MODULES = [
    "variantcomp_normalization",
    "variantcomp_harmonization",
    "vep",
    "mosdepth",
    "variantcomp_input_checks",
    "core-unit-bioinformatics-workflow-nxf-somatic-variantcomp-methods-description",
    "multiqc_software_versions",
    "core-unit-bioinformatics-workflow-nxf-somatic-variantcomp-summary",
]

CSS = """
<style>
.vc-cards { display: flex; flex-wrap: wrap; gap: 10px; margin: 8px 0 14px; }
.vc-card { background: var(--bs-tertiary-bg, #f3f6f9); border-radius: 8px; padding: 8px 14px;
           min-width: 130px; border-left: 4px solid var(--vc-accent, #4e79a7); }
.vc-card .v { font-size: 24px; font-weight: 700; line-height: 1.2; }
.vc-card .l { font-size: 12px; color: var(--bs-secondary-color, #5c6b7a); }
.vc-card .s { font-size: 11px; color: var(--bs-secondary-color, #5c6b7a); }
.vc-card.metric { min-width: 150px; }
.vc-group { margin: 6px 0 2px; font-weight: 600; }
.vc-note { font-size: 13px; color: var(--bs-secondary-color, #5c6b7a); }
.vc-badge { color: #fff; border-radius: 4px; padding: 1px 6px; font-size: 12px; white-space: nowrap; }
.vc-banner { padding: 8px 12px; border-radius: 6px; margin: 6px 0 12px; font-weight: 600; }
.vc-banner.ok { background: rgba(89,161,79,0.15); color: #2e7d32; }
.vc-banner.fail { background: rgba(225,87,89,0.15); color: #c62828; }
.vc-links a { margin: 0 6px 6px 0; }
.vc-flow { display: flex; flex-wrap: wrap; align-items: center; gap: 6px; margin: 10px 0; }
.vc-step { background: var(--bs-tertiary-bg, #f3f6f9); border-left: 4px solid #4e79a7;
           border-radius: 4px; padding: 4px 10px; font-size: 12px; }
.vc-step b { display: block; font-size: 13px; }
.vc-arrow { color: var(--bs-secondary-color, #5c6b7a); font-size: 18px; }
.vc-matrix td, .vc-matrix th { text-align: center; white-space: nowrap; }
.vc-matrix td.pass { background: rgba(89,161,79,0.25); }
.vc-matrix td.fail { background: rgba(225,87,89,0.3); font-weight: 600; }
.vc-plot-buttons .btn { margin: 0 4px 4px 0; }
.vc-hidden { display: none; }
.vc-venns { display: flex; flex-wrap: wrap; gap: 24px; margin: 6px 0; }
.vc-venn svg { max-width: 100%; height: auto; }
.vc-legend { font-size: 13px; } .vc-legend code { font-size: 11px; }
.vc-dot { display: inline-block; width: 11px; height: 11px; border-radius: 50%; margin-right: 6px; opacity: 0.8; }
/* the MultiQC welcome banner (tutorial video) is not needed in this report */
#mqc_welcome { display: none !important; }
/* the report as a book (decision D39): larger text, chapters with a heading line and an
   introduction box; the tool modules after a chapter (anchors in TOOL_MODULES) look like
   its sections, in the text and in the table of contents */
.mqc-module-section p, .mqc-module-section li, .mqc-section-description { font-size: 15px; line-height: 1.55; }
.mqc-module-section .vc-note { font-size: 14px; }
.mqc-module-section[data-module-anchor^="vc_ch_"] h2 { font-size: 2.1rem; font-weight: 700;
    border-bottom: 3px solid #4e79a7; padding-bottom: 6px; width: 100%; }
.mqc-module-section[data-module-anchor^="vc_ch_"] { margin-top: 3rem; }
.mqc-module-section[data-module-anchor^="vc_ch_"] .mqc-module-section-first > p,
.mqc-module-section[data-module-anchor^="vc_ch_"] .mqc-module-section-first > ul,
.mqc-module-section[data-module-anchor^="vc_ch_"] .mqc-module-section-first > ol {
    background: var(--bs-tertiary-bg, #f3f6f9); margin: 0; padding: 6px 16px 6px 16px; font-size: 15.5px; }
.mqc-module-section[data-module-anchor^="vc_ch_"] .mqc-module-section-first > ul,
.mqc-module-section[data-module-anchor^="vc_ch_"] .mqc-module-section-first > ol { padding-left: 36px; }
.mqc-module-section[data-module-anchor^="vc_ch_"] .mqc-module-section-first > p:first-of-type { padding-top: 12px; border-radius: 8px 8px 0 0; }
.mqc-module-section[data-module-anchor^="vc_ch_"] .mqc-module-section-first > :last-child { padding-bottom: 12px; border-radius: 0 0 8px 8px; margin-bottom: 12px; }
.mqc-section h3 { font-size: 1.45rem; margin-top: 1.6rem; }
/* filtering funnel: centred boxes, the changes of every step on the right */
.vc-funnel { display: flex; flex-direction: column; gap: 5px; margin: 8px 0 16px; }
.vc-frow { display: flex; flex-wrap: wrap; align-items: center; gap: 4px 12px; }
.vc-fbarwrap { flex: 0 0 62%; display: flex; justify-content: center; }
.vc-fbar { background: rgba(78,121,167,0.18); border: 1px solid rgba(78,121,167,0.6);
           border-radius: 6px; padding: 5px 10px; text-align: center; font-size: 13px; line-height: 1.35; }
.vc-fbar.final { background: rgba(89,161,79,0.25); border-color: #59a14f; }
.vc-fbar .vc-fn { font-size: 17px; font-weight: 700; margin-left: 4px; }
.vc-fbar .vc-fpct, .vc-fbar .vc-fkind { display: block; font-size: 11.5px; color: var(--bs-secondary-color, #5c6b7a); }
.vc-fchanges { flex: 1 1 30%; font-size: 12.5px; display: flex; flex-wrap: wrap; gap: 3px 8px; }
.vc-fremoved, .vc-fadded { border-left: 4px solid #e15759; padding-left: 5px; white-space: nowrap; }
.vc-fadded { border-color: #59a14f; }
.vc-fchanges code { font-size: 11.5px; }
.mqc-nav a.nav-l1[href^="#vc_ch_"] { font-weight: 700; }
TOOL_MODULE_CSS
</style>
"""
# tool modules: heading of the size of a section, indented in the table of contents
CSS = CSS.replace(
    "TOOL_MODULE_CSS",
    "\n".join(
        # plain concatenation: nf-core lint takes double braces for Jinja placeholders
        '.mqc-module-section[data-module-anchor="' + a + '"] h2 { font-size: 1.45rem; }\n'
        '.mqc-nav a.nav-l1[href="#' + a + '"] { padding-left: 1.4em; font-size: 0.92em; }'
        for a in TOOL_MODULES
    ),
)

# the plot library: draws every registered plot (window.vcPlotSpecs) with the Plotly of the
# MultiQC report, adds view buttons, and redraws on the toolbox events (hide, highlight and
# rename samples) with the same matching rules as MultiQC (substring or regex)
JS_LIB = r"""
<script>
(function () {
  if (window.vcLib) { return; }
  window.vcPlotSpecs = window.vcPlotSpecs || [];
  // does a name match a toolbox filter text (substring, or regex in regex mode)?
  function matches(name, f, regex) {
    if (!f || !name) { return false; }
    if (regex) { try { return new RegExp(f).test(name); } catch (e) { return false; } }
    return name.indexOf(f) > -1;
  }
  // is a sample hidden by the toolbox ("hide" or "show only" mode)?
  function hidden(name) {
    var fs = window.mqc_hide_f_texts || [];
    if (!name || fs.length === 0) { return false; }
    var m = fs.some(function (f) { return matches(name, f, window.mqc_hide_regex_mode); });
    return window.mqc_hide_mode === "show" ? !m : m;
  }
  // highlight colour of a sample (last matching filter), or null
  function highlight(name) {
    var fs = window.mqc_highlight_f_texts || [], cs = window.mqc_highlight_f_cols || [];
    var col = null;
    fs.forEach(function (f, i) { if (matches(name, f, window.mqc_highlight_regex_mode)) { col = cs[i]; } });
    return col;
  }
  // sample name after the toolbox renaming
  function rename(name) {
    var f = window.mqc_rename_f_texts || [], t = window.mqc_rename_t_texts || [];
    f.forEach(function (p, i) {
      if (!p) { return; }
      name = window.mqc_rename_regex_mode ? name.replace(new RegExp(p, "g"), t[i]) : name.split(p).join(t[i]);
    });
    return name;
  }
  var ARRAY_KEYS = ["x", "y", "text", "customdata", "hovertext", "s"];
  // a copy of a trace with the hidden samples removed and the highlights applied;
  // trace.s: sample of every point; trace.sample: sample of the whole trace (lines)
  function prepare(trace) {
    var t = JSON.parse(JSON.stringify(trace));
    if (t.sample !== undefined) {
      if (hidden(t.sample)) { return null; }
      var c = highlight(t.sample);
      if (c) { t.line = Object.assign({}, t.line || {}, { color: c, width: 4 }); }
      t.name = rename(t.name);
      delete t.sample;
    }
    if (t.s) {
      var keep = t.s.map(function (s) { return !s || !hidden(s); });
      ARRAY_KEYS.forEach(function (k) {
        if (Array.isArray(t[k])) { t[k] = t[k].filter(function (_, i) { return keep[i]; }); }
      });
      if (t.renameX) { t.x = t.x.map(rename); }
      var cols = t.s.map(function (s) { return s ? highlight(s) : null; });
      if (cols.some(function (c) { return c; })) {
        t.marker = Object.assign({}, t.marker || {});
        t.marker.line = { color: cols.map(function (c) { return c || "rgba(0,0,0,0)"; }),
                          width: cols.map(function (c) { return c ? 3 : 0; }) };
      }
      delete t.s;
    }
    delete t.renameX;
    return t;
  }
  // layout with the colours of the current MultiQC theme (light or dark)
  function themed(layout) {
    var fg = getComputedStyle(document.body).color;
    var l = JSON.parse(JSON.stringify(layout || {}));
    l.paper_bgcolor = "rgba(0,0,0,0)";
    l.plot_bgcolor = "rgba(0,0,0,0)";
    l.font = Object.assign({ color: fg, size: 12 }, l.font || {});
    Object.keys(l).filter(function (k) { return /^[xy]axis\d*$/.test(k); }).forEach(function (a) {
      l[a] = Object.assign({ gridcolor: "rgba(128,128,128,0.25)", zerolinecolor: "rgba(128,128,128,0.5)", automargin: true }, l[a]);
    });
    l.margin = Object.assign({ t: 30, r: 20, b: 40, l: 60 }, l.margin || {});
    l.legend = Object.assign({ orientation: "h", yanchor: "bottom", y: 1.02, x: 0 }, l.legend || {});
    // vertical toolbar on the right: it does not cover a legend of several rows
    l.modebar = { orientation: "v" };
    return l;
  }
  function draw(spec) {
    var el = document.getElementById(spec.id);
    if (!el || typeof Plotly === "undefined") { return; }
    var view = spec.views[spec.current || 0];
    var traces = view.traces.map(prepare).filter(function (t) { return t; });
    Plotly.react(el, traces, themed(view.layout),
      { displaylogo: false, responsive: true, toImageButtonOptions: { filename: spec.id, scale: 2 } });
  }
  // one button per view above the plot
  function buttons(spec) {
    var bar = document.getElementById(spec.id + "-buttons");
    if (!bar || spec.views.length < 2) { return; }
    spec.views.forEach(function (v, i) {
      var b = document.createElement("button");
      b.type = "button";
      b.className = "btn btn-sm btn-outline-secondary" + (i === 0 ? " active" : "");
      b.textContent = v.label;
      b.onclick = function () {
        spec.current = i;
        Array.prototype.forEach.call(bar.children, function (c) { c.classList.remove("active"); });
        b.classList.add("active");
        draw(spec);
      };
      bar.appendChild(b);
    });
  }
  // toggle groups of HTML elements (e.g. the cards of one matching method)
  window.vcToggle = function (group, key, button) {
    document.querySelectorAll("[data-vc-group='" + group + "']").forEach(function (el) {
      el.classList.toggle("vc-hidden", el.getAttribute("data-vc-key") !== key);
    });
    if (button) {
      Array.prototype.forEach.call(button.parentNode.children, function (c) { c.classList.remove("active"); });
      button.classList.add("active");
    }
  };
  function drawAll() { window.vcPlotSpecs.forEach(draw); }
  window.vcLib = { draw: draw };
  document.addEventListener("DOMContentLoaded", function () {
    // the software versions module has a fixed name in MultiQC: number it as section 6.4
    // of the chapter Methods (heading and table of contents)
    var name = "6.4 Software versions";
    var nav = document.querySelector('.mqc-nav a[href="#multiqc_software_versions"]');
    if (nav) { nav.textContent = name; }
    var head = document.getElementById("multiqc_software_versions");
    if (head) { head.textContent = name; }
    window.vcPlotSpecs.forEach(function (s) { buttons(s); draw(s); });
    if (typeof $ !== "undefined") {
      $(document).on("mqc_hidesamples mqc_highlights mqc_renamesamples", drawAll);
    }
    // redraw when the light/dark theme changes (colours of axes and text)
    new MutationObserver(drawAll).observe(document.documentElement, { attributes: true, attributeFilter: ["data-bs-theme"] });
  });
})();
</script>
"""


# ----------------------------------------------------------------------------- chapters

# The MultiQC report reads like a book (decision D39): six chapters, each one MultiQC
# module (custom content grouped by parent_id) with an introduction; chapters and
# sections are numbered, and MultiQC sorts the sections of a chapter by their numbered
# names. The tool modules (VEP, mosdepth), the step tables (MQC_TABLES), the methods
# description, the software versions and the workflow summary can't be nested in a
# custom-content module: assets/multiqc_config.yml gives them the number of their chapter
# and places them right after it (report_section_order).
CHAPTERS = {
    "summary": ("1", "Summary"),
    "known": ("2", "Known variants"),
    "discordance": ("3", "Discordance"),
    "filtering": ("4", "Variant filtering"),
    "regions": ("5", "Coverage and genomic regions"),
    "methods": ("6", "Methods and quality control"),
}
# section id -> (chapter, number of the section within its chapter)
SECTIONS = {
    "summary": ("summary", 1),
    "introduction": ("summary", 2),
    "recall_precision": ("summary", 3),
    "metrics": ("summary", 4),
    "metrics_table": ("summary", 5),
    "overlap": ("summary", 6),
    "truth_found": ("summary", 7),
    "known": ("known", 1),
    "why": ("discordance", 1),
    "characterization": ("discordance", 2),
    "filtering": ("filtering", 1),
    "regions": ("regions", 1),
    "coverage": ("regions", 2),
    "checks": ("methods", 1),
}


def section_title(section_id, name):
    """Numbered section name, e.g. '1.3 Recall and precision per sample'."""
    chapter, number = SECTIONS[section_id]
    return f"{CHAPTERS[chapter][0]}.{number} {name}"


class Mqc:
    """Writes the MultiQC custom-content files of the variantcomp sections."""

    def __init__(self, outdir, tables_link, intros=None):
        self.outdir = outdir
        self.tables = tables_link.rstrip("/")
        self.n_plots = 0
        # introduction (markdown) of every chapter, see chapter_intros()
        self.intros = intros or {}
        os.makedirs(outdir, exist_ok=True)

    def link(self, view):
        """Path of a datavzrd view relative to the MultiQC report."""
        return f"{self.tables}/{view}/index_1.html"

    def parent(self, section_id):
        """The chapter of a section as MultiQC custom-content keys. Every section of a
        chapter carries the same name and introduction: MultiQC takes them from whichever
        file of the chapter it reads first."""
        chapter, _ = SECTIONS[section_id]
        number, title = CHAPTERS[chapter]
        return {
            "parent_id": f"vc_ch_{chapter}",
            "parent_name": f"{number} {title}",
            "parent_description": self.intros.get(chapter, ""),
        }

    def html(self, section_id, name, description, body):
        """An HTML section: vc_<id>_mqc.html (name and description: markdown)."""
        keys = {
            "id": f"vc_{section_id}",
            "section_name": section_title(section_id, name),
            "description": description,
            **self.parent(section_id),
        }
        # YAML header in an HTML comment; JSON strings are valid YAML
        header = "<!--\n" + "".join(f"{k}: {json.dumps(v)}\n" for k, v in keys.items()) + "-->\n"
        with open(os.path.join(self.outdir, f"vc_{section_id}_mqc.html"), "w") as fh:
            fh.write(header + body + "\n")

    def json(self, section_id, content):
        """A native MultiQC section (plot or table): vc_<id>_mqc.json."""
        content = {
            "id": f"vc_{section_id}",
            **content,
            "section_name": section_title(section_id, content["section_name"]),
            **self.parent(section_id),
        }
        with open(os.path.join(self.outdir, f"vc_{section_id}_mqc.json"), "w") as fh:
            json.dump(content, fh, indent=1)

    def plot(self, views, height=420, caption=""):
        """HTML of a Plotly plot with one or more views (buttons switch between them).

        views: list of {label, traces, layout}; a trace may carry s (sample of every
        point) or sample (sample of the whole trace) for the toolbox.
        """
        views = [v for v in views if v and v.get("traces")]
        if not views:
            return "<p class='text-muted'>No data.</p>"
        self.n_plots += 1
        plot_id = f"vc_plot_{self.n_plots}"
        spec = {"id": plot_id, "views": views}
        return (
            f"<div class='vc-plot-buttons' id='{plot_id}-buttons'></div>"
            f"<div id='{plot_id}' style='height:{height}px'></div>"
            + (f"<p class='vc-note'>{caption}</p>" if caption else "")
            + "<script>(window.vcPlotSpecs = window.vcPlotSpecs || []).push("
            + json.dumps(spec, separators=(",", ":"), default=float).replace("</", "<\\/")
            + ");</script>"
        )


def cards(items, accent=None):
    """A row of cards: items = [(value, label[, subtext[, accent colour]])]."""
    out = []
    for item in items:
        value, label = item[0], item[1]
        sub = f"<div class='s'>{esc(item[2])}</div>" if len(item) > 2 and item[2] else ""
        # colour of the left border: of the card, else of the row
        color = item[3] if len(item) > 3 and item[3] else accent
        style = f" style='--vc-accent:{color}'" if color else ""
        out.append(
            f"<div class='vc-card'{style}><div class='v'>{esc(value)}</div>"
            f"<div class='l'>{esc(label)}</div>{sub}</div>"
        )
    return "<div class='vc-cards'>" + "".join(out) + "</div>"


def toggle_buttons(group, keys_labels):
    """Buttons that show the elements of one key of a toggle group (window.vcToggle)."""
    buttons = "".join(
        f"<button type='button' class='btn btn-sm btn-outline-secondary{' active' if i == 0 else ''}' "
        f"onclick=\"vcToggle('{group}', '{key}', this)\">{esc(label)}</button>"
        for i, (key, label) in enumerate(keys_labels)
    )
    return f"<div class='vc-plot-buttons'>{buttons}</div>"


def bar_layout(ytitle, ymax=None, barmode="stack", xtitle=None, extra=None):
    """Layout of a vertical bar plot with the samples on the x-axis."""
    layout = {
        "barmode": barmode,
        "yaxis": {"title": {"text": ytitle}, "rangemode": "tozero"},
        # no vertical grid lines between the categories (samples)
        "xaxis": {
            "type": "category",
            "showgrid": False,
            **({"title": {"text": xtitle}} if xtitle else {}),
        },
        "hovermode": "closest",
    }
    if ymax:
        layout["yaxis"]["range"] = [0, ymax]
    layout.update(extra or {})
    return layout


# ----------------------------------------------------------------------------- sections


def section_summary(mqc, d):
    """Summary: consistency banner, recall/precision/F1 cards per matching method,
    context cards, key findings, links."""
    parts = [CSS, JS_LIB]
    checks = d["checks"]
    if checks.empty:
        parts.append("<div class='vc-banner fail'>No consistency checks available.</div>")
    else:
        failed = checks[checks["status"] != "ok"]
        if failed.empty:
            parts.append(
                f"<div class='vc-banner ok'>All {len(checks)} consistency checks passed: the "
                "counts of all steps add up.</div>"
            )
        else:
            parts.append(
                f"<div class='vc-banner fail'>{len(failed)} of {len(checks)} consistency "
                "checks FAILED: the results may be wrong (section 6.1 Consistency checks).</div>"
            )
    agg = d["metrics"]
    pairs = comparisons(d["per_pair"])
    # recall, precision, F1 per comparison (summed over samples), one card set per method
    if pairs and not agg.empty:
        parts.append(toggle_buttons("summary_method", [(m[0], m[3]) for m in METHODS]))
        for key, stype, stratum, label, text in METHODS:
            rows = overall_rows(agg, stype, stratum)
            block = [
                f"<div data-vc-group='summary_method' data-vc-key='{key}' "
                f"class='{'' if key == 'vcfeval' else 'vc-hidden'}'>"
                f"<p class='vc-note'><b>{esc(label)}</b>: {esc(text)} Counts summed over the "
                "samples (micro average).</p>"
            ]
            for truth, query in pairs:
                row = rows[(rows["truth"] == truth) & (rows["query"] == query)]
                if row.empty:
                    continue
                r = row.iloc[0]
                block.append(
                    f"<div class='vc-group'>{esc(query)} compared with {esc(truth)} "
                    f"(truth), {fmt_int(r['n_samples'])} samples</div>"
                )
                block.append(
                    cards(
                        [
                            (
                                fmt_num(r["micro_recall"]),
                                "Recall",
                                f"{fmt_int(r['tp_baseline'])} of "
                                f"{fmt_int(r['tp_baseline'] + r['fn'])} {truth} variants found",
                            ),
                            (
                                fmt_num(r["micro_precision"]),
                                "Precision",
                                f"{fmt_int(r['tp_call'])} of {fmt_int(r['tp_call'] + r['fp'])} "
                                f"{query} calls confirmed",
                            ),
                            (fmt_num(r["micro_f1"]), "F1", "harmonic mean"),
                            (fmt_int(r["fn"]), f"in {truth} only (FN)", f"missed by {query}"),
                            (fmt_int(r["fp"]), f"in {query} only (FP)", f"not in {truth}"),
                        ]
                    )
                )
            block.append("</div>")
            parts.append("".join(block))
    # context of the run
    samples = sorted(d["per_pair"]["sample"].unique()) if not d["per_pair"].empty else []
    context = [(len(samples), "samples (tumor/normal pairs)")]
    if pairs:
        context.append((pairs[0][0], "truth technology"))
        context.append((", ".join(sorted({q for _, q in pairs})), "query technologies"))
    steps = d["steps"]
    if not steps.empty:
        context.append(
            (
                fmt_bp(steps["bp"].iloc[-1]),
                "compared region",
                f"{fmt_pct(steps['bp'].iloc[-1], steps['bp'].iloc[0], 2)} of the genome",
            )
        )
    if not d["truth"].empty:
        context.append((fmt_int(len(d["truth"])), "truth variants compared", "all samples"))
    req = d["required_summary"]
    if not req.empty:
        total = req[req["sample"] == "all"]
        context.append(
            (
                f"{int(total['passed'].sum())}/{int(total['required'].sum())}",
                "required variants found",
                "summed over the technologies",
            )
        )
    parts.append("<div class='vc-group'>This run</div>" + cards(context, accent="#9c755f"))
    findings = key_findings(d)
    if findings:
        parts.append(
            "<div class='vc-group'>Key findings</div><ul>"
            + "".join(f"<li>{f}</li>" for f in findings)
            + "</ul>"
        )
    links = [
        ("truth-variants", "Truth variants"),
        ("discordant-variants", "Discordant variants"),
        ("compared-variants", "All compared variants"),
        ("genes", "Genes"),
    ]
    if not d["required_variants"].empty:
        links.append(("required-variants", "Required variants"))
    parts.append(
        "<div class='vc-links'>Variant tables (datavzrd): "
        + "".join(
            f"<a class='btn btn-sm btn-outline-primary' href='{esc(mqc.link(v))}' "
            f"target='_blank'>{esc(t)}</a>"
            for v, t in links
        )
        + "</div>"
    )
    mqc.html(
        "summary",
        "Results at a glance",
        "**What is shown:** the main numbers of the run on one screen.\n\n"
        "- **Banner:** do the counts of all steps add up (section 6.1)?\n"
        "- **Cards per query technology:** recall, precision and F1, and the numbers of "
        "variants in one technology only (FN, FP), summed over all samples (micro average). "
        "The buttons switch the matching method.\n"
        "- **This run:** samples, technologies, the size of the compared region and the "
        "number of compared truth variants.\n"
        "- **Key findings:** the range of recall over the samples and the most common reason "
        "of discordant variants.\n"
        "- **Variant tables:** links to the table browser (datavzrd) with every variant.",
        "\n".join(parts),
    )


def key_findings(d):
    """Short sentences about the most important results (HTML)."""
    findings = []
    pp = overall_rows(d["per_pair"])
    for (_truth, query), rows in pp.groupby(["truth", "query"]) if not pp.empty else []:
        rows = rows.dropna(subset=["recall"])
        if rows.empty:
            continue
        low = rows.loc[rows["recall"].idxmin()]
        high = rows.loc[rows["recall"].idxmax()]
        findings.append(
            f"<b>{esc(query)}</b> recall per sample ranges from {fmt_num(low['recall'])} "
            f"({esc(low['sample'])}) to {fmt_num(high['recall'])} ({esc(high['sample'])})."
        )
    variants = d["variants"]
    if not variants.empty and "class" in variants.columns:
        # with several query technologies, every sentence names its comparison
        several = variants["query"].nunique() > 1
        for (truth, query), rows in variants.groupby(["truth", "query"]):
            prefix = f"<b>{esc(query)}</b> vs {esc(truth)}: " if several else ""
            for status, where in (("FN", f"in {truth} only"), ("FP", f"in {query} only")):
                sub = rows[rows["status"] == status]
                if sub.empty:
                    continue
                top = sub["class"].value_counts()
                cls = top.index[0]
                findings.append(
                    f"{prefix}{fmt_int(len(sub))} variants {esc(where)} ({status}); most common "
                    f"reason: <b>{esc(CLASS_LABEL.get(cls, cls))}</b> "
                    f"({fmt_int(top.iloc[0])}, {fmt_pct(top.iloc[0], len(sub), 0)})."
                )
    req = d["required_variants"]
    if not req.empty:
        failed = req[req["status"] != "PASS"]
        if failed.empty:
            findings.append("Every required variant was found by every technology.")
        else:
            names = ", ".join(
                sorted(
                    {
                        f"{value_of(r, 'gene') or value_of(r, 'label') or r['pos']} "
                        f"({r['sample']}, {r['technology']})"
                        for _, r in failed.iterrows()
                    }
                )[:5]
            )
            findings.append(
                f"<b>{len(failed)} required variant check(s) failed</b>: {esc(names)}"
                + (" ..." if len(failed) > 5 else "")
            )
    return findings


def section_introduction(mqc, d, args):
    """Introduction: this run and the comparison methods as lists, technologies table."""
    params = d["params"]
    pairs = comparisons(d["per_pair"])
    truth = pairs[0][0] if pairs else param(params, "truth_technology")
    queries = sorted({q for _, q in pairs})
    samples = sorted(d["per_pair"]["sample"].unique()) if not d["per_pair"].empty else []
    sheet = d["samplesheet"]
    run = [
        f"<b>Question</b>: are the somatic SNVs and indels of the same tumor/normal pairs "
        f"found with every sequencing technology? <b>{esc(truth)}</b> is the truth for all "
        f"samples; compared with it: <b>{esc(', '.join(queries) or 'NA')}</b>.",
        f"<b>Samples</b>: {len(samples)} ({esc(', '.join(samples[:12]))}"
        f"{' ...' if len(samples) > 12 else ''}).",
    ]
    if not sheet.empty and "caller" in sheet.columns:
        callers = ", ".join(sorted(sheet["caller"].astype(str).unique()))
        run.append(f"<b>Callers</b>: {esc(callers)}; input are the unfiltered VCFs.")
    if not d["steps"].empty:
        steps = d["steps"]
        run.append(
            f"<b>Compared region</b>: {fmt_bp(steps['bp'].iloc[-1])} "
            f"({fmt_pct(steps['bp'].iloc[-1], steps['bp'].iloc[0], 2)} of the genome): the "
            "intersection of all target BEDs"
            + (
                ", restricted per sample to the bases with enough depth in every alignment"
                if args.callable_regions == "true"
                else ""
            )
            + "."
        )
    run.append(
        "<b>Filters</b>: "
        + (
            "only PASS calls are compared"
            if args.filter_pass_only != "false"
            else "all calls are compared"
        )
        + f"; uniform thresholds for every technology: VAF &ge; {esc(args.min_vaf)}, "
        f"depth &ge; {esc(args.min_depth)}, alt reads &ge; {esc(args.min_alt_reads)}."
    )
    methods = [
        "<b>Same treatment of every VCF</b>: sample selection, contig renaming, split of "
        "multiallelic records, left-alignment, atomization of MNVs, deduplication "
        "(bcftools norm); VAF, depth and alt reads computed the same way for every caller.",
        "<b>Matching</b>: RTG vcfeval (haplotype-aware, main result); bcftools isec with the "
        "same allele and with the same position as a check of representation differences. "
        "The buttons of the Summary and of the recall/precision plot switch between them.",
        "<b>Discordant variants</b> (FN: in the truth only; FP: in the query only) are explained "
        "from the view of the technology that misses them: exactly one main reason (non-PASS "
        "call, removed by a filter, or not called with low coverage, with alt reads, without "
        "alt reads), then further reasons: germline evidence in the normal and the context "
        "(indel, low VAF, repeat, GIAB and other stratifications). Read evidence is counted "
        "for the called alleles, indels by their haplotype sequence.",
        "<b>Annotation</b>: Ensembl VEP (gene, consequence, impact, HGVS), repeat regions "
        "(msisensor-pro) and GIAB stratifications; enrichment of features among discordant "
        "variants with Fisher's exact test (Benjamini-Hochberg q).",
        "<b>Checks</b>: every input is checked against the reference first; at the end the "
        "counts of all steps must add up (consistency checks).",
    ]
    parts = [
        "<div class='vc-group'>This run</div><ul>",
        *[f"<li>{r}</li>" for r in run],
        "</ul><div class='vc-group'>How the technologies are compared</div><ul>",
        *[f"<li>{m}</li>" for m in methods],
        "</ul>",
    ]
    # one row per technology: library type, role, caller, samples, target BED and size
    if not sheet.empty:
        targets = d["targets"]
        rows = []
        for tech, sub in sheet.groupby("technology", sort=False):
            bed = (
                sub["target_bed"].astype(str).str.split("/").str[-1].iloc[0]
                if "target_bed" in sub
                else ""
            )
            size = ""
            if not targets.empty:
                hit = targets[targets["region"] == tech]
                if not hit.empty:
                    size = fmt_bp(hit["bp"].iloc[0])
            rows.append(
                {
                    "technology": tech,
                    "role": "truth" if tech == truth else "query",
                    "library type": value_of(sub.iloc[0], "library_type"),
                    "caller": (
                        ", ".join(sorted(sub["caller"].astype(str).unique()))
                        if "caller" in sub
                        else ""
                    ),
                    "samples": len(sub),
                    "target BED": bed if bed not in ("", "nan", "n/a", "NA") else "(genome)",
                    "target size": size,
                }
            )
        parts.append(
            html_table(pd.DataFrame(rows), link=mqc.link("samplesheet"), link_text="Samplesheet")
        )
    parts.append(
        f"<p class='vc-note'>Pipeline version {esc(args.pipeline_version)}, run "
        f"{esc(args.run_name)}, started {esc(args.start)}; report content written "
        f"{datetime.now().isoformat(timespec='seconds')}.</p>"
    )
    mqc.html(
        "introduction",
        "What was compared and how",
        "**What is shown:** the question of the run, its inputs and settings, and how the "
        "pipeline compares the technologies. The table lists every technology with its "
        "library type, role (truth or query), caller and target BED.",
        "\n".join(parts),
    )


def section_recall_precision(mqc, d):
    """Recall vs. precision per sample: native MultiQC scatter plot, one dataset per
    matching method; and the metrics table per sample."""
    pp = d["per_pair"]
    if pp.empty:
        return
    queries = sorted(pp["query"].unique())
    colors = {q: PALETTE[i % len(PALETTE)] for i, q in enumerate(queries)}
    datasets, labels = [], []
    for _key, stype, stratum, label, _ in METHODS:
        rows = overall_rows(pp, stype, stratum).dropna(subset=["recall", "precision"])
        if rows.empty:
            continue
        datasets.append(
            {
                f"{r['sample']} ({r['query']})": {
                    "x": round(float(r["recall"]), 4),
                    "y": round(float(r["precision"]), 4),
                    "color": colors[r["query"]],
                }
                for _, r in rows.iterrows()
            }
        )
        labels.append({"name": label, "xlab": "Recall", "ylab": "Precision"})
    if datasets:
        legend = ", ".join(
            f"{q}: <span style='color:{c}'>&#9679;</span>" for q, c in colors.items()
        )
        mqc.json(
            "recall_precision",
            {
                "section_name": "Recall and precision per sample",
                "description": (
                    "**What is shown:** one dot per sample and query technology (hover for "
                    f"the sample name); colours: {legend}.\n\n"
                    "- **x-axis: recall** = TP / (TP + FN), the share of the truth variants "
                    "found; **y-axis: precision** = TP / (TP + FP), the share of the query "
                    "calls confirmed. Both axes span 0 to 1, so runs can be compared.\n"
                    "- **Top right is best.** A dot far left: the query misses truth variants "
                    "(chapter 3, FN); a dot low: the query has calls the truth lacks (FP).\n"
                    "- **Buttons:** the matching method: RTG vcfeval (haplotype-aware, main "
                    "result), bcftools isec with the same allele or the same position."
                ),
                "plot_type": "scatter",
                "pconfig": {
                    "id": "vc_recall_precision_plot",
                    "title": "Recall and precision per sample",
                    "xlab": "Recall",
                    "ylab": "Precision",
                    "xmin": 0,
                    "xmax": 1,
                    "ymin": 0,
                    "ymax": 1,
                    "square": True,
                    "data_labels": labels,
                },
                "data": datasets,
            },
        )
    # metrics table per sample and query technology (all three methods)
    data = {}
    for (sample, _truth, query), _rows in overall_rows(pp).groupby(["sample", "truth", "query"]):
        row = {}
        for key, stype, stratum, _, _ in METHODS:
            sub = overall_rows(pp, stype, stratum)
            sub = sub[(sub["sample"] == sample) & (sub["query"] == query)]
            if sub.empty:
                continue
            r = sub.iloc[0]
            if key == "vcfeval":
                for c in ("tp_baseline", "fn", "fp"):
                    row[c] = None if pd.isna(r[c]) else int(r[c])
            for c in ("recall", "precision", "f1"):
                row[f"{c}_{key}"] = None if pd.isna(r[c]) else round(float(r[c]), 4)
        data[f"{sample} ({query})"] = row
    headers = {
        "tp_baseline": {
            "title": "TP",
            "description": "truth variants found (vcfeval)",
            "format": "{:,.0f}",
            "scale": "Greens",
        },
        "fn": {
            "title": "FN",
            "description": "truth variants missed (in the truth only)",
            "format": "{:,.0f}",
            "scale": "Oranges",
        },
        "fp": {
            "title": "FP",
            "description": "query calls without truth variant (in the query only)",
            "format": "{:,.0f}",
            "scale": "Purples",
        },
    }
    for key, _, _, label, _ in METHODS:
        for c in ("recall", "precision", "f1"):
            headers[f"{c}_{key}"] = {
                "title": f"{c.capitalize()} ({key})",
                "description": f"{c} with {label}",
                "min": 0,
                "max": 1,
                "format": "{:.3f}",
                "scale": "RdYlGn",
                "hidden": key != "vcfeval" and c == "f1",
            }
    if data:
        mqc.json(
            "metrics_table",
            {
                "section_name": "Metrics table per sample",
                "description": (
                    "**What is shown:** the numbers behind 1.3 and 1.4, one row per sample "
                    "and query technology.\n\n"
                    "- **TP, FN, FP** (vcfeval): truth variants found, truth variants missed, "
                    "query calls without a truth variant.\n"
                    "- **Recall, precision and F1** with the three matching methods (F1 of the "
                    "isec methods hidden; *Configure columns* shows it).\n"
                    "- Rows can be sorted, hidden and highlighted with the toolbox."
                ),
                "plot_type": "table",
                "pconfig": {
                    "id": "vc_metrics_table",
                    "col1_header": "Sample (query)",
                    "no_violin": True,
                },
                "headers": headers,
                "data": data,
            },
        )


def section_truth_found(mqc, d):
    """Truth variants found per sample (vertical bars, label found/total), per method."""
    pp = d["per_pair"]
    if pp.empty:
        return
    views = []
    for _key, stype, stratum, label, _ in METHODS:
        rows = overall_rows(pp, stype, stratum)
        if rows.empty:
            continue
        traces, top = [], 0
        for i, (query, sub) in enumerate(rows.groupby("query")):
            sub = sub.sort_values("sample")
            total = sub["tp_baseline"] + sub["fn"]
            truth = sub["truth"].iloc[0]
            top = max(top, float(total.max()))
            traces.append(
                {
                    "type": "bar",
                    "name": f"found by {query}",
                    "x": sub["sample"].tolist(),
                    "y": sub["tp_baseline"].tolist(),
                    "s": sub["sample"].tolist(),
                    "renameX": True,
                    "text": [
                        f"{int(f)}/{int(t)}" for f, t in zip(sub["tp_baseline"], total, strict=True)
                    ],
                    "textposition": "outside",
                    "cliponaxis": False,
                    "marker": {"color": PALETTE[i % len(PALETTE)]},
                    "hovertext": [
                        f"{s}: {query} found {int(f)} of {int(t)} {truth} variants ({fmt_pct(f, t)})"
                        for s, f, t in zip(sub["sample"], sub["tp_baseline"], total, strict=True)
                    ],
                    "hoverinfo": "text",
                }
            )
        views.append(
            {
                "label": label,
                "traces": traces,
                "layout": bar_layout(
                    "Truth variants found", ymax=top * 1.15 or None, barmode="group"
                ),
            }
        )
    body = [
        mqc.plot(
            views,
            caption="Bars: truth variants also called by the query technology; label: found / "
            "all compared truth variants of the sample. Buttons: matching method.",
        )
    ]
    mqc.html(
        "truth_found",
        "Truth variants found",
        "**What is shown:** how many of the truth variants of every sample the query "
        "technology found.\n\n"
        "- **Bar:** truth variants also called by the query; **label:** found / all compared "
        "truth variants of the sample.\n"
        "- Compare the bar heights between samples to see samples with few variants, where "
        "one missed variant changes the recall a lot.",
        "\n".join(body),
    )


def section_overlap(mqc, d):
    """Venn diagrams of truth and query per sample (TP in both, FN truth only, FP query
    only), one set of diagrams per matching method (toggle)."""
    pp = d["per_pair"]
    if pp.empty:
        return
    parts = [toggle_buttons("overlap_method", [(m[0], m[3]) for m in METHODS])]
    for key, stype, stratum, label, text in METHODS:
        rows = overall_rows(pp, stype, stratum).sort_values(["query", "sample"])
        block = [
            f"<div data-vc-group='overlap_method' data-vc-key='{key}' "
            f"class='{'' if key == 'vcfeval' else 'vc-hidden'}'>"
            f"<p class='vc-note'><b>{esc(label)}</b>: {esc(text)}</p><div class='vc-venns'>"
        ]
        for _, r in rows.iterrows():
            tpb, fn, tpc, fp = (int(r[c]) for c in ("tp_baseline", "fn", "tp_call", "fp"))
            block.append(
                f"<div class='vc-venn' data-sample='{esc(r['sample'])}'>"
                f"<div class='vc-group'>{esc(r['sample'])}: {esc(r['truth'])} vs "
                f"{esc(r['query'])}</div>"
                + venn2_svg(r["truth"], r["query"], fn, tpb, fp, tpc)
                + "</div>"
            )
        block.append("</div></div>")
        parts.append("".join(block))
    parts.append(
        "<p class='vc-note'>Left: variants only in the truth technology (missed by the "
        "query, FN); middle: truth variants called by both technologies (TP); right: calls "
        "only in the query technology (FP). With RTG vcfeval one truth record can match "
        "several query records (e.g. a complex variant against two simple ones): the number "
        "of query calls is then shown in brackets. Hover a number for the details.</p>"
    )
    mqc.html(
        "overlap",
        "Overlap of the technologies per sample",
        "**What is shown:** one Venn diagram per sample and query technology.\n\n"
        "- **Middle: truth variants** called by both technologies (TP).\n"
        "- **Left:** in the truth technology only (FN, missed by the query); **right:** in the "
        "query technology only (FP).\n"
        "- **Buttons:** the matching method; the overlap grows when variants that are "
        "represented differently are matched by their haplotype (vcfeval).",
        "\n".join(parts),
    )


def section_known(mqc, d):
    """Known variants: required variants (PASS/FAIL matrix) and the benchmark against an
    external truth set."""
    req, req_var, bench = d["required_summary"], d["required_variants"], d["known_benchmark"]
    if req.empty and bench.empty:
        return
    parts = []
    if not req.empty:
        total = req[req["sample"] == "all"]
        items = [
            (
                f"{int(r['passed'])}/{int(r['required'])}",
                f"{r['technology']}: required variants found",
                f"list {r['list']}" + ("" if r["all_passed"] == "yes" else " - not all found"),
                # green: all found, red: at least one missing
                "#59a14f" if r["all_passed"] == "yes" else "#e15759",
            )
            for _, r in total.iterrows()
        ]
        parts.append("<div class='vc-group'>Required variants (clinical checklist)</div>")
        parts.append(cards(items))
        parts.append(
            "<p class='vc-note'>Every required variant must be found by every technology in "
            "every sample: PASS = called with FILTER PASS inside the comparison region; FAIL "
            "states why not (filtered, outside the region, not called). Hover a cell for the "
            "details.</p>"
        )
        if not req_var.empty:
            parts.append(required_matrix(req_var))
            parts.append(
                f"<p class='vc-note'><a href='{esc(mqc.link('required-variants'))}' "
                "target='_blank'>Required variants in the table browser &rarr;</a></p>"
            )
    if not bench.empty:
        agg = overall_rows(bench)
        parts.append("<div class='vc-group'>Benchmark against an external truth set</div>")
        if not agg.empty and "micro_recall" in agg.columns:
            for _, r in agg.iterrows():
                parts.append(
                    f"<div class='vc-note'>{esc(r['query'])} vs {esc(r['truth'])}</div>"
                    + cards(
                        [
                            (
                                fmt_num(r["micro_recall"]),
                                "Recall",
                                f"{fmt_int(r['tp_baseline'])} of {fmt_int(r['tp_baseline'] + r['fn'])} truth-set variants",
                            ),
                            (
                                fmt_num(r["micro_precision"]),
                                "Precision",
                                f"{fmt_int(r['fp'])} calls not in the truth set",
                            ),
                            (fmt_num(r["micro_f1"]), "F1"),
                        ]
                    )
                )
        per = overall_rows(d["known_benchmark_pairs"])
        if not per.empty:
            views = []
            for metric in ("recall", "precision", "f1"):
                traces = []
                for i, (query, sub) in enumerate(per.groupby("query")):
                    sub = sub.sort_values("sample")
                    traces.append(
                        {
                            "type": "bar",
                            "name": query,
                            "x": sub["sample"].tolist(),
                            "y": sub[metric].tolist(),
                            "s": sub["sample"].tolist(),
                            "renameX": True,
                            "marker": {"color": PALETTE[i % len(PALETTE)]},
                            "text": [fmt_num(v) for v in sub[metric]],
                            "textposition": "outside",
                            "cliponaxis": False,
                        }
                    )
                views.append(
                    {
                        "label": metric.capitalize(),
                        "traces": traces,
                        "layout": bar_layout(metric.capitalize(), ymax=1.1, barmode="group"),
                    }
                )
            parts.append(
                mqc.plot(
                    views,
                    height=360,
                    caption="Per sample, within the confident regions of the truth set.",
                )
            )
    mqc.html(
        "known",
        "Required variants and external truth sets",
        "**What is shown:** the validation with variants known beforehand.\n\n"
        "- **Required variants:** cards per technology (found / required; green: all found), "
        "and a matrix variant × technology × sample: **PASS** = called with FILTER PASS inside "
        "the comparison region; **FAIL** states why not (filtered, outside the region, not "
        "called). Hover a cell for FILTER, VAF and depth.\n"
        "- **Benchmark:** recall, precision and F1 of every technology against the external "
        "truth set, summed over the samples (cards) and per sample (bars), within the "
        "confident regions of the set and the comparison region.",
        "\n".join(parts),
    )


def required_matrix(req_var):
    """HTML matrix of the required variants: rows variants, columns technology x sample,
    cells PASS / FAIL with the details as tooltip."""
    df = req_var.copy()
    df["key"] = (
        df["chrom"].astype(str) + ":" + df["pos"].astype(str) + " " + df["ref"] + ">" + df["alt"]
    )
    techs = sorted(df["technology"].unique())
    samples = sorted(df["sample"].unique())
    head1 = "<th rowspan='2'>variant</th><th rowspan='2'>position</th>" + "".join(
        f"<th colspan='{len(samples)}'>{esc(t)}</th>" for t in techs
    )
    head2 = "".join(f"<th>{esc(s)}</th>" for _ in techs for s in samples)
    rows = []
    for (label, gene, key), sub in df.groupby(["label", "gene", "key"], sort=False):
        cells = []
        for tech in techs:
            for sample in samples:
                hit = sub[(sub["technology"] == tech) & (sub["sample"] == sample)]
                if hit.empty:
                    cells.append("<td class='text-muted'>-</td>")
                    continue
                r = hit.iloc[0]
                ok = r["status"] == "PASS"
                tip = (
                    f"{r['status']}; FILTER {value_of(r, 'filter') or '-'}; VAF "
                    f"{value_of(r, 'vaf') or '-'}; depth {value_of(r, 'depth') or '-'}; "
                    f"in region {value_of(r, 'in_region') or '-'}"
                )
                text = "PASS" if ok else str(r["status"]).replace("FAIL_", "FAIL: ")
                cells.append(
                    f"<td class='{'pass' if ok else 'fail'}' title='{esc(tip)}'>{esc(text)}</td>"
                )
        rows.append(
            # gene in bold (left out if the list gives none), then the label of the list
            f"<tr><td><b>{esc(gene if str(gene) not in ('.', 'nan', '') else '')}</b> {esc(label)}</td>"
            f"<td>{esc(key)}</td>{''.join(cells)}</tr>"
        )
    return (
        "<div class='table-responsive'><table class='table table-sm table-bordered vc-matrix'>"
        f"<thead><tr>{head1}</tr><tr>{head2}</tr></thead><tbody>{''.join(rows)}</tbody>"
        "</table></div>"
    )


def section_why(mqc, d):
    """Why variants differ: discordant variants per sample (vertical bars) stacked by
    class; views FN + FP, FN only, FP only (per query technology)."""
    variants = d["variants"]
    if variants.empty or "class" not in variants.columns:
        return
    disc = variants[variants["status"].isin(["FN", "FP"])]
    samples = sorted(variants["sample"].unique())
    views = []
    for (truth, query), _rows in variants.groupby(["truth", "query"]):
        sub_disc = disc[(disc["truth"] == truth) & (disc["query"] == query)]
        prefix = f"{query}: " if variants["query"].nunique() > 1 else ""
        for statuses, label in (
            (["FN", "FP"], f"{prefix}all discordant"),
            (["FN"], f"{prefix}in {truth} only (FN)"),
            (["FP"], f"{prefix}in {query} only (FP)"),
        ):
            part = sub_disc[sub_disc["status"].isin(statuses)]
            counts = part.groupby(["sample", "class"]).size().unstack(fill_value=0)
            counts = counts.reindex(samples, fill_value=0)
            totals = counts.sum(axis=1) if not counts.empty else pd.Series(0, index=samples)
            traces = []
            for cls in CLASS_ORDER:
                if cls not in counts.columns or counts[cls].sum() == 0:
                    continue
                traces.append(
                    {
                        "type": "bar",
                        "name": CLASS_LABEL[cls],
                        "x": samples,
                        "y": counts[cls].tolist(),
                        "s": samples,
                        "renameX": True,
                        "marker": {"color": CLASS_COLOR[cls]},
                        "hovertemplate": "%{x}: %{y} variants<br>" + cls + "<extra></extra>",
                    }
                )
            # total above every bar
            if traces:
                traces.append(
                    {
                        "type": "scatter",
                        "mode": "text",
                        "x": samples,
                        "y": totals.tolist(),
                        "s": samples,
                        "renameX": True,
                        "text": [str(int(t)) for t in totals],
                        "textposition": "top center",
                        "showlegend": False,
                        "hoverinfo": "skip",
                        "cliponaxis": False,
                    }
                )
            views.append(
                {
                    "label": label,
                    "traces": traces,
                    "layout": bar_layout(
                        "Discordant variants", ymax=float(totals.max() or 1) * 1.15
                    ),
                }
            )
    # further reasons: per sample, how many discordant variants have each (not stacked:
    # a variant can have several)
    further_views = []
    for (truth, query), _rows in variants.groupby(["truth", "query"]):
        sub_disc = disc[(disc["truth"] == truth) & (disc["query"] == query)]
        prefix = f"{query}: " if variants["query"].nunique() > 1 else ""
        reason_lists = (
            sub_disc["reasons"].astype(str).str.split(",") if "reasons" in sub_disc else []
        )
        names = further_reasons(sub_disc)
        traces = []
        for i, name in enumerate(names):
            has = [name in r for r in reason_lists]
            per_sample = (
                sub_disc.assign(has=has)
                .groupby("sample")["has"]
                .sum()
                .reindex(samples, fill_value=0)
            )
            traces.append(
                {
                    "type": "bar",
                    "name": CLASS_LABEL.get(name, name),
                    "x": samples,
                    "y": [int(v) for v in per_sample],
                    "s": samples,
                    "renameX": True,
                    "marker": {"color": CLASS_COLOR.get(name, PALETTE[i % len(PALETTE)])},
                    "hovertemplate": "%{x}: %{y} discordant variants<br>"
                    + name
                    + "<extra></extra>",
                }
            )
        further_views.append(
            {
                "label": f"{prefix}FN + FP",
                "traces": traces,
                "layout": bar_layout("Discordant variants with the reason", barmode="group"),
            }
        )
    main_items = "".join(
        f"<li>{badge(CLASS_LABEL[c], CLASS_COLOR[c])} <code>{c}</code>: {esc(CLASS_TEXT[c])}</li>"
        for c in CLASS_ORDER
    )
    further_items = (
        "".join(
            f"<li>{badge(CLASS_LABEL[c], CLASS_COLOR[c])} <code>{c}</code>: {esc(FURTHER_TEXT[c])}</li>"
            for c in FURTHER_ORDER
        )
        + f"<li><code>giab_&lt;name&gt;</code>, stratifications: {esc(STRATUM_TEXT)}</li>"
    )
    body = [
        "<p class='vc-note'>Every discordant variant has exactly one <b>main reason</b>: the "
        "call and the reads of the technology that misses it (called but filtered, or not "
        "called with low coverage, with alt reads, without alt reads, or no alignment). "
        "<b>Further reasons</b> can come in addition: germline evidence in the normal and the "
        "context (indel, low VAF, repeat, GIAB and other stratifications). The column "
        "<code>reasons</code> of the tables lists all of them.</p>",
        "<div class='vc-group'>Main reason</div>",
        mqc.plot(
            views,
            caption="Discordant variants per sample by their main reason (exactly one per "
            "variant, so the bars add up). FN: in the truth technology only (missed by the "
            "query); FP: in the query technology only. Number above a bar: all discordant "
            "variants of the sample.",
        ),
        "<div class='vc-group'>Further reasons</div>",
        mqc.plot(
            further_views,
            height=380,
            caption="Number of discordant variants (FN + FP) per sample with each further "
            "reason. A variant can have several further reasons, so the bars don't add up.",
        ),
        f"<details><summary>What the reasons mean</summary><b>Main reasons</b><ul>{main_items}</ul>"
        f"<b>Further reasons</b><ul>{further_items}</ul></details>",
    ]
    mqc.html(
        "why",
        "Why variants differ",
        "**What is shown:** the discordant variants of every sample by the reasons why one "
        "technology misses them.\n\n"
        "- **Main reason** (stacked bars, they add up to all discordant variants): the call "
        "and the reads of the technology that misses the variant. Buttons: FN + FP, FN only "
        "(in the truth only), FP only (in the query only).\n"
        "- **Further reasons** (grouped bars, they don't add up): germline evidence in the "
        "normal and the context (indel, low VAF, repeat, GIAB and other stratifications).\n"
        "- **How to read it:** *filtered* reasons point to the caller settings, *no or "
        "subthreshold evidence* to the data (depth, library), *low coverage* and *absent* to "
        "the target design.",
        "\n".join(body),
    )


def further_reasons(disc):
    """The further reasons that occur among discordant variants, in a fixed order:
    FURTHER_ORDER first, then the stratifications (alphabetically)."""
    if disc.empty or "reasons" not in disc:
        return []
    seen = {r for rs in disc["reasons"].astype(str) for r in rs.split(",")}
    seen -= set(CLASS_ORDER) | {"concordant"}
    return [r for r in FURTHER_ORDER if r in seen] + sorted(seen - set(FURTHER_ORDER))


def vaf_bin(values):
    """VAF values as the bins of VAF_LABELS ('NA' for missing)."""
    return pd.cut(values, VAF_BINS, labels=VAF_LABELS, right=False).astype(str).replace("nan", "NA")


def section_characterization(mqc, d, fdr):
    """The discordant variants in detail, per comparison: what characterizes the variants
    found in one technology only (lists, profile plot, VAF, reason combinations, top
    protein-altering variants, genes, enrichment)."""
    variants = d["variants"]
    if variants.empty or "class" not in variants.columns:
        return
    variants = add_variant_columns(variants)
    parts = []
    for (truth, query), rows in variants.groupby(["truth", "query"]):
        tp, fn, fp = (rows[rows["status"] == s] for s in ("TP", "FN", "FP"))
        sets = {
            "TP": ("both (TP)", tp, "truth_vaf"),
            "FN": (f"{truth} only (FN)", fn, "truth_vaf"),
            "FP": (f"{query} only (FP)", fp, "query_vaf"),
        }
        parts.append(f"<h4>{esc(query)} compared with {esc(truth)}</h4>")
        # summary list
        bullets = [
            f"{fmt_int(len(tp))} variants found by both, {fmt_int(len(fn))} in "
            f"<b>{esc(truth)} only</b> (missed by {esc(query)}, FN) and {fmt_int(len(fp))} in "
            f"<b>{esc(query)} only</b> (FP)."
        ]
        for status in ("FN", "FP"):
            name, sub, vaf_col = sets[status]
            if sub.empty:
                continue
            top = sub["class"].value_counts().head(3)
            reasons_text = ", ".join(f"{CLASS_LABEL.get(c, c)} ({n})" for c, n in top.items())
            text = f"<b>{esc(name)}</b>: main reasons {esc(reasons_text)}"
            if "reasons" in sub.columns:
                multi = (sub["reasons"].astype(str).str.count(",") > 0).sum()
                text += f"; {fmt_int(multi)} of {fmt_int(len(sub))} have more than one reason"
            indel_share = (sub["type"] != "SNV").mean() if "type" in sub else None
            tp_indel = (tp["type"] != "SNV").mean() if "type" in tp and len(tp) else None
            if indel_share is not None:
                text += (
                    f"; {fmt_pct(indel_share, 1, 0)} indels (TP: "
                    f"{fmt_pct(tp_indel, 1, 0) if tp_indel is not None else 'NA'})"
                )
            vaf = pd.to_numeric(sub[vaf_col], errors="coerce").median() if vaf_col in sub else None
            tp_vaf = pd.to_numeric(tp["truth_vaf"], errors="coerce").median() if len(tp) else None
            if vaf is not None and not pd.isna(vaf):
                text += f"; median VAF {fmt_num(vaf, 2)} (TP: {fmt_num(tp_vaf, 2)})"
            if "consequence_group" in sub.columns:
                prot = (sub["consequence_group"] == "non_synonymous").sum()
                text += f"; {fmt_int(prot)} protein-altering"
            if "gene" in sub.columns:
                genes = sub[sub["gene"].astype(str) != "."]["gene"].value_counts().head(3)
                if not genes.empty:
                    text += "; most affected genes " + ", ".join(
                        f"{esc(g)} ({n})" for g, n in genes.items()
                    )
            bullets.append(text + ".")
        parts.append("<ul>" + "".join(f"<li>{b}</li>" for b in bullets) + "</ul>")
        # profile: share of every category of a feature among TP, FN and FP
        features = [("type", "Variant type")]
        if "consequence_group" in rows.columns:
            features += [("consequence_group", "Consequence"), ("impact", "Impact")]
        features += [("vaf_bin", "VAF")]
        if "in_repeat" in rows.columns:
            features.append(("in_repeat", "Repeat"))
        features += [
            (c, c.replace("giab_", "GIAB ")) for c in rows.columns if c.startswith("giab_")
        ]
        features.append(("class", "Main reason"))
        views = []
        for col, label in features:
            names, frames = [], []
            for status in ("TP", "FN", "FP"):
                name, sub, vaf_col = sets[status]
                if col == "vaf_bin":
                    values = (
                        vaf_bin(pd.to_numeric(sub[vaf_col], errors="coerce"))
                        if vaf_col in sub
                        else pd.Series(dtype=str)
                    )
                else:
                    values = sub[col].astype(str) if col in sub else pd.Series(dtype=str)
                names.append(name)
                frames.append(values.value_counts())
            table = pd.DataFrame(frames, index=names).fillna(0)
            if table.empty or table.values.sum() == 0:
                continue
            if col == "class":
                order = [c for c in ["concordant", *CLASS_ORDER] if c in table.columns]
            elif col == "vaf_bin":
                order = [c for c in [*VAF_LABELS, "NA"] if c in table.columns]
            else:
                order = sorted(table.columns)
            share = table.div(table.sum(axis=1).replace(0, 1), axis=0) * 100
            traces = []
            for i, cat in enumerate(order):
                color = (
                    CLASS_COLOR.get(cat, PALETTE[i % len(PALETTE)])
                    if col == "class"
                    else PALETTE[i % len(PALETTE)]
                )
                traces.append(
                    {
                        "type": "bar",
                        "name": CLASS_LABEL.get(cat, cat) if col == "class" else cat,
                        "x": names,
                        "y": share[cat].round(1).tolist(),
                        "customdata": table[cat].astype(int).tolist(),
                        "marker": {"color": color},
                        "hovertemplate": "%{x}<br>"
                        + esc(cat)
                        + ": %{customdata} variants (%{y} %)<extra></extra>",
                    }
                )
            views.append(
                {
                    "label": label,
                    "traces": traces,
                    "layout": bar_layout("% of the variants of the set", ymax=100),
                }
            )
        parts.append("<div class='vc-group'>What characterizes the discordant variants?</div>")
        parts.append(
            mqc.plot(
                views,
                height=380,
                caption="Share of every category among the variants found by both technologies "
                f"(TP), found in {esc(truth)} only (FN) and in {esc(query)} only (FP); hover for "
                "the counts. A feature that is much more frequent among FN or FP than among TP "
                "points to a reason for the difference (tested in the enrichment table below).",
            )
        )
        # VAF distribution: truth VAF of TP and FN, query VAF of FP (one point per variant)
        traces = []
        for status in ("TP", "FN", "FP"):
            name, sub, vaf_col = sets[status]
            if sub.empty or vaf_col not in sub:
                continue
            sub = sub.dropna(subset=[vaf_col])
            traces.append(
                {
                    "type": "box",
                    "name": name,
                    "y": pd.to_numeric(sub[vaf_col], errors="coerce").round(4).tolist(),
                    "s": sub["sample"].tolist(),
                    "text": (sub["variant"] + " · " + sub["sample"]).tolist(),
                    "boxpoints": "all" if len(sub) <= 3000 else "outliers",
                    "jitter": 0.4,
                    "pointpos": 0,
                    "marker": {"color": STATUS_COLOR[status], "size": 4, "opacity": 0.6},
                    "line": {"color": STATUS_COLOR[status]},
                    "hovertemplate": "%{text}<br>VAF %{y}<extra></extra>",
                }
            )
        parts.append("<div class='vc-group'>VAF of concordant and discordant variants</div>")
        parts.append(
            mqc.plot(
                [
                    {
                        "label": "VAF",
                        "traces": traces,
                        "layout": {
                            "yaxis": {"title": {"text": "Tumor VAF"}, "range": [0, 1.02]},
                            "showlegend": False,
                        },
                    }
                ],
                height=360,
                caption=f"Tumor VAF of the {esc(truth)} call for TP and FN, of the "
                f"{esc(query)} call for FP. Low VAFs among the discordant variants point to "
                "the detection limit; hover a point for the variant.",
            )
        )
        # combinations of reasons (a variant can have several)
        disc = rows[rows["status"].isin(["FN", "FP"])]
        if "reasons" in disc.columns and not disc.empty:
            combos = (
                disc.assign(combo=disc["reasons"].astype(str).str.replace(",", " + "))
                .groupby(["combo", "status"])
                .size()
                .unstack(fill_value=0)
            )
            combos = combos.loc[combos.sum(axis=1).sort_values(ascending=False).index].head(12)
            traces = [
                {
                    "type": "bar",
                    "name": sets[s][0],
                    "x": combos.index.tolist(),
                    "y": combos[s].tolist(),
                    "marker": {"color": STATUS_COLOR[s]},
                }
                for s in ("FN", "FP")
                if s in combos.columns
            ]
            parts.append("<div class='vc-group'>Combinations of reasons</div>")
            parts.append(
                mqc.plot(
                    [
                        {
                            "label": "reasons",
                            "traces": traces,
                            "layout": bar_layout(
                                "Discordant variants",
                                extra={
                                    "xaxis": {
                                        "type": "category",
                                        "showgrid": False,
                                        "tickangle": -30,
                                    }
                                },
                            ),
                        }
                    ],
                    height=400,
                    caption="Combinations of the reasons of the discordant variants (column "
                    "reasons: the main reason first, then the further reasons), the 12 most "
                    "frequent.",
                )
            )
        # protein-altering discordant variants (top list)
        if "consequence_group" in disc.columns:
            impact_rank = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
            prot = disc[disc["consequence_group"] == "non_synonymous"].copy()
            prot["rank"] = prot["impact"].map(impact_rank).fillna(4) if "impact" in prot else 4
            prot = prot.sort_values(["rank", "gene", "sample"])
            shown = pd.DataFrame(
                {
                    "variant": prot["variant"],
                    "consequence": prot["consequence"],
                    "impact": prot.get("impact", ""),
                    "sample": prot["sample"],
                    "found in": [
                        (
                            badge(f"{truth} only", STATUS_COLOR["FN"])
                            if s == "FN"
                            else badge(f"{query} only", STATUS_COLOR["FP"])
                        )
                        for s in prot["status"]
                    ],
                    "class": [
                        badge(CLASS_LABEL.get(c, c), CLASS_COLOR.get(c, "#888"))
                        for c in prot["class"]
                    ],
                    "truth VAF": prot["truth_vaf"],
                    "query VAF": prot["query_vaf"],
                }
            )
            parts.append("<div class='vc-group'>Protein-altering discordant variants</div>")
            parts.append(
                html_table(
                    shown,
                    link=mqc.link("discordant-variants"),
                    link_text="All discordant variants in the table browser",
                )
            )
        # genes and enrichment (top lists)
        genes = d["genes"]
        if not genes.empty:
            g = genes[(genes["truth"] == truth) & (genes["query"] == query)].copy()
            g = g[(g["fn"] + g["fp"]) > 0].sort_values(["fn", "fp"], ascending=False)
            if not g.empty:
                parts.append("<div class='vc-group'>Genes with discordant variants</div>")
                cols = [
                    c
                    for c in ("gene", "tp", "fn", "fp", "discordance_rate", "samples", "q_value")
                    if c in g
                ]
                parts.append(
                    html_table(
                        g[cols], link=mqc.link("genes"), link_text="All genes in the table browser"
                    )
                )
        enrich = d["enrichment"]
        if not enrich.empty:
            e = enrich[(enrich["truth"] == truth) & (enrich["query"] == query)].copy()
            if not e.empty:
                e = e.sort_values("p_value")
                e["significant"] = [
                    badge(f"q < {fdr:g}", "#e15759") if q < fdr else ""
                    for q in e["q_value"].fillna(1)
                ]
                cols = [
                    "discordant",
                    "feature",
                    "fraction_discordant",
                    "fraction_concordant",
                    "odds_ratio",
                    "p_value",
                    "q_value",
                    "significant",
                ]
                parts.append(
                    "<div class='vc-group'>Enrichment of features among discordant variants</div>"
                    "<p class='vc-note'>Fisher's exact test of the discordant set (FN, FP or both) "
                    "against the concordant variants (TP), Benjamini-Hochberg q over all tests; "
                    "sorted by p-value.</p>"
                )
                parts.append(
                    html_table(
                        e[[c for c in cols if c in e]],
                        link=mqc.link("enrichment"),
                        link_text="All tests in the table browser",
                    )
                )
    mqc.html(
        "characterization",
        "Discordant variants in detail",
        "**What is shown:** which variants are found in one technology only, and what they "
        "have in common, per comparison.\n\n"
        "- **Lists:** FN (in the truth technology only) and FP (in the query technology only), "
        "with gene, protein change and reasons (at most 15 rows; all in the table browser).\n"
        "- **Category profiles:** the share of indels, low-VAF variants, repeats and "
        "stratifications among TP, FN and FP; **VAF** of the three sets.\n"
        "- **Combinations of reasons**, **protein-altering variants** and **genes** with "
        "discordant variants.\n"
        "- **Enrichment:** features more frequent among discordant than among concordant "
        "variants (Fisher's exact test, Benjamini-Hochberg q).",
        "\n".join(parts),
    )


def section_metrics(mqc, d):
    """Metrics per sample: recall/precision/F1 overall (per method), by variant type and
    by truth VAF bin."""
    pp = d["per_pair"]
    if pp.empty:
        return
    multi = pp["query"].nunique() > 1
    parts = []
    # by variant type: x = sample, bars SNV / INDEL (per query technology)
    views = []
    for metric in ("recall", "precision", "f1"):
        traces = []
        rows = pp[pp["stratum_type"] == "variant_type"]
        for i, ((query, stratum), sub) in enumerate(rows.groupby(["query", "stratum"])):
            sub = sub.sort_values("sample")
            traces.append(
                {
                    "type": "bar",
                    "name": f"{query} {stratum}" if multi else stratum,
                    "x": sub["sample"].tolist(),
                    "y": sub[metric].tolist(),
                    "s": sub["sample"].tolist(),
                    "renameX": True,
                    "marker": {"color": PALETTE[i % len(PALETTE)]},
                    "customdata": [
                        f"TP {fmt_int(a)}, FN {fmt_int(b)}, FP {fmt_int(c)}"
                        for a, b, c in zip(sub["tp_baseline"], sub["fn"], sub["fp"], strict=True)
                    ],
                    "hovertemplate": "%{x}: %{y:.3f}<br>%{customdata}<extra>%{fullData.name}</extra>",
                }
            )
        views.append(
            {
                "label": metric.capitalize(),
                "traces": traces,
                "layout": bar_layout(metric.capitalize(), ymax=1.05, barmode="group"),
            }
        )
    parts.append("<div class='vc-group'>By variant type</div>")
    parts.append(
        mqc.plot(
            views,
            height=380,
            caption="RTG vcfeval per sample; SNV includes atomized MNVs, INDEL insertions, "
            "deletions and complex variants. Hover for the counts.",
        )
    )
    # by truth VAF: one line per sample
    rows = pp[(pp["stratum_type"] == "truth_vaf") & (pp["stratum"] != "unknown")]
    if not rows.empty:
        order = sorted(
            rows["stratum"].unique(),
            key=lambda s: float(re.sub(r"[\[\(\]\)]", "", s).split(",")[0]),
        )
        traces = []
        for i, ((sample, query), sub) in enumerate(rows.groupby(["sample", "query"])):
            sub = sub.set_index("stratum").reindex(order)
            name = f"{sample} ({query})" if multi else sample
            traces.append(
                {
                    "type": "scatter",
                    "mode": "lines+markers",
                    "name": name,
                    "sample": sample,
                    "x": order,
                    "y": sub["recall"].tolist(),
                    "customdata": [
                        f"{fmt_int(a)} of {fmt_int(a + b)}" if not pd.isna(a) else "no variants"
                        for a, b in zip(sub["tp_baseline"], sub["fn"], strict=True)
                    ],
                    "line": {"color": PALETTE[i % len(PALETTE)]},
                    "hovertemplate": esc(name)
                    + "<br>VAF %{x}: recall %{y:.3f} (%{customdata})<extra></extra>",
                }
            )
        parts.append("<div class='vc-group'>Recall by truth VAF</div>")
        parts.append(
            mqc.plot(
                [
                    {
                        "label": "recall",
                        "traces": traces,
                        "layout": {
                            "xaxis": {
                                "type": "category",
                                "title": {"text": "Truth VAF bin (--vaf_bins)"},
                            },
                            "yaxis": {"title": {"text": "Recall"}, "range": [0, 1.05]},
                        },
                    }
                ],
                height=380,
                caption="One line per sample: recall of the truth variants in every VAF bin. "
                "Low recall in the low-VAF bins shows the detection limit of the query.",
            )
        )
    # power-adjusted recall only if it was computed (needs the query tumor alignments)
    if not d["power"].empty and d["power"]["micro_recall"].notna().any():
        p = d["power"]
        parts.append(
            "<p class='vc-note'>Power-adjusted recall (only truth variants the query could "
            "detect at its depth, binomial model, --detection_power): "
            + ", ".join(f"{esc(r['query'])} {fmt_num(r['micro_recall'])}" for _, r in p.iterrows())
            + ".</p>"
        )
    parts.append(
        f"<p class='vc-note'><a href='{esc(mqc.link('metrics-per-sample'))}' target='_blank'>"
        "All metrics and strata per sample in the table browser &rarr;</a></p>"
    )
    mqc.html(
        "metrics",
        "Comparison metrics per sample",
        "**What is shown:** recall, precision and F1 of every sample split by what makes "
        "variants hard to find (RTG vcfeval matching).\n\n"
        "- **By variant type:** SNVs and indels separately; indels are usually harder "
        "(alignment, representation).\n"
        "- **By VAF of the truth call:** recall in VAF bins; low recall in the low bins shows "
        "the detection limit of the query technology (depth, purity).\n"
        "- **Power-adjusted recall** (with alignments): only the truth variants that the "
        "query could detect at its depth.",
        "\n".join(parts),
    )


def flow_html(steps):
    """The region steps as a flow of boxes (sizes in Mbp)."""
    boxes = []
    genome = steps["bp"].iloc[0]
    for i, (_, r) in enumerate(steps.iterrows()):
        removed = (
            f"<span class='text-danger' style='font-size:11px'>-{fmt_bp(r['removed_bp'])}</span>"
            if i and not pd.isna(r.get("removed_bp"))
            else ""
        )
        boxes.append(
            f"<div class='vc-step'><b>{esc(step_label(r['name']))}</b>{fmt_bp(r['bp'])} "
            f"({fmt_pct(r['bp'], genome, 2)}) {removed}</div>"
        )
    return "<div class='vc-flow'>" + "<span class='vc-arrow'>&rarr;</span>".join(boxes) + "</div>"


def section_regions(mqc, d):
    """Compared genomic region: cards (genome, every target BED, common region, callable),
    flow of the steps, vertical bars in Mbp."""
    steps, targets, callable_rows, sheet = d["steps"], d["targets"], d["callable"], d["samplesheet"]
    if steps.empty:
        return
    genome = steps["bp"].iloc[0]
    pairs = comparisons(d["per_pair"])
    truth = pairs[0][0] if pairs else ""
    items = [(fmt_bp(genome), "compared genome", "without chrM, decoy and alt contigs")]
    labels, values, colors, s = [], [], [], []
    labels.append("genome")
    values.append(genome)
    colors.append(PALETTE[0])
    s.append("")
    # every target BED (one per targeted technology)
    if not targets.empty:
        for _, r in targets[targets["type"] == "target"].iterrows():
            tech = r["region"]
            lib = ""
            if not sheet.empty:
                hit = sheet[sheet["technology"] == tech]
                lib = value_of(hit.iloc[0], "library_type") if not hit.empty else ""
            role = "truth" if tech == truth else "query"
            items.append(
                (
                    fmt_bp(r["bp"]),
                    f"{tech} target regions ({role})",
                    f"{lib + ' library, ' if lib else ''}{fmt_pct(r['bp'], genome, 2)} of the genome",
                )
            )
            labels.append(f"{tech} targets")
            values.append(r["bp"])
            colors.append(PALETTE[1])
            s.append("")
    for _, r in steps.iloc[1:].iterrows():
        labels.append(step_label(r["name"]))
        values.append(r["bp"])
        colors.append(PALETTE[0])
        s.append("")
    final = steps["bp"].iloc[-1]
    items.append(
        (fmt_bp(final), "common comparison region", f"{fmt_pct(final, genome, 2)} of the genome")
    )
    if not callable_rows.empty:
        med = callable_rows["callable_bp"].median()
        items.append(
            (
                fmt_bp(med),
                "callable per sample (median)",
                f"{fmt_bp(callable_rows['callable_bp'].min())} - {fmt_bp(callable_rows['callable_bp'].max())}",
            )
        )
        for _, r in callable_rows.sort_values("sample").iterrows():
            labels.append(f"callable {r['sample']}")
            values.append(r["callable_bp"])
            colors.append(PALETTE[2])
            s.append(r["sample"])
    texts = [f"{fmt_bp(v)} ({fmt_pct(v, genome, 1)})" for v in values]
    unit, unit_name = (1e6, "Mbp") if genome >= 1e7 else (1e3, "kbp")
    trace = {
        "type": "bar",
        "x": labels,
        # Mbp for genomes, kbp for small test references (readable log ticks)
        "y": [round(v / unit, 4) for v in values],
        "s": s,
        "text": texts,
        "textposition": "outside",
        "cliponaxis": False,
        "marker": {"color": colors},
        "hovertemplate": "%{x}: %{text}<extra></extra>",
        "showlegend": False,
    }
    views = [
        {
            "label": "log scale",
            "traces": [trace],
            "layout": bar_layout(
                f"Size ({unit_name}, log scale)",
                extra={
                    "yaxis": {
                        "type": "log",
                        "dtick": 1,
                        "title": {"text": f"Size ({unit_name}, log scale)"},
                    }
                },
            ),
        },
        {
            "label": "linear scale",
            "traces": [trace],
            "layout": bar_layout(f"Size ({unit_name})", ymax=genome / unit * 1.15),
        },
    ]
    parts = [
        cards(items, accent="#f28e2b"),
        flow_html(steps),
        mqc.plot(
            views,
            height=400,
            caption="Blue: compared genome and the steps of the region; orange: target BED of "
            "every targeted technology; green: callable region per sample.",
        ),
    ]
    if not callable_rows.empty:
        table = callable_rows.copy()
        table["callable"] = table["callable_bp"].map(fmt_bp)
        table["common region"] = table["common_bp"].map(fmt_bp)
        cols = [
            c
            for c in (
                "sample",
                "alignments",
                "mode",
                "min_coverage",
                "common region",
                "callable",
                "fraction_callable",
            )
            if c in table
        ]
        parts.append(
            html_table(
                table[cols],
                link=mqc.link("region-steps"),
                link_text="Region steps in the table browser",
            )
        )
    mqc.html(
        "regions",
        "Compared genomic region",
        "**What is shown:** the genomic region in which the technologies are compared, and "
        "how it was built.\n\n"
        "- **Cards:** the compared genome (without chrM, decoy and alt contigs), the target "
        "BED of every targeted technology, the common region and the callable region per "
        "sample.\n"
        "- **Flow:** the steps from the genome to the common region with the size after every "
        "step (intersection of the target BEDs, `--regions_bed`, `--exclude_bed`).\n"
        "- **Bars:** the same sizes in Mbp (log or linear scale); green: callable per sample "
        "(bases with enough depth in every alignment of the sample).\n"
        "- **Table:** the callable region per sample (mosdepth mode and depth thresholds).",
        "\n".join(parts),
    )


# the boxes of the filtering funnel: step name of the accounting tables -> (title, what
# the step does); in the order of the normalization chain (variant_accounting.py)
FUNNEL_STEPS = {
    "input": ("Input records", "all records of the input VCFs"),
    "pass_split": ("PASS filter", "only calls with FILTER PASS (or '.') are compared"),
    "contigs": ("Compared contigs", "records on chrM, decoy, alt and unplaced contigs removed"),
    "gvcf_blocks": ("gVCF blocks", "reference blocks (only ALT <NON_REF> or <*>) removed"),
    "normalization": (
        "Normalization",
        "split per ALT, left-aligned, MNVs atomized, REF checked, deduplicated (bcftools norm)",
    ),
    "harmonization": (
        "Harmonization",
        "VAF, depth and alt reads set the same way for every caller; leftover symbolic and "
        "'*' alleles removed",
    ),
    "filters": (
        "Uniform filters and region",
        "the same VAF, depth and alt-read thresholds for every VCF; only the comparison region",
    ),
}

# preprocessing steps that always run but are written to the accounting tables only if they
# removed records (variant_accounting.py): step number -> step name
ALWAYS_STEPS = {3: "contigs", 4: "gvcf_blocks"}


def funnel_steps(acc):
    """The records left after every preprocessing step, summed over the given accounting
    rows (one or more VCFs).

    acc: accounting rows (variant_accounting.py: step, step_name, outcome, reason,
    records, snv, indel). Returns a list of dicts in step order: name, records (left
    after the step), removed ({reason: records}), added (records), and for the input and
    the compared records snv/indel. The last entry is the compared records ('compared',
    the kept records of the last step); 'expected' is the count from the steps before,
    which must be the same (consistency check of the accounting).
    """
    if acc.empty:
        return []
    acc = acc.copy()
    acc["records"] = pd.to_numeric(acc["records"], errors="coerce").fillna(0)
    steps = []
    left = 0
    acc["step"] = acc["step"].astype(int)
    # steps 3 and 4 always run, but have rows only if they removed records
    step_nos = sorted(set(acc["step"]) | set(ALWAYS_STEPS))
    for step_no in step_nos:
        rows = acc[acc["step"] == step_no]
        if rows.empty:
            steps.append(
                {"name": ALWAYS_STEPS[step_no], "records": left, "removed": {}, "added": 0}
            )
            continue
        name = str(rows["step_name"].iloc[0])
        # step 1: the records read from the input VCF
        inp = rows[rows["outcome"] == "in"]
        if not inp.empty:
            left = int(inp["records"].sum())
            steps.append({"name": name, "records": left, "removed": {}, "added": 0, **types(inp)})
            continue
        removed = rows[rows["outcome"] == "removed"].groupby("reason")["records"].sum()
        removed = {str(k): int(v) for k, v in removed.items() if v > 0}
        added = int(rows.loc[rows["outcome"] == "added", "records"].sum())
        left = left - sum(removed.values()) + added
        steps.append({"name": name, "records": left, "removed": removed, "added": added})
    # the compared records: kept by the last step
    last = acc[acc["step"] == acc["step"].max()]
    kept = last[last["outcome"] == "kept"]
    steps.append(
        {
            "name": "compared",
            "records": int(kept["records"].sum()),
            "expected": left,
            "removed": {},
            "added": 0,
            **types(kept),
        }
    )
    return steps


def types(rows):
    """SNV and indel counts of accounting rows (NA if a table gives none)."""
    out = {}
    for col in ("snv", "indel"):
        values = (
            pd.to_numeric(rows[col], errors="coerce") if col in rows else pd.Series(dtype=float)
        )
        out[col] = int(values.sum()) if values.notna().any() else None
    return out


def funnel_html(steps):
    """The steps as a funnel of boxes: width = share of the input records (at least 18 %,
    so that small numbers stay readable), removed records per reason next to every box."""
    if not steps:
        return "<p class='text-muted'>No data.</p>"
    total = max(steps[0]["records"], 1)
    rows = []
    for i, s in enumerate(steps):
        title, what = FUNNEL_STEPS.get(s["name"], ("Compared variants", ""))
        if s["name"] == "compared":
            what = "PASS calls in the comparison region, compared between the technologies"
        width = max(18.0, 100.0 * s["records"] / total)
        kind = ""
        if s.get("snv") is not None and s.get("indel") is not None:
            kind = f"<span class='vc-fkind'>{fmt_int(s['snv'])} SNV · {fmt_int(s['indel'])} indel</span>"
        # what the step changed: removed per reason, added split records
        changes = [
            f"<span class='vc-fremoved' style='border-color:{REASON_COLOR.get(r, '#e15759')}'>"
            f"&minus;{fmt_int(n)} <code>{esc(r)}</code></span>"
            for r, n in sorted(s["removed"].items(), key=lambda kv: -kv[1])
        ]
        if s["added"]:
            changes.append(
                f"<span class='vc-fadded'>+{fmt_int(s['added'])} <code>split</code></span>"
            )
        if i and not changes and s["name"] != "compared":
            changes.append("<span class='vc-note'>no change</span>")
        if s["name"] == "compared" and s.get("expected") not in (None, s["records"]):
            changes.append(
                f"<span class='vc-fremoved'>expected {fmt_int(s['expected'])} from the steps "
                "above: see the consistency checks</span>"
            )
        final = " final" if s["name"] == "compared" else ""
        rows.append(
            f"<div class='vc-frow'><div class='vc-fbarwrap'>"
            f"<div class='vc-fbar{final}' style='width:{width:.1f}%' title='{esc(what)}'>"
            f"<b>{esc(title)}</b> <span class='vc-fn'>{fmt_int(s['records'])}</span>"
            f"<span class='vc-fpct'>{fmt_pct(s['records'], total, 1)}</span>{kind}</div></div>"
            f"<div class='vc-fchanges'>{''.join(changes)}</div></div>"
        )
    return "<div class='vc-funnel'>" + "".join(rows) + "</div>"


def section_filtering(mqc, d):
    """Input variant sets and their filtering: funnel of the preprocessing steps (all VCFs,
    per technology), then the records per VCF, compared and removed per reason (vertical
    stacked bars)."""
    acc = d["accounting"]
    if acc.empty:
        return
    acc = acc.copy()
    # funnel: all VCFs and every technology (summed over the samples)
    groups = [("all", "all VCFs", acc)] + [
        (f"tech_{i}", f"{tech} VCFs", sub)
        for i, (tech, sub) in enumerate(
            sorted(acc.groupby("technology"), key=lambda kv: str(kv[0]))
        )
    ]
    funnel = [toggle_buttons("funnel", [(k, label) for k, label, _ in groups])]
    for key, label, sub in groups:
        n_vcfs = sub[["sample", "technology"]].drop_duplicates().shape[0]
        funnel.append(
            f"<div data-vc-group='funnel' data-vc-key='{key}' class='{'' if key == 'all' else 'vc-hidden'}'>"
            f"<p class='vc-note'>{esc(label)}: {n_vcfs} VCF(s), records summed over the samples.</p>"
            + funnel_html(funnel_steps(sub))
            + "</div>"
        )
    acc["vcf"] = acc["sample"].astype(str) + " " + acc["technology"].astype(str)
    order = sorted(acc["vcf"].unique())
    sample_of = acc.drop_duplicates("vcf").set_index("vcf")["sample"].astype(str)
    removed = acc[acc["outcome"] == "removed"].pivot_table(
        index="vcf", columns="reason", values="records", aggfunc="sum", fill_value=0
    )
    last = acc["step"].max()
    kept = acc[(acc["step"] == last) & (acc["outcome"] == "kept")].groupby("vcf")["records"].sum()
    data = removed.join(kept.rename("compared"), how="outer").fillna(0).reindex(order, fill_value=0)
    data = data[["compared"] + [c for c in data.columns if c != "compared"]]
    views = []
    for label, table in (
        ("records", data),
        ("% of the records", data.div(data.sum(axis=1).replace(0, 1), axis=0) * 100),
    ):
        traces = [
            {
                "type": "bar",
                "name": col,
                "x": order,
                "y": table[col].round(2).tolist(),
                "s": [sample_of[v] for v in order],
                "marker": {"color": REASON_COLOR.get(col, PALETTE[i % len(PALETTE)])},
                "hovertemplate": "%{x}<br>" + esc(col) + ": %{y}<extra></extra>",
            }
            for i, col in enumerate(table.columns)
        ]
        views.append(
            {
                "label": label,
                "traces": traces,
                "layout": bar_layout(
                    label.capitalize(),
                    extra={"xaxis": {"type": "category", "showgrid": False, "tickangle": -45}},
                ),
            }
        )
    reasons = "".join(f"<li><code>{esc(k)}</code>: {esc(v)}</li>" for k, v in REASON_TEXT.items())
    body = [
        "<div class='vc-group'>From the input records to the compared variants</div>",
        *funnel,
        "<div class='vc-group'>Records of every VCF</div>",
        mqc.plot(
            views,
            height=440,
            caption="Input records per VCF (sample and technology): compared (green) and removed "
            "per reason; 'split' additions are not shown. Buttons: records or percent.",
        ),
        f"<details><summary>What the reasons mean</summary><ul>{reasons}</ul></details>",
        f"<p class='vc-note'><a href='{esc(mqc.link('preprocessing'))}' target='_blank'>Counts per "
        f"step</a> and <a href='{esc(mqc.link('removed-records'))}' target='_blank'>every removed "
        "record</a> in the table browser &rarr;</p>",
    ]
    mqc.html(
        "filtering",
        "Input variant sets and their filtering",
        "**What is shown:** how many records of the input VCFs reach the comparison, and why "
        "the others were removed.\n\n"
        "- **Funnel:** one box per preprocessing step, in the order the pipeline runs them; "
        "the box width is the share of the input records that is left (at least 18 % wide so "
        "that small numbers stay readable). Next to every box: the records the step removed, "
        "per reason, and the records it added (`split`: one record per ALT allele). The "
        "buttons show all VCFs or the VCFs of one technology.\n"
        "- **Records of every VCF:** the same counts for every input VCF as stacked bars.\n"
        "- **How to read it:** the comparison only covers the green part. Most input records "
        "are usually removed by the **PASS filter** (unfiltered VCFs carry many rejected "
        "candidates) and by the **comparison region** (calls outside the common targets or the "
        "callable bases).",
        "\n".join(body),
    )


def section_coverage(mqc, d):
    """Coverage: mean depth and percent of bases >= thresholds per alignment, in the common
    region and in the target region (vertical bars)."""
    cov = d["coverage"]
    if cov.empty:
        return
    cov = cov.sort_values(["sample", "technology", "role"])
    order = list(dict.fromkeys(cov["alignment"]))
    sample_of = cov.drop_duplicates("alignment").set_index("alignment")["sample"].astype(str)
    views = []
    metrics = [("mean_depth", "Mean depth")] + [
        (c, c.replace("pct_ge_", "% bases >= ")) for c in cov.columns if c.startswith("pct_ge_")
    ]
    for col, label in metrics:
        traces = []
        for i, (region, sub) in enumerate(cov.groupby("region")):
            sub = sub.set_index("alignment").reindex(order)
            traces.append(
                {
                    "type": "bar",
                    "name": f"{region} region",
                    "x": order,
                    "y": sub[col].tolist(),
                    "s": [sample_of[a] for a in order],
                    "marker": {"color": PALETTE[i % len(PALETTE)]},
                }
            )
        views.append(
            {
                "label": label,
                "traces": traces,
                "layout": bar_layout(
                    label,
                    ymax=105 if col.startswith("pct_") else None,
                    barmode="group",
                    extra={"xaxis": {"type": "category", "showgrid": False, "tickangle": -45}},
                ),
            }
        )
    mqc.html(
        "coverage",
        "Coverage",
        "**What is shown:** the depth of every alignment (mosdepth), the main reason why a "
        "technology misses variants.\n\n"
        "- **Two regions per alignment:** the common comparison region (where variants are "
        "compared) and the target region of its technology (the genome for WGS).\n"
        "- **Buttons:** mean depth, or the percent of bases with at least the depths of "
        "`--coverage_thresholds`.\n"
        "- A low depth in the comparison region explains *low coverage* and *subthreshold "
        "evidence* reasons in chapter 3.",
        mqc.plot(views, height=420)
        + f"<p class='vc-note'><a href='{esc(mqc.link('coverage'))}' target='_blank'>Coverage table "
        "in the table browser &rarr;</a></p>",
    )


def section_checks(mqc, d):
    """Consistency checks: counts per category and status, failed checks listed."""
    checks = d["checks"]
    if checks.empty:
        return
    summary = (
        checks.assign(ok=checks["status"] == "ok")
        .groupby("category")
        .agg(checks=("status", "size"), passed=("ok", "sum"))
        .reset_index()
    )
    summary["failed"] = summary["checks"] - summary["passed"]
    failed = checks[checks["status"] != "ok"]
    parts = [
        html_table(
            summary,
            link=mqc.link("consistency-checks"),
            link_text="All checks in the table browser",
        ),
    ]
    if not failed.empty:
        parts.append("<div class='vc-group text-danger'>Failed checks</div>")
        parts.append(
            html_table(failed, link=mqc.link("consistency-checks"), link_text="All checks")
        )
    mqc.html(
        "checks",
        "Consistency checks",
        "**What is shown:** whether the counts of all steps add up, as a table of check "
        "categories (all checks in the table browser).\n\n"
        "- **Accounting:** records per VCF and step, input = removed + kept + added.\n"
        "- **Metrics and tables:** the metrics equal the counts of the variant tables, the "
        "aggregated metrics the sum over the samples, every isec status its metrics.\n"
        "- **Discordance and reasons:** every discordant variant has exactly one main reason, "
        "and every reason rule is recomputed; at least `--evidence_check_min_fraction` of "
        "the TP variants have alt reads in the alignments (read evidence works).\n"
        "- **Regions and required variants:** sizes and statuses agree between the steps.\n\n"
        "A failed check stops the pipeline after the reports were written, so that the "
        "failure can be inspected.",
        "\n".join(parts),
    )


def chapter_intros(d, args):
    """The introduction of every chapter (markdown): what the chapter shows, what was done,
    and which section answers which question. Each ends with a sentence (MultiQC appends a
    full stop to a module introduction)."""
    pairs = comparisons(d["per_pair"])
    truth = pairs[0][0] if pairs else param(d["params"], "truth_technology")
    queries = ", ".join(sorted({q for _, q in pairs})) or "NA"
    n_samples = d["per_pair"]["sample"].nunique() if not d["per_pair"].empty else 0
    known = not (d["required_summary"].empty and d["known_benchmark"].empty)
    pass_only = args.filter_pass_only != "false"
    return {
        "summary": (
            "This chapter gives the **main result of the run**: how well every query technology "
            f"(**{queries}**) finds the somatic variants of the truth technology **{truth}** in "
            f"the same tumor/normal pairs ({n_samples} samples).\n\n"
            "- **Recall** = TP / (TP + FN): the share of the truth variants that the query "
            "technology also called.\n"
            "- **Precision** = TP / (TP + FP): the share of the query calls that the truth "
            "technology confirms.\n"
            "- **F1**: the harmonic mean of recall and precision.\n"
            "- **Matching methods**: RTG vcfeval (haplotype-aware, the main result), bcftools "
            "isec with the same allele and with the same position; the buttons switch between "
            "them. Different numbers show variants that the technologies represent "
            "differently.\n\n"
            "Read **1.1** for the overall numbers and the key findings, **1.2** for what was "
            "compared, **1.3-1.5** for every sample and **1.6-1.7** for the overlap of the "
            "variant sets. Chapter 3 explains the variants found by one technology only."
        ),
        "known": (
            (
                "Validation with **variants known beforehand**, independent of the comparison "
                "of the technologies (`--known_variants`):\n\n"
                "- **Required variants** (a clinical checklist): every variant must be called "
                "with FILTER PASS inside the comparison region by every technology in every "
                "sample; FAIL states why not.\n"
                "- **Benchmark** (an external truth set such as SEQC2 or GIAB): recall and "
                "precision of every technology within the confident regions of the set.\n\n"
                "These results are as important for a validation as the comparison of the "
                "technologies."
            )
            if known
            else (
                "Validation with **variants known beforehand** (a clinical checklist of "
                "required variants, or an external truth set such as SEQC2 or GIAB). No "
                "known-variant list was given in this run; see `--known_variants` in the "
                "documentation."
            )
        ),
        "discordance": (
            "Variants called by only one technology are **discordant**: **FN** are in the "
            f"truth technology {truth} only (missed by the query), **FP** are in the query "
            "technology only. Every discordant variant is explained from the view of the "
            "technology that **misses** it:\n\n"
            "- **Main reason** (exactly one): called there but not compared (FILTER not PASS, "
            "or removed by the uniform filters or the region), or not called: no alignment "
            "given, low coverage, alt reads below the caller's threshold, or no alt reads.\n"
            "- **Further reasons** (any number): germline evidence in the normal, and the "
            "context: indel, low VAF, repeat region, GIAB and further stratifications.\n"
            "- **Read evidence**: reads with the variant allele are counted in the alignments "
            "of the technology that misses it (indels by their haplotype sequence).\n\n"
            "**3.1** shows the reasons per sample, **3.2** the variants themselves and what "
            "characterizes them (type, VAF, genes, enriched features)."
            + (
                f"\n\n**IGV views**: every compared variant can be inspected in the reads of the "
                f"tumor and normal alignments of every technology ([IGV views]({args.igv_link}); "
                "column *igv* of the variant tables)."
                if args.igv_links
                else ""
            )
        ),
        "filtering": (
            "Before the comparison, **every input VCF runs through the same steps**, so that "
            "differences between the technologies come from the data and not from the formats "
            "of the callers:\n\n"
            + (
                "1. **PASS filter**: only PASS calls are compared; the other calls are kept "
                "aside to explain discordant variants.\n"
                if pass_only
                else "1. **No PASS filter** (`--filter_pass_only false`): all calls are compared.\n"
            )
            + "2. **Contigs and gVCF blocks**: only the compared reference contigs (no chrM, "
            "decoy or alt contigs); reference blocks are not variants.\n"
            "3. **Normalization** (bcftools norm): one record per ALT allele, left-aligned, "
            "MNVs split into SNVs, REF checked against the FASTA, duplicates removed.\n"
            "4. **Harmonization**: VAF, depth and alt reads computed the same way for every "
            "caller; leftover symbolic and '*' alleles removed.\n"
            f"5. **Uniform filters and region**: VAF &ge; {args.min_vaf}, depth &ge; "
            f"{args.min_depth}, alt reads &ge; {args.min_alt_reads} (0 = off), and only the "
            "comparison region (chapter 5).\n\n"
            "**4.1** shows how many records every step removed, **4.2-4.3** the normalization "
            "and harmonization of every VCF, **4.4** the annotation of the compared variants "
            "by Ensembl VEP (if it was run)."
        ),
        "regions": (
            "Variants are only compared **where every technology could have called them**:\n\n"
            "- the **intersection of the target BEDs** of all targeted technologies "
            "(whole-genome sequencing has none), optionally only `--regions_bed` and without "
            "`--exclude_bed`;\n"
            + (
                "- per sample, only the **callable bases**: enough depth in the tumor (and "
                "normal) alignment of every technology (mosdepth).\n\n"
                if args.callable_regions == "true"
                else "- callable regions were not used in this run (`--callable_regions false`).\n\n"
            )
            + "**5.1** shows how the region was built and how large it is, **5.2** the depth of "
            "every alignment in the compared and in its target region, **5.3** the depth "
            "distributions (mosdepth)."
        ),
        "methods": (
            "How the results were produced, and **whether the counts can be trusted**:\n\n"
            "- **6.1 Consistency checks**: the counts of all steps must add up; a failed check "
            "stops the run after the reports were written.\n"
            "- **6.2 Input checks**: every input was checked against the reference at the start "
            "(contig names and lengths, genome build, REF alleles).\n"
            "- **6.3 Methods description** with the tools and their citations, **6.4 software "
            "versions** of every step, **6.5 workflow summary** with the parameters of the "
            "run.\n\n"
            "Cite the tools of 6.3 when you use these results."
        ),
    }


def section_known_placeholder(mqc):
    """Known variants without a list: a short note, so that the chapter numbers stay the
    same in every run."""
    mqc.html(
        "known",
        "Required variants and external truth sets",
        "No known-variant list was given (`--known_variants`).",
        "<p class='vc-note'>A known-variant sheet adds a clinical checklist of required "
        "variants (PASS/FAIL per sample and technology) or a benchmark against an external "
        "truth set (docs/usage.md, section Known variants).</p>",
    )


def section_coverage_placeholder(mqc):
    """Coverage without alignments: a short note (stable chapter numbers)."""
    mqc.html(
        "coverage",
        "Coverage",
        "No alignments were given, so the coverage was not computed.",
        "<p class='vc-note'>With tumor (and normal) alignments in the samplesheet, the "
        "pipeline computes the depth of every alignment, the callable regions and the read "
        "evidence of discordant variants.</p>",
    )


def write_multiqc(args, d):
    """Write all variantcomp MultiQC sections, chapter by chapter."""
    mqc = Mqc(args.multiqc, args.tables_link, chapter_intros(d, args))
    # 1 Summary
    section_summary(mqc, d)
    section_introduction(mqc, d, args)
    section_recall_precision(mqc, d)
    section_metrics(mqc, d)
    section_overlap(mqc, d)
    section_truth_found(mqc, d)
    # 2 Known variants
    if d["required_summary"].empty and d["known_benchmark"].empty:
        section_known_placeholder(mqc)
    else:
        section_known(mqc, d)
    # 3 Discordance
    section_why(mqc, d)
    section_characterization(mqc, d, args.fdr)
    # 4 Variant filtering
    section_filtering(mqc, d)
    # 5 Coverage and genomic regions
    section_regions(mqc, d)
    if d["coverage"].empty:
        section_coverage_placeholder(mqc)
    else:
        section_coverage(mqc, d)
    # 6 Methods and quality control
    section_checks(mqc, d)
    n = len([f for f in os.listdir(args.multiqc) if f.endswith(("_mqc.html", "_mqc.json"))])
    print(f"INFO: MultiQC: {n} sections, {mqc.n_plots} plots in {args.multiqc}/", file=sys.stderr)


# ----------------------------------------------------------------------------- datavzrd

# yes/no columns: green/red where yes is the goal (found, passed), blue/grey for
# descriptive flags (in a repeat, in a stratum)
YES_NO = {"scale": "ordinal", "domain": ["yes", "no"], "range": ["#a1d99b", "#fcae91"]}
FLAG = {"scale": "ordinal", "domain": ["yes", "no"], "range": ["#9ecae1", "#eeeeee"]}
IMPACT = {
    "scale": "ordinal",
    "domain": ["HIGH", "MODERATE", "LOW", "MODIFIER"],
    "range": ["#e15759", "#f28e2b", "#edc948", "#d9d9d9"],
}
# columns shown in the detail view of a row (expand a row) instead of the table
DETAIL_COLUMNS = [
    "hgvsc",
    "transcript",
    "biotype",
    "consequence_group",
    "variant_class",
    "repeat_unit",
    "repeat_times",
    "truth_alt",
    "query_alt",
    "other_in_normalized",
    "other_vaf_normalized",
    "normal_alt",
    "query_power",
    "chrom",
    "pos",
]
# first columns of the variant tables: sample, gene and protein change, status and
# technologies, VAF/depth of the truth, query and normal, then the reasons and the
# consequence; the column variant (gene + protein change in one text, the same as gene and
# hgvsp) comes last (VARIANT_BACK)
VARIANT_FRONT = [
    "sample",
    "gene",
    "hgvsp",
    "status",
    "class",
    "truth",
    "query",
    "truth_vaf",
    "query_vaf",
    "normal_vaf",
    "truth_dp",
    "query_dp",
    "normal_dp",
    "reasons",
    "consequence",
    "impact",
    "locus",
    "ref",
    "alt",
    "type",
    "length",
    "query_filter",
    "other_filter",
    "other_tumor_dp",
    "other_tumor_alt",
    "detectable",
    "in_repeat",
    "dbsnp",
    "cosmic",
    "clin_sig",
]
VARIANT_BACK = ["variant", "variant_id"]


def heatmap(spec):
    """A datavzrd heatmap column."""
    return {"plot": {"heatmap": dict(spec)}}


def column_specs(df, genome):
    """datavzrd column configuration of a table: colours, links, detail view."""
    ensembl = (
        "https://grch37.ensembl.org"
        if genome.upper() in ("GRCH37", "HG19")
        else "https://www.ensembl.org"
    )
    specs = {}
    classes = ["concordant", *CLASS_ORDER, "."]
    for col in df.columns:
        values = set(df[col].dropna().astype(str).unique())
        if col == "status" and values <= {"TP", "FN", "FP"}:
            specs[col] = heatmap(
                {
                    "scale": "ordinal",
                    "domain": ["TP", "FN", "FP"],
                    "range": ["#a1d99b", "#fdae6b", "#d4b9da"],
                }
            )
        elif col == "status":
            # PASS / FAIL_* (required variants), ok / failed (checks)
            good = {"PASS", "ok", "yes"}
            domain = sorted(values)
            specs[col] = heatmap(
                {
                    "scale": "ordinal",
                    "domain": domain,
                    "range": ["#a1d99b" if v in good else "#fc9272" for v in domain],
                }
            )
        elif col == "class" or col.startswith("class_"):
            domain = [c for c in classes if c in values] or classes
            specs[col] = heatmap(
                {
                    "scale": "ordinal",
                    "domain": domain,
                    "range": [CLASS_COLOR.get(c, "#eeeeee") for c in domain],
                }
            )
        elif col == "reasons" or col.startswith("reasons_"):
            # one coloured pill per reason; the domain is every reason in the table (main
            # reasons, further reasons, stratifications)
            domain = sorted({r for v in values for r in v.split(",") if r not in ("", ".")})
            specs[col] = {
                "plot": {
                    "pills": {
                        "separator": ",",
                        "domain": domain or ["concordant"],
                        "range": [CLASS_COLOR.get(r, STRATUM_COLOR) for r in domain]
                        or [CLASS_COLOR["concordant"]],
                    }
                }
            }
        elif col.startswith("found_") or col in ("detectable", "all_passed", "called", "in_region"):
            # yes is the goal: green / red
            if values <= {"yes", "no", "NA", "."}:
                specs[col] = heatmap(YES_NO)
        elif values and values <= {"yes", "no", "NA", "."}:
            # descriptive flags (in_repeat, giab_*, stratifications): blue / grey
            specs[col] = heatmap(FLAG)
        elif col == "impact":
            specs[col] = heatmap(IMPACT)
        elif col == "gene":
            specs[col] = {
                "link-to-url": {
                    "GeneCards": {
                        "url": "https://www.genecards.org/cgi-bin/carddisp.pl?gene={value}"
                    },
                    "Ensembl": {"url": ensembl + "/Homo_sapiens/Gene/Summary?g={value}"},
                }
            }
        elif col == "igv":
            # the IGV page of the variant, opened at its row (#variant_id); the cell shows
            # 'IGV' instead of the file name (empty for rows without a page)
            specs[col] = {
                "link-to-url": {
                    "IGV": {"url": "../../igv/{value}#{variant_id}", "new-window": True},
                    "custom-content": "function(value, row) { return value ? 'IGV' : ''; }",
                }
            }
        elif col == "locus":
            specs[col] = {
                "link-to-url": {
                    "Ensembl": {"url": ensembl + "/Homo_sapiens/Location/View?r={value}"},
                }
            }
        elif col == "dbsnp":
            specs[col] = {
                "plot": {"pills": {"separator": ",", "range": ["#4e79a7"]}},
                "link-to-url": {"dbSNP": {"url": "https://www.ncbi.nlm.nih.gov/snp/{pill-value}"}},
            }
        elif col == "cosmic":
            specs[col] = {
                "plot": {"pills": {"separator": ",", "range": ["#b07aa1"]}},
                "link-to-url": {
                    "COSMIC": {"url": "https://cancer.sanger.ac.uk/cosmic/search?q={pill-value}"}
                },
            }
        elif col in ("truth_vaf", "query_vaf", "vaf", "normal_vaf", "other_vaf_normalized"):
            specs[col] = {"plot": {"ticks": {"scale": "linear", "domain": [0, 1]}}}
        elif (
            col in ("recall", "precision", "f1")
            or col.startswith(("micro_", "macro_"))
            and not col.endswith("_sd")
        ):
            specs[col] = heatmap(
                {
                    "scale": "linear",
                    "domain": [0, 0.5, 1],
                    "range": ["#fc9272", "#ffffbf", "#a1d99b"],
                }
            )
        elif col == "q_value":
            specs[col] = heatmap(
                {
                    "scale": "linear",
                    "domain": [0, 0.05, 1],
                    "range": ["#fc9272", "#fee0d2", "#ffffff"],
                }
            )
        # less important columns of wide tables go into the detail view of a row
        if col in DETAIL_COLUMNS or (col.startswith("giab_") and len(df.columns) > 25):
            specs.setdefault(col, {})["display-mode"] = "detail"
        # the variant name stays visible when scrolling to the right
        if col == "variant":
            specs.setdefault(col, {})["display-mode"] = "pinned"
    return specs


def view_desc(title, what, how=(), links=(), multiqc_link=None, df=None, details=None):
    """Markdown description at the top of a datavzrd view: what the table is, how to read
    it, links to the MultiQC report and the related tables, then collapsed: further
    explanations (details = (summary, items)) and the column glossary."""
    lines = [f"## {title}", "", what]
    if how:
        lines += ["", *[f"- {h}" for h in how]]
    nav = []
    if multiqc_link:
        nav.append(f"[MultiQC report]({multiqc_link})")
    nav += [f"[{text}](../{view}/index_1.html)" for view, text in links]
    if nav:
        lines += ["", "See also: " + " · ".join(nav)]
    # long explanations (e.g. of the discordance classes) are collapsed
    if details:
        summary, items = details
        lines += ["", f"<details><summary>{summary}</summary>", ""]
        lines += [f"- {item}" for item in items]
        lines += ["", "</details>"]
    if df is not None:
        known = [c for c in df.columns if c in COLUMN_TEXT]
        if known:
            lines += ["", "<details><summary>Columns</summary>", ""]
            lines += [f"- **{c}**: {COLUMN_TEXT[c]}" for c in known]
            lines += ["", "</details>"]
    return "\n".join(lines)


MULTIQC_PAGE_JS = """// variantcomp: show the MultiQC report inside the datavzrd report (datavzrd render-html)
var canvas = document.getElementById("canvas");
var src = "%s";
// the canvas centres its content in a row: show the text above the report instead
canvas.style.width = "100%%";
canvas.style.display = "block";
canvas.innerHTML =
  "<p style='margin:8px 0'>The MultiQC report of this run (main report: summary, plots, QC). " +
  "<a href='" + src + "' target='_blank'>Open it in a new tab</a>.</p>" +
  "<iframe src='" + src + "' style='width:100%%;height:85vh;border:1px solid #ddd;border-radius:4px'></iframe>";
"""


def add_igv_column(df, links):
    """The column igv (IGV page of the variant, from the links table of --igv_report) after
    the gene (or the sample); unchanged without links or without a variant_id column. Rows whose
    variant_id has no page (e.g. removed records) get an empty cell."""
    if links.empty or "variant_id" not in df.columns:
        return df
    pages = dict(zip(links["variant_id"].astype(str), links["igv_file"].astype(str), strict=True))
    df = df.copy()
    igv = df["variant_id"].astype(str).map(pages).fillna("")
    # right after the gene (or the sample), so that it is visible without scrolling
    after = [c for c in ("sample", "gene") if c in df.columns]
    pos = df.columns.get_loc(after[-1]) + 1 if after else 0
    df.insert(pos, "igv", igv)
    return df


def write_datavzrd(args, d):
    """Write the datavzrd tables, the MultiQC page script and the configuration."""
    tables, datasets, views = {}, {}, {}
    mq = args.multiqc_link

    def add(
        name,
        df,
        title,
        what,
        how=(),
        links=None,
        see=(),
        front=VARIANT_FRONT,
        page_size=25,
        excel=False,
        details=None,
    ):
        """One table: file, dataset and view (columns ordered, configured)."""
        if df.empty:
            return
        df = order_columns(df, front, back=VARIANT_BACK)
        df = add_igv_column(df, d["igv_links"])
        # the main reason is the first entry of reasons: no separate column class
        df = df.drop(columns=[c for c in ("class",) if c in df.columns and "reasons" in df.columns])
        for col in df.columns:
            values = df[col].dropna()
            # counts that pandas read as floats (because of missing values) as integers
            if df[col].dtype.kind == "f" and len(values) and (values % 1 == 0).all():
                df[col] = df[col].astype("Int64")
            # '.' = no ID: an empty cell instead of a pill with a dot
            if col in ("dbsnp", "cosmic", "clin_sig"):
                df[col] = df[col].replace(".", "")
        tables[f"{name}.tsv"] = df
        datasets[name] = {
            "path": f"{name}.tsv",
            "separator": "\t",
            **({"links": links} if links else {}),
            **({"offer-excel": True} if excel else {}),
        }
        views[name] = {
            "dataset": name,
            "desc": view_desc(title, what, how, see, mq, df, details),
            "page-size": page_size,
            "render-table": {"columns": column_specs(df, args.genome)},
        }

    variants = add_variant_columns(d["variants"])
    truth = add_variant_columns(d["truth"])
    truth_link = {
        "column": "variant_id",
        "table-row": "truth-variants/variant_id",
        "optional": True,
    }
    gene_link = {"column": "gene", "table-row": "genes/gene", "optional": True}
    # the MultiQC report as the first page
    reports = pd.DataFrame([{"report": "MultiQC", "path": mq}])
    tables["reports.tsv"] = reports
    datasets["reports"] = {"path": "reports.tsv", "separator": "\t"}
    views["multiqc-report"] = {
        "dataset": "reports",
        "desc": "## MultiQC report\n\nThe main report of the run: summary, plots and QC of every step.",
        "render-html": {"script-path": "multiqc_page.js"},
    }
    status_how = [
        "**status**: TP found by both technologies (green), FN in the truth technology only = missed by the query (orange), FP in the query technology only (purple).",
        "**reasons**: every reason of a variant. A discordant variant has its main reason first (exactly one: the call and the reads of the technology that misses the variant), a concordant one `concordant`; then germline evidence in the normal and the context (indel, low VAF, repeat, stratifications). Filter reasons for e.g. low_vaf or giab_ to see which hard variants were found by both technologies (TP) and which were missed.",
        "Click a row to see the detail columns (HGVS.c, transcript, biotype, evidence counts).",
        "Links: gene to GeneCards/Ensembl, locus to the Ensembl genome browser, dbSNP and COSMIC IDs (with a VEP cache).",
    ]
    add(
        "truth-variants",
        truth,
        "Truth variants",
        "Every variant of the truth technology that was compared, per sample: gene and protein "
        "change first, then whether each query technology found it (found_&lt;query&gt;, green "
        "yes / red no) and its reasons (reasons_&lt;query&gt;: if missed, the main reason "
        "first; if found, concordant first; then the context such as low_vaf or giab_*).",
        how=[
            "Sort by found_by to see the variants found by no query technology: they may be artefacts of the truth technology.",
            "Filter a column with the funnel icon in its header, e.g. impact = HIGH.",
            *status_how[1:],
        ],
        links={
            "discordant variant": {
                "column": "variant_id",
                "table-row": "discordant-variants/variant_id",
                "optional": True,
            },
            "gene table": gene_link,
        },
        see=[("discordant-variants", "Discordant variants"), ("genes", "Genes")],
        front=[
            "sample",
            "gene",
            "hgvsp",
            "found_by",
            *[c for c in truth.columns if c.startswith(("found_", "reasons_"))],
            *VARIANT_FRONT,
        ],
        excel=True,
    )
    if not variants.empty:
        disc = variants[variants["status"].isin(["FN", "FP"])]
        # explanations of the main reasons and the further reasons (collapsed)
        class_list = [
            f"main reason **{c}** ({CLASS_LABEL[c]}): {CLASS_TEXT[c]}" for c in CLASS_ORDER
        ]
        class_list += [f"further reason **{c}**: {FURTHER_TEXT[c]}" for c in FURTHER_ORDER]
        class_list += [f"further reason **giab_&lt;name&gt;** and stratifications: {STRATUM_TEXT}"]
        add(
            "discordant-variants",
            disc,
            "Discordant variants",
            "Every variant found by one technology only: FN = in the truth technology only "
            "(missed by the query), FP = in the query technology only. Each is explained from "
            "the view of the technology that misses it.",
            how=status_how,
            details=("Reasons", class_list),
            links={"truth variant": truth_link, "gene table": gene_link},
            see=[
                ("truth-variants", "Truth variants"),
                ("compared-variants", "All compared variants"),
                ("genes", "Genes"),
            ],
            excel=True,
        )
        add(
            "compared-variants",
            variants,
            "All compared variants",
            "Every variant of every comparison with its status (TP, FN, FP), the discordance "
            "columns and the annotation.",
            how=status_how,
            links={"truth variant": truth_link, "gene table": gene_link},
            see=[
                ("truth-variants", "Truth variants"),
                ("discordant-variants", "Discordant variants"),
            ],
            excel=True,
        )
    add(
        "genes",
        d["genes"],
        "Genes",
        "Concordant (tp) and discordant (fn, fp) variants per comparison and gene (Ensembl VEP), "
        "with Fisher's exact test of the discordance rate of the gene against all other genes "
        "(Benjamini-Hochberg q, red: low q).",
        see=[("discordant-variants", "Discordant variants"), ("enrichment", "Enrichment")],
        front=["gene", "truth", "query", "tp", "fn", "fp", "discordance_rate", "q_value"],
    )
    add(
        "enrichment",
        d["enrichment"],
        "Enrichment",
        "Are discordant variants (FN, FP, FN+FP) enriched for a feature compared with the "
        "concordant ones? Fisher's exact test per comparison and feature, Benjamini-Hochberg q "
        "over all tests (red: low q).",
        see=[("genes", "Genes")],
        front=["truth", "query", "discordant", "feature", "q_value"],
    )
    add(
        "metrics-per-sample",
        d["per_pair"],
        "Metrics per sample",
        "Recall, precision and F1 of every sample and technology pair: overall (RTG vcfeval), "
        "per matching method (match: exact = bcftools isec same allele, position = same "
        "position), per variant type and per VAF bin of the truth and the query call.",
        front=["sample", "truth", "query", "stratum_type", "stratum", "recall", "precision", "f1"],
    )
    req_var = add_variant_columns(d["required_variants"])
    add(
        "required-variants",
        req_var,
        "Required variants",
        "PASS/FAIL of every required (clinical checklist) variant per sample and technology. "
        "PASS: called with FILTER PASS in the comparison region. FAIL_filtered: called but "
        "not PASS; FAIL_outside_region; FAIL_not_called (with read evidence counts if "
        "alignments were given).",
        links={"truth variant": truth_link},
        front=[
            "label",
            "gene",
            "status",
            "sample",
            "technology",
            "locus",
            "ref",
            "alt",
            "filter",
            "vaf",
            "depth",
            "alt_reads",
            "in_region",
            "called",
        ],
    )
    add(
        "known-benchmark",
        d["known_benchmark_pairs"],
        "Benchmark against the external truth set",
        "Recall, precision and F1 of every technology and sample against the external truth "
        "set (known variants, mode benchmark), within its confident regions.",
        front=["sample", "truth", "query", "stratum_type", "stratum", "recall", "precision", "f1"],
    )
    removed = d["removed"]
    add(
        "removed-records",
        add_variant_columns(removed) if not removed.empty else removed,
        "Records removed during preprocessing",
        "Every record of the input VCFs that was not compared, with the step and the reason.",
        how=[f"**{k}**: {v}" for k, v in REASON_TEXT.items() if k != "split"],
        see=[("preprocessing", "Counts per step")],
        front=["sample", "technology", "step_name", "reason", "locus", "ref", "alt", "filter"],
    )
    add(
        "preprocessing",
        d["accounting"],
        "Preprocessing accounting",
        "Records per input VCF and step: kept, removed (with the reason) and added (splitting "
        "of multiallelic records).",
        see=[("removed-records", "Removed records")],
        front=["sample", "technology", "step", "step_name", "outcome", "reason", "records"],
    )
    steps = d["steps"].copy()
    if not steps.empty:
        steps.insert(2, "size", steps["bp"].map(fmt_bp))
    add(
        "region-steps",
        steps,
        "Comparison region",
        "Size of the comparison region after every step: the compared genome, the "
        "intersection of the target BEDs, --regions_bed, --exclude_bed.",
        front=["step", "name", "size"],
    )
    add(
        "coverage",
        d["coverage"],
        "Coverage",
        "Mean depth and percent of the bases with at least the depth of each threshold, per "
        "alignment, in the common region and in its target region (the genome for WGS).",
        front=["alignment", "sample", "technology", "role", "region"],
    )
    add(
        "consistency-checks",
        d["checks"],
        "Consistency checks",
        "Checks that the counts of the pipeline steps add up (green ok, red failed).",
        front=["status", "category", "scope", "check", "expected", "observed"],
    )
    sheet = d["samplesheet"]
    add(
        "samplesheet",
        sheet,
        "Samplesheet",
        "The validated input: one row per sample and technology.",
        front=["sample", "technology", "library_type", "caller", "is_truth"],
    )
    os.makedirs(args.datavzrd, exist_ok=True)
    for name, df in tables.items():
        # missing values as empty cells (datavzrd shows them blank, not as NaN)
        df.to_csv(os.path.join(args.datavzrd, name), sep="\t", index=False, na_rep="")
    with open(os.path.join(args.datavzrd, "multiqc_page.js"), "w") as fh:
        fh.write(MULTIQC_PAGE_JS % mq)
    config = {
        "name": "variantcomp: variant tables",
        # the first page is the MultiQC report; the tables are in the menu
        "default-view": "multiqc-report",
        "datasets": datasets,
        "views": views,
    }
    # JSON is valid YAML: no YAML library is needed
    with open(os.path.join(args.datavzrd, "datavzrd.yaml"), "w") as fh:
        json.dump(config, fh, indent=2)
    print(f"INFO: datavzrd: {len(tables)} tables, {len(views)} views", file=sys.stderr)


# ----------------------------------------------------------------------------- main


def load(args):
    """Read every input table."""
    d = {
        "per_pair": read(args.per_pair),
        "metrics": read(args.metrics),
        "power": read(args.power),
        "variants": read(args.variants),
        "truth": read(args.truth_variants),
        "classes": read(args.classes),
        "enrichment": read(args.enrichment),
        "genes": read(args.genes),
        "annotation_counts": read(args.annotation_counts),
        # the accounting folder also holds the *.checks.tsv tables: keep accounting tables
        "accounting": read_many(args.accounting, required=["step", "outcome", "reason", "records"]),
        "removed": read_many(args.removed),
        "steps": read(args.region_steps),
        "targets": read(args.targets),
        "callable": read(args.callable),
        "coverage": read(args.coverage),
        "required_summary": read(args.required_summary),
        "required_variants": read(args.required_variants),
        "known_benchmark": read(args.known_benchmark),
        "known_benchmark_pairs": read(args.known_benchmark_per_pair),
        "checks": read(args.checks),
        "samplesheet": read(args.samplesheet, sep=","),
        "params": read(args.params),
        "igv_links": read(args.igv_links),
    }
    # VAFs and metric values as numbers
    for col in ("truth_vaf", "query_vaf"):
        if col in d["variants"].columns:
            d["variants"][col] = pd.to_numeric(d["variants"][col], errors="coerce")
    return d


def main():
    """Parse the command line and write the requested outputs."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for name in (
        "per-pair",
        "metrics",
        "power",
        "variants",
        "truth-variants",
        "classes",
        "enrichment",
        "genes",
        "annotation-counts",
        "region-steps",
        "targets",
        "callable",
        "coverage",
        "required-summary",
        "required-variants",
        "known-benchmark",
        "known-benchmark-per-pair",
        "checks",
        "samplesheet",
        "params",
        "igv-links",
    ):
        parser.add_argument(f"--{name}")
    parser.add_argument("--accounting", nargs="*", default=[])
    parser.add_argument("--removed", nargs="*", default=[])
    parser.add_argument("--multiqc", help="output: folder of the MultiQC custom-content files")
    parser.add_argument("--datavzrd", help="output: folder of the datavzrd tables and config")
    parser.add_argument("--pipeline-version", default="unknown")
    parser.add_argument("--run-name", default="unknown")
    parser.add_argument("--start", default="unknown")
    parser.add_argument("--genome", default="GRCh38", help="build of the Ensembl links")
    parser.add_argument("--fdr", type=float, default=0.05, help="q-value threshold (marks rows)")
    # settings of the run shown in the Introduction (pipeline params of the same name; the
    # parameter table only lists parameters that differ from their defaults)
    parser.add_argument("--filter-pass-only", default="true", help="--filter_pass_only")
    parser.add_argument("--min-vaf", default="0", help="--min_vaf")
    parser.add_argument("--min-depth", default="0", help="--min_depth")
    parser.add_argument("--min-alt-reads", default="0", help="--min_alt_reads")
    parser.add_argument("--callable-regions", default="false", help="--callable_regions")
    parser.add_argument(
        "--tables-link",
        default="../report/tables",
        help="datavzrd report folder relative to the MultiQC report",
    )
    parser.add_argument(
        "--multiqc-link",
        default="../../../multiqc/multiqc_report.html",
        help="MultiQC report relative to a datavzrd view page (tables/<view>/index_1.html)",
    )
    parser.add_argument(
        "--igv-link",
        default="../report/igv/index.html",
        help="IGV overview page relative to the MultiQC report (with --igv-links)",
    )
    args = parser.parse_args()
    if not args.multiqc and not args.datavzrd:
        sys.exit("ERROR: nothing to write: give --multiqc and/or --datavzrd")
    if not 0 < args.fdr < 1:
        sys.exit(f"ERROR: --fdr {args.fdr} must be between 0 and 1")
    d = load(args)
    if d["per_pair"].empty:
        print("WARNING: no per-pair metrics: the report has no comparison results", file=sys.stderr)
    if args.datavzrd:
        write_datavzrd(args, d)
    if args.multiqc:
        write_multiqc(args, d)


if __name__ == "__main__":
    main()
