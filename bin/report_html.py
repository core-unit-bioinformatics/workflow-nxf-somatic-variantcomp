#!/usr/bin/env python3
"""The variantcomp report: a narrative HTML report and the datavzrd tables it links to.

Outputs (one call writes both)
  --html      the narrative report (one offline HTML file): Summary, Introduction, Results,
              Methods; plots are inline SVG, there is no external JavaScript
  --datavzrd  folder with the tables and the configuration (datavzrd.yaml) of the datavzrd
              report: browsable, linked tables of every variant, truth variant, discordance
              class, removed record, region, coverage and check

The explanations of the discordance classes, filtering reasons and columns (TEXTS below)
are shared by both outputs, so the two reports describe the data the same way.

INPUT (all tables are TSV with header; optional inputs may be left out)
  --per-pair, --metrics       benchmark/metrics.per_pair.tsv, metrics.aggregated.tsv
  --power                     benchmark/power.aggregated.tsv (optional)
  --variants                  discordance/variants.annotated.tsv (every compared variant)
  --truth-variants            report/<prefix>.truth_variants.tsv (variantcomp_report.py)
  --classes, --strata, --enrichment, --genes, --annotation-counts   report tables
  --accounting, --removed     vcf_preprocessed/accounting/*.accounting.tsv / *.removed.tsv
  --region-steps, --targets   regions/common_targets.steps.tsv / .summary.tsv
  --callable, --coverage      coverage/callable_summary.tsv, coverage/coverage.coverage.tsv
  --required-summary, --required-variants, --known-benchmark   known_variants/ (optional)
  --checks                    report/consistency_checks.tsv
  --samplesheet               pipeline_info/samplesheet.validated.csv
  --params                    parameter table (group, parameter, value)
  --versions                  pipeline_info/*_software_mqc_versions.yml
  --pipeline-version, --run-name, --start   run information
OUTPUT
  --html      e.g. variantcomp_report.html
  --datavzrd  <folder>/datavzrd.yaml and the tables it names (<folder>/<view>.tsv; the paths
              in the configuration are file names: datavzrd runs next to the tables)

Why a script: the report combines the tables of all steps with explanatory text and plots;
no tool writes such a report from these tables (MultiQC shows QC tables and plots, datavzrd
browsable tables).
"""

import argparse
import html
import io
import json
import os
import re
import sys
from datetime import datetime

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

# ----------------------------------------------------------------------------- texts

CLASS_ORDER = [
    "filtered_non_pass",
    "present_but_filtered",
    "germline_in_normal",
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
    "germline_in_normal": (
        "The normal alignments of the sample show the alternative allele (VAF >= "
        "--germline_vaf with at least --evidence_min_alt_reads reads). The variant is probably "
        "germline (or a contamination of the normal); one caller rejected it as germline, the "
        "other did not."
    ),
    "low_coverage": (
        "The tumor alignment of the technology that misses the variant has a depth below "
        "--min_coverage_tumor at this position: the variant could not be detected reliably. "
        "Typical for positions at the edge of or outside a capture target."
    ),
    "subthreshold_evidence": (
        "The tumor alignment of the technology that misses the variant has reads with the "
        "alternative allele, but the caller did not call it (or the call failed a threshold). "
        "Typical for low-VAF variants near the detection limit of the technology."
    ),
    "no_evidence": (
        "The tumor alignment of the technology that misses the variant has enough depth but "
        "no read with the alternative allele (counted with MQ >= --evidence_min_mapq and "
        "BQ >= --evidence_min_baseq). For a false negative this suggests a variant below the "
        "detection limit or a sampling difference; for a false positive an artefact of the "
        "technology that called it (sequencing or alignment error)."
    ),
    "absent": (
        "Not called and no read evidence available to look further (no alignment given, or "
        "--read_evidence false)."
    ),
}
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
    "sample": "sample (one tumor/normal pair)",
    "truth": "truth (baseline) technology",
    "query": "query technology compared with the truth",
    "status": "TP: found by both; FN: truth variant missed by the query; FP: query call without truth variant",
    "class": "discordance class (see the class sections); concordant for TP",
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
    "gene": "gene symbol (Ensembl VEP)",
    "consequence": "most severe consequence (VEP, Sequence Ontology term)",
    "consequence_group": "synonymous, non_synonymous or other (from the consequence)",
    "impact": "VEP impact: HIGH, MODERATE, LOW, MODIFIER",
    "in_repeat": "inside a repeat (microsatellite) region of the reference",
    "found_by": "number of query technologies that found the truth variant",
    "variant_id": "sample:chrom:pos:ref>alt (link key between the tables)",
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
PALETTE = ["#2a6f97", "#e07a5f", "#81b29a", "#f2cc8f", "#6d597a", "#3d405b", "#b56576", "#4d908e"]

# ----------------------------------------------------------------------------- helpers


def read(path, sep="\t"):
    """A table as a data frame ('NA' as missing); an empty frame if path is None."""
    if not path:
        return pd.DataFrame()
    return pd.read_csv(path, sep=sep, na_values=["NA"], keep_default_na=False, low_memory=False)


def read_many(paths):
    """Concatenation of several tables with the same header (empty frame if none)."""
    frames = [read(p) for p in paths or []]
    frames = [f for f in frames if not f.empty]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def esc(value):
    """HTML-escaped text."""
    return html.escape(str(value))


def fmt_int(value):
    """Integer with thousands separators."""
    return f"{int(value):,}"


def fmt_number(value):
    """A float of a table: whole numbers as integers (counts that pandas turned into floats
    because of missing values), others with 4 significant digits."""
    return f"{value:.0f}" if float(value).is_integer() else f"{value:.4g}"


def fmt_pct(part, whole, digits=1):
    """part / whole as a percentage string ('NA' if whole is 0)."""
    return f"{100 * part / whole:.{digits}f} %" if whole else "NA"


def fmt_bp(value):
    """A size in bp, kb or Mb."""
    value = int(value)
    if value >= 1e6:
        return f"{value / 1e6:.2f} Mb"
    return f"{value / 1e3:.1f} kb" if value >= 1e3 else f"{value} bp"


def step_label(name):
    """Readable label of a region step (names like 2_targets_intersection)."""
    return STEP_TEXT.get(re.sub(r"^\d+_", "", str(name)), str(name))


def add_variant_id(df):
    """The column variant_id (sample:chrom:pos:ref>alt): link key between the tables."""
    if df.empty or not {"sample", "chrom", "pos", "ref", "alt"} <= set(df.columns):
        return df
    df = df.copy()
    ids = df["sample"].astype(str) + ":" + df["chrom"].astype(str) + ":" + df["pos"].astype(str)
    df.insert(0, "variant_id", ids + ":" + df["ref"].astype(str) + ">" + df["alt"].astype(str))
    return df


def svg(fig):
    """A matplotlib figure as an inline SVG string."""
    buf = io.StringIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    plt.close(fig)
    # drop the XML header and the fixed size, so the figure scales with the page
    text = re.sub(r"^.*?<svg", "<svg", buf.getvalue(), flags=re.S)
    return re.sub(r'(<svg[^>]*?) width="[^"]*" height="[^"]*"', r'\1 width="100%"', text, count=1)


def table_html(df, max_rows=None, glossary=True, caption=None):
    """An HTML table (first max_rows rows) with an optional column glossary."""
    if df is None or df.empty:
        return "<p class='na'>No rows.</p>"
    shown = df if max_rows is None else df.head(max_rows)
    parts = []
    if caption:
        parts.append(f"<p class='caption'>{caption}</p>")
    parts.append(
        "<div class='tablewrap'>"
        + shown.to_html(index=False, na_rep="NA", border=0, classes="data", float_format=fmt_number)
        + "</div>"
    )
    if max_rows is not None and len(df) > max_rows:
        parts.append(
            f"<p class='na'>Showing {max_rows} of {fmt_int(len(df))} rows; all rows in the table browser.</p>"
        )
    if glossary:
        known = [c for c in df.columns if c in COLUMN_TEXT]
        if known:
            items = "".join(f"<dt>{esc(c)}</dt><dd>{esc(COLUMN_TEXT[c])}</dd>" for c in known)
            parts.append(
                f"<details><summary>Columns</summary><dl class='glossary'>{items}</dl></details>"
            )
    return "".join(parts)


# ----------------------------------------------------------------------------- plots


def plot_pr_scatter(per_pair):
    """Recall vs. precision, one dot per sample, coloured by query technology."""
    overall = per_pair[(per_pair["stratum_type"] == "overall") & (per_pair["stratum"] == "all")]
    if overall.empty:
        return ""
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    for i, (query, rows) in enumerate(overall.groupby("query")):
        ax.scatter(
            rows["recall"],
            rows["precision"],
            s=60,
            color=PALETTE[i % len(PALETTE)],
            label=f"{rows['truth'].iloc[0]} vs {query}",
            edgecolor="black",
            linewidth=0.5,
        )
        for _, r in rows.iterrows():
            ax.annotate(
                str(r["sample"]),
                (r["recall"], r["precision"]),
                fontsize=7,
                xytext=(4, 3),
                textcoords="offset points",
            )
    low = min(overall["recall"].min(), overall["precision"].min(), 0.8)
    ax.set_xlim(max(0, low - 0.05), 1.01)
    ax.set_ylim(max(0, low - 0.05), 1.01)
    ax.set_xlabel("Recall (truth variants found)")
    ax.set_ylabel("Precision (query calls that are true)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")
    return svg(fig)


def plot_truth_found(truth):
    """Percent of the truth variants found by each query technology, per sample."""
    found_cols = [c for c in truth.columns if c.startswith("found_") and c != "found_by"]
    if truth.empty or not found_cols:
        return ""
    rows = []
    for col in found_cols:
        for sample, sub in truth.groupby("sample"):
            known = sub[col].isin(["yes", "no"])
            if known.any():
                rows.append(
                    (
                        col.removeprefix("found_"),
                        sample,
                        100 * (sub.loc[known, col] == "yes").mean(),
                    )
                )
    data = pd.DataFrame(rows, columns=["query", "sample", "pct"])
    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    queries = sorted(data["query"].unique())
    for i, q in enumerate(queries):
        sub = data[data["query"] == q]
        ax.bar(i, sub["pct"].mean(), color=PALETTE[i % len(PALETTE)], alpha=0.6)
        ax.scatter([i] * len(sub), sub["pct"], color="black", s=15, zorder=3)
    ax.set_xticks(range(len(queries)), queries)
    ax.set_ylim(0, 105)
    ax.set_ylabel("Truth variants found (%)")
    ax.set_title("bars: mean over samples, dots: samples", fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    return svg(fig)


def plot_classes(classes):
    """Discordant variants per comparison and status, stacked by class."""
    if classes.empty:
        return ""
    labels = (classes["truth"] + " vs " + classes["query"] + " " + classes["status"]).tolist()
    fig, ax = plt.subplots(figsize=(7, 0.5 * len(labels) + 1.2))
    left = [0] * len(labels)
    for i, cls in enumerate(c for c in CLASS_ORDER if c in classes.columns):
        values = classes[cls].fillna(0).tolist()
        ax.barh(labels, values, left=left, color=PALETTE[i % len(PALETTE)], label=cls)
        left = [a + b for a, b in zip(left, values, strict=True)]
    ax.set_xlabel("Discordant variants (all samples)")
    ax.legend(fontsize=7, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.invert_yaxis()
    return svg(fig)


def plot_region_steps(steps, callable_rows):
    """bp of the comparison region after every step, and the callable bp per sample."""
    if steps.empty:
        return ""
    names = [f"{r['step']}. {step_label(r['name'])}" for _, r in steps.iterrows()]
    values = steps["bp"].tolist()
    for _, r in callable_rows.iterrows():
        names.append(f"callable: {r['sample']}")
        values.append(r["callable_bp"])
    fig, ax = plt.subplots(figsize=(7, 0.45 * len(names) + 1))
    colors = [PALETTE[0]] * len(steps) + [PALETTE[2]] * len(callable_rows)
    ax.barh(names, values, color=colors)
    ax.set_xscale("log")
    genome = steps["bp"].iloc[0]
    for i, v in enumerate(values):
        ax.text(v, i, f"  {fmt_int(v)} bp ({fmt_pct(v, genome, 2)})", va="center", fontsize=7)
    ax.set_xlabel("bp (log scale); in brackets: percent of the compared genome")
    ax.invert_yaxis()
    ax.set_xlim(right=genome * 30)
    return svg(fig)


def plot_filtering(accounting):
    """Records per VCF: compared and removed per reason (stacked bars)."""
    if accounting.empty:
        return ""
    acc = accounting.copy()
    acc["vcf"] = acc["sample"].astype(str) + " " + acc["technology"].astype(str)
    removed = acc[acc["outcome"] == "removed"].pivot_table(
        index="vcf", columns="reason", values="records", aggfunc="sum", fill_value=0
    )
    kept = acc[(acc["step"] == 7) & (acc["outcome"] == "kept")].set_index("vcf")["records"]
    data = removed.join(kept.rename("compared"), how="outer").fillna(0)
    data = data[["compared"] + [c for c in data.columns if c != "compared"]]
    fig, ax = plt.subplots(figsize=(7.5, 0.45 * len(data) + 1.2))
    left = [0] * len(data)
    for i, col in enumerate(data.columns):
        ax.barh(data.index, data[col], left=left, color=PALETTE[i % len(PALETTE)], label=col)
        left = [a + b for a, b in zip(left, data[col], strict=True)]
    ax.set_xlabel("Records (input VCF; 'split' additions not shown)")
    ax.legend(fontsize=7, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.invert_yaxis()
    return svg(fig)


def plot_metrics_by_type(metrics):
    """Micro recall and precision by variant type per comparison."""
    rows = metrics[metrics["stratum_type"] == "variant_type"]
    if rows.empty:
        return ""
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.2), sharey=True)
    for ax, col, title in (
        (axes[0], "micro_recall", "Recall"),
        (axes[1], "micro_precision", "Precision"),
    ):
        labels = (rows["query"] + " " + rows["stratum"]).tolist()
        ax.bar(
            range(len(rows)),
            rows[col],
            color=[PALETTE[0] if s == "SNV" else PALETTE[1] for s in rows["stratum"]],
        )
        ax.set_xticks(range(len(rows)), labels, rotation=45, ha="right", fontsize=8)
        ax.set_title(title, fontsize=9)
        ax.set_ylim(0, 1.05)
        ax.grid(axis="y", alpha=0.3)
    return svg(fig)


def plot_recall_by_vaf(metrics):
    """Micro recall by truth VAF bin per comparison."""
    rows = metrics[(metrics["stratum_type"] == "truth_vaf") & (metrics["stratum"] != "unknown")]
    if rows.empty:
        return ""
    fig, ax = plt.subplots(figsize=(6, 3.2))
    for i, (query, sub) in enumerate(rows.groupby("query")):
        ax.plot(
            sub["stratum"],
            sub["micro_recall"],
            marker="o",
            color=PALETTE[i % len(PALETTE)],
            label=query,
        )
    ax.set_xlabel("Truth VAF bin")
    ax.set_ylabel("Recall")
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right", fontsize=8)
    return svg(fig)


def plot_annotation(counts, annotation):
    """Proportion of each value of an annotation per technology and variant set."""
    rows = counts[counts["annotation"] == annotation]
    if rows.empty:
        return ""
    data = rows.groupby(["set", "technology", "value"])["count"].sum().reset_index()
    data["label"] = data["set"] + ": " + data["technology"]
    table = data.pivot_table(index="label", columns="value", values="count", fill_value=0)
    table = table.div(table.sum(axis=1), axis=0) * 100
    fig, ax = plt.subplots(figsize=(7.5, 0.4 * len(table) + 1.2))
    left = [0] * len(table)
    for i, col in enumerate(table.columns):
        ax.barh(table.index, table[col], left=left, color=PALETTE[i % len(PALETTE)], label=col)
        left = [a + b for a, b in zip(left, table[col], strict=True)]
    ax.set_xlabel(f"{annotation} (% of the variants of the set)")
    ax.legend(fontsize=7, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.invert_yaxis()
    return svg(fig)


def plot_coverage(coverage):
    """Mean depth per alignment in the common region and the target region."""
    if coverage.empty:
        return ""
    data = coverage.pivot_table(index="alignment", columns="region", values="mean_depth")
    fig, ax = plt.subplots(figsize=(7, 0.4 * len(data) + 1.2))
    height = 0.8 / len(data.columns)
    for i, col in enumerate(data.columns):
        ax.barh(
            [j + i * height for j in range(len(data))],
            data[col],
            height=height,
            color=PALETTE[i % len(PALETTE)],
            label=col,
        )
    ax.set_yticks([j + 0.4 - height / 2 for j in range(len(data))], data.index, fontsize=8)
    ax.set_xlabel("Mean depth")
    ax.legend(fontsize=8)
    ax.invert_yaxis()
    return svg(fig)


# ----------------------------------------------------------------------------- datavzrd


def view(dataset, desc, page_size=25, columns=None):
    """A datavzrd table view."""
    spec = {"dataset": dataset, "desc": desc, "page-size": page_size}
    spec["render-table"] = {"columns": columns or {}}
    return spec


def glossary_md(df):
    """Markdown list of the known columns of a table."""
    items = [f"- **{c}**: {COLUMN_TEXT[c]}" for c in df.columns if c in COLUMN_TEXT]
    return "\n\nColumns:\n\n" + "\n".join(items) if items else ""


def datavzrd(args):
    """Write the datavzrd tables and configuration."""
    tables, datasets, views = {}, {}, {}
    variants = add_variant_id(read(args.variants))
    truth = add_variant_id(read(args.truth_variants))
    per_pair = read(args.per_pair)
    link_to_truth = {
        "truth variant": {
            "column": "variant_id",
            "table-row": "truth-variants/variant_id",
            "optional": True,
        }
    }

    def add(name, df, desc, links=None, page_size=25):
        """One table: file, dataset and view."""
        if df.empty:
            return
        tables[f"{name}.tsv"] = df
        datasets[name] = {
            "path": f"{name}.tsv",
            "separator": "\t",
            **({"links": links} if links else {}),
        }
        views[name] = view(name, desc + glossary_md(df), page_size)

    add(
        "truth-variants",
        truth,
        "## Truth variants\n\nEvery variant of the truth technology that was compared, per "
        "sample, with its annotation and, for every query technology, whether it was found "
        "(found_<query>) and, if missed, the discordance class (class_<query>).",
    )
    add(
        "metrics-per-pair",
        per_pair,
        "## Metrics per sample\n\nRecall, precision and F1 of every sample and technology "
        "pair, overall and per stratum (variant type, VAF bins, match type).",
    )
    add(
        "compared-variants",
        variants,
        "## All compared variants\n\nEvery variant of every comparison with its status (TP, "
        "FN, FP), the discordance columns and the annotation. FN and FP rows link to the "
        "truth variant table.",
        links=link_to_truth,
    )
    if not variants.empty:
        for cls in CLASS_ORDER:
            sub = variants[variants["class"] == cls]
            add(
                f"discordant-{cls.replace('_', '-')}",
                sub,
                f"## Discordant: {cls}\n\n{CLASS_TEXT[cls]}",
                links=link_to_truth,
            )
    add(
        "removed-records",
        read_many(args.removed),
        "## Records removed during preprocessing\n\nEvery record of the input VCFs that was "
        "not compared, with the step and the reason:\n\n"
        + "\n".join(f"- **{k}**: {v}" for k, v in REASON_TEXT.items() if k != "split"),
    )
    add(
        "preprocessing",
        read_many(args.accounting),
        "## Preprocessing accounting\n\nRecords per input VCF and step: kept, removed (with "
        "the reason) and added (splitting of multiallelic records).",
    )
    add(
        "region-steps",
        read(args.region_steps),
        "## Comparison region\n\nbp of the comparison region after every step: the compared "
        "genome, the intersection of the target BEDs, --regions_bed, --exclude_bed.",
    )
    add(
        "coverage",
        read(args.coverage),
        "## Coverage\n\nMean depth and percent of the bases with at least the depth of each "
        "threshold, per alignment, in the common region and in its target region (the genome "
        "for WGS).",
    )
    add(
        "enrichment",
        read(args.enrichment),
        "## Enrichment\n\nAre discordant variants (FN, FP, FN+FP) enriched for a feature "
        "compared with the concordant ones? Fisher's exact test per comparison and feature, "
        "Benjamini-Hochberg q over all tests.",
    )
    add(
        "genes",
        read(args.genes),
        "## Genes\n\nConcordant and discordant variants per comparison and gene (Ensembl "
        "VEP), with Fisher's exact test of the gene against all other genes.",
    )
    add(
        "required-variants",
        read(args.required_variants),
        "## Required known variants\n\nPASS/FAIL of every required variant per sample and "
        "technology.",
    )
    add(
        "consistency-checks",
        read(args.checks),
        "## Consistency checks\n\nChecks that the counts of the pipeline steps add up.",
    )
    os.makedirs(args.datavzrd, exist_ok=True)
    for name, df in tables.items():
        df.to_csv(os.path.join(args.datavzrd, name), sep="\t", index=False, na_rep="NA")
    config = {
        "name": "variantcomp: tables",
        "default-view": "truth-variants" if "truth-variants" in views else next(iter(views)),
        "datasets": datasets,
        "views": views,
    }
    # JSON is valid YAML: no YAML library is needed
    with open(os.path.join(args.datavzrd, "datavzrd.yaml"), "w") as fh:
        json.dump(config, fh, indent=2)
    print(f"INFO: datavzrd: {len(tables)} tables, {len(views)} views", file=sys.stderr)


# ----------------------------------------------------------------------------- html

CSS = """
:root { --fg: #1d2731; --muted: #5c6b7a; --accent: #2a6f97; --bg: #ffffff; --soft: #f3f6f9;
        --ok: #2e7d32; --fail: #c62828; }
body { font-family: -apple-system, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; color: var(--fg);
       background: var(--bg); margin: 0; line-height: 1.5; }
nav { position: fixed; top: 0; left: 0; width: 230px; height: 100vh; overflow-y: auto;
      background: var(--soft); padding: 16px; box-sizing: border-box; font-size: 13px; }
nav a { color: var(--fg); text-decoration: none; display: block; padding: 2px 0; }
nav a.l2 { font-weight: 600; margin-top: 6px; }
nav a.l3 { padding-left: 12px; color: var(--muted); }
main { margin-left: 250px; max-width: 1050px; padding: 16px 24px 80px; }
h1 { color: var(--accent); margin-bottom: 4px; } h2 { border-bottom: 2px solid var(--soft);
padding-bottom: 4px; margin-top: 40px; } h3 { margin-top: 28px; }
.meta { color: var(--muted); font-size: 14px; }
.banner { padding: 10px 14px; border-radius: 6px; margin: 16px 0; font-weight: 600; }
.banner.ok { background: #e8f5e9; color: var(--ok); } .banner.fail { background: #ffebee; color: var(--fail); }
.cards { display: flex; flex-wrap: wrap; gap: 12px; margin: 12px 0; }
.card { background: var(--soft); border-radius: 8px; padding: 10px 14px; min-width: 150px; }
.card .v { font-size: 22px; font-weight: 700; color: var(--accent); } .card .l { font-size: 12px; color: var(--muted); }
.figs { display: flex; flex-wrap: wrap; gap: 16px; } .fig { flex: 1 1 420px; }
.fig p { font-size: 13px; color: var(--muted); }
table.data { border-collapse: collapse; font-size: 12px; } .tablewrap { overflow-x: auto; max-height: 520px; overflow-y: auto; }
table.data th { background: var(--soft); position: sticky; top: 0; text-align: left; }
table.data th, table.data td { padding: 3px 8px; border-bottom: 1px solid #e3e8ee; white-space: nowrap; }
.na, .caption { color: var(--muted); font-size: 13px; } dl.glossary dt { font-weight: 600; }
dl.glossary dd { margin: 0 0 6px 16px; color: var(--muted); font-size: 13px; }
.flow { display: flex; flex-wrap: wrap; align-items: center; gap: 6px; margin: 12px 0; }
.step { background: var(--soft); border-left: 4px solid var(--accent); border-radius: 4px; padding: 6px 10px; font-size: 12px; }
.step b { display: block; font-size: 13px; } .arrow { color: var(--muted); font-size: 18px; }
.removed { color: var(--fail); font-size: 11px; display: block; }
code { background: var(--soft); padding: 1px 4px; border-radius: 3px; }
"""


class Report:
    """Collects the sections of the HTML report."""

    def __init__(self):
        self.parts, self.toc = [], []

    def h(self, level, anchor, title):
        """A heading with an entry in the navigation."""
        self.toc.append((level, anchor, title))
        self.parts.append(f"<h{level} id='{anchor}'>{esc(title)}</h{level}>")

    def p(self, text):
        """A paragraph (text may contain HTML)."""
        self.parts.append(f"<p>{text}</p>")

    def add(self, fragment):
        """Raw HTML."""
        self.parts.append(fragment)

    def fig(self, svg_text, caption):
        """A figure with its caption."""
        if svg_text:
            self.parts.append(f"<div class='fig'>{svg_text}<p>{caption}</p></div>")

    def render(self, title):
        """The complete HTML page."""
        nav = "".join(
            f"<a class='l{lvl}' href='#{a}'>{esc(t)}</a>" for lvl, a, t in self.toc if lvl <= 3
        )
        return (
            f"<!doctype html><html lang='en'><head><meta charset='utf-8'>"
            f"<meta name='viewport' content='width=device-width, initial-scale=1'>"
            f"<title>{esc(title)}</title><style>{CSS}</style></head><body>"
            f"<nav><b>Contents</b>{nav}</nav><main>{''.join(self.parts)}</main></body></html>"
        )


def overall_rows(df):
    """The overall/all rows of a metrics table."""
    if df.empty:
        return df
    return df[(df["stratum_type"] == "overall") & (df["stratum"] == "all")]


def read_versions(path):
    """{process: {tool: version}} of the collated versions YAML of the pipeline (two levels:
    'PROCESS:' lines, then indented 'tool: version' lines)."""
    versions, process = {}, None
    if not path:
        return versions
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            key, _, value = line.strip().partition(":")
            # an unindented line starts a process; indented lines are its tools
            if not line[0].isspace():
                process = key.strip("'\"")
                versions[process] = {}
            elif process is not None:
                versions[process][key.strip("'\"")] = value.strip().strip("'\"")
    return versions


def section_summary(rep, d):
    """Summary: key numbers and key plots."""
    rep.h(2, "summary", "Summary")
    samples = sorted(d["per_pair"]["sample"].unique()) if not d["per_pair"].empty else []
    overall = overall_rows(d["metrics"])
    steps = d["steps"]
    cards = [(len(samples), "samples (tumor/normal pairs)")]
    if not overall.empty:
        cards.append((overall["truth"].iloc[0], "truth technology"))
        cards.append((", ".join(sorted(overall["query"].unique())), "query technologies"))
    if not steps.empty:
        cards.append(
            (
                fmt_bp(steps["bp"].iloc[-1]),
                f"compared region ({fmt_pct(steps['bp'].iloc[-1], steps['bp'].iloc[0], 2)} of the genome)",
            )
        )
    if not d["truth"].empty:
        cards.append((fmt_int(len(d["truth"])), "truth variants compared"))
    rep.add(
        "<div class='cards'>"
        + "".join(
            f"<div class='card'><div class='v'>{esc(v)}</div><div class='l'>{esc(label)}</div></div>"
            for v, label in cards
        )
        + "</div>"
    )
    if not overall.empty:
        summary = overall[
            [
                "truth",
                "query",
                "n_samples",
                "tp_baseline",
                "fn",
                "tp_call",
                "fp",
                "micro_recall",
                "micro_precision",
                "micro_f1",
            ]
        ].rename(
            columns={"micro_recall": "recall", "micro_precision": "precision", "micro_f1": "f1"}
        )
        rep.add(
            table_html(
                summary,
                caption="Overall result over all samples (counts summed over the samples; recall, precision and F1 from the summed counts).",
            )
        )
    rep.add("<div class='figs'>")
    rep.fig(
        plot_pr_scatter(d["per_pair"]),
        "<b>Recall and precision per sample.</b> Each dot is one sample compared with the truth technology; the closer to the upper right corner, the more similar the technologies.",
    )
    rep.fig(
        plot_truth_found(d["truth"]),
        "<b>Truth variants found.</b> Percent of the compared truth variants found by each query technology (bars: mean over the samples, dots: samples).",
    )
    rep.add("</div>")
    rep.fig(
        plot_classes(d["classes"]),
        "<b>Why variants differ.</b> Discordant variants (FN: missed by the query, FP: only called by the query) by discordance class, over all samples. The classes are explained in the section Discordant variants.",
    )


def section_introduction(rep, d):
    """Introduction: question, data, metrics."""
    rep.h(2, "introduction", "Introduction")
    truth_tech = (
        overall_rows(d["metrics"])["truth"].iloc[0]
        if not overall_rows(d["metrics"]).empty
        else "the truth technology"
    )
    rep.p(
        "This report compares somatic SNV and indel calls of the <b>same tumor/normal pairs</b> "
        "sequenced with <b>different technologies</b> (e.g. whole-genome, whole-exome and panel "
        f"sequencing). The calls of <b>{esc(truth_tech)}</b> serve as the truth for every sample; "
        "every other technology (query) is compared with it. Because the data, not the callers, "
        "differ, all VCFs are normalized and filtered in the same way and compared only in the "
        "genomic region that every technology covers. The question is how many of the truth "
        "variants each technology finds, how many of its calls are confirmed, and why the "
        "remaining variants differ."
    )
    rep.p(
        "<b>Recall</b> is the fraction of the truth variants that the query technology also "
        "called (TP / (TP + FN)); <b>precision</b> the fraction of the query calls that are "
        "truth variants (TP / (TP + FP)); <b>F1</b> their harmonic mean. Variants are matched "
        "haplotype-aware with RTG vcfeval, so different representations of the same variant "
        "match. Every discordant variant (FN, FP) gets a <b>class</b> that explains the "
        "difference from the point of view of the technology that misses it."
    )
    if not d["samplesheet"].empty:
        cols = [
            c
            for c in (
                "sample",
                "technology",
                "library_type",
                "caller",
                "tumor_id",
                "normal_id",
                "target_bed",
            )
            if c in d["samplesheet"].columns
        ]
        sheet = d["samplesheet"][cols].copy()
        if "target_bed" in sheet:
            sheet["target_bed"] = sheet["target_bed"].astype(str).str.split("/").str[-1]
        rep.add(
            table_html(
                sheet,
                glossary=False,
                caption="Input: one row per sample and technology (samplesheet).",
            )
        )


def flow_html(steps):
    """The region steps as a flow of boxes."""
    boxes = []
    genome = steps["bp"].iloc[0]
    for i, (_, r) in enumerate(steps.iterrows()):
        name = step_label(r["name"])
        removed = (
            f"<span class='removed'>-{fmt_int(r['removed_bp'])} bp</span>"
            if i and str(r.get("removed_bp", "NA")) != "NA"
            else ""
        )
        boxes.append(
            f"<div class='step'><b>{esc(name)}</b>{fmt_int(r['bp'])} bp ({fmt_pct(r['bp'], genome, 2)}){removed}</div>"
        )
    return "<div class='flow'>" + "<span class='arrow'>&rarr;</span>".join(boxes) + "</div>"


def section_regions(rep, d):
    """Results: the compared region."""
    rep.h(3, "regions", "Compared genomic region")
    rep.p(
        "Variants are only compared where every technology could have called them. The "
        "comparison region is built step by step: it starts with the compared contigs of the "
        "reference genome (without the mitochondrial genome, decoy and alternative contigs), "
        "is restricted to the intersection of the target BEDs of all targeted technologies "
        "(whole-genome sequencing has no target BED and does not restrict it), optionally to "
        "--regions_bed and without --exclude_bed. With alignments, each sample is further "
        "restricted to the bases with enough depth in every alignment of the sample "
        "(callable region)."
    )
    if not d["steps"].empty:
        rep.add(flow_html(d["steps"]))
    rep.fig(
        plot_region_steps(d["steps"], d["callable"]),
        "<b>Size of the region after each step</b> (blue) and of the callable region per sample (green), in bp (log scale) and in percent of the compared genome.",
    )
    if not d["targets"].empty:
        rep.add(
            table_html(
                d["targets"],
                glossary=False,
                caption="Target BEDs: size and the part inside the common region (fraction_in_common).",
            )
        )
    if not d["callable"].empty:
        rep.add(
            table_html(
                d["callable"],
                glossary=False,
                caption="Callable region per sample: bases of the common region with depth >= the threshold in all alignments of the sample.",
            )
        )


def section_variants(rep, d):
    """Results: input variant sets and their filtering."""
    rep.h(3, "inputs", "Input variant sets and their filtering")
    acc = d["accounting"]
    rep.p(
        "Each input VCF is processed identically: only PASS calls are compared, records on "
        "contigs outside the compared genome and gVCF reference blocks are removed, the "
        "records are split into one ALT allele each, left-aligned and deduplicated, VAF, "
        "depth and alt reads are computed in the same way for every caller, and finally the "
        "uniform filters and the comparison region are applied. The table and the plot show "
        "how many records each step removed and why."
    )
    if acc.empty:
        rep.p("<span class='na'>No accounting available.</span>")
        return
    inputs = acc[acc["step"] == 1][
        ["sample", "technology", "records", "snv", "indel", "other"]
    ].rename(columns={"records": "input_records"})
    kept = acc[(acc["step"] == 7) & (acc["outcome"] == "kept")][
        ["sample", "technology", "records", "snv", "indel"]
    ].rename(columns={"records": "compared", "snv": "compared_snv", "indel": "compared_indel"})
    overview = inputs.merge(kept, on=["sample", "technology"], how="left")
    rep.add(
        table_html(
            overview,
            glossary=False,
            caption="Input records per VCF (SNV, indel, other = symbolic/mixed) and the records compared after all steps.",
        )
    )
    rep.fig(
        plot_filtering(acc),
        "<b>Fate of the input records.</b> Records compared and removed per reason, per VCF.",
    )
    reasons = (
        acc[acc["outcome"].isin(["removed", "added"])]
        .pivot_table(
            index=["sample", "technology"],
            columns="reason",
            values="records",
            aggfunc="sum",
            fill_value=0,
        )
        .rename_axis(columns=None)
        .reset_index()
    )
    rep.add(
        table_html(
            reasons, glossary=False, caption="Records removed (and added by splitting) per reason."
        )
    )
    items = "".join(f"<dt>{esc(k)}</dt><dd>{esc(v)}</dd>" for k, v in REASON_TEXT.items())
    rep.add(f"<details><summary>Reasons</summary><dl class='glossary'>{items}</dl></details>")


def section_metrics(rep, d):
    """Results: comparison metrics."""
    rep.h(3, "metrics", "Comparison metrics")
    rep.p(
        "Recall and precision per sample, and over all samples by variant type and by the "
        "VAF of the truth call. Low recall in the low-VAF bins shows the detection limit of "
        "the query technology. Power-adjusted recall counts only truth variants the query "
        "could detect given its depth (binomial model)."
    )
    rep.add("<div class='figs'>")
    rep.fig(
        plot_metrics_by_type(d["metrics"]),
        "<b>By variant type</b> (SNV includes atomized MNVs; INDEL insertions, deletions and complex variants).",
    )
    rep.fig(plot_recall_by_vaf(d["metrics"]), "<b>Recall by truth VAF</b> (bins of --vaf_bins).")
    rep.add("</div>")
    pp = overall_rows(d["per_pair"])
    if not pp.empty:
        rep.add(
            table_html(
                pp.drop(columns=["stratum_type", "stratum"]),
                caption="Metrics per sample (overall).",
            )
        )
    if not d["power"].empty:
        rep.add(
            table_html(
                d["power"],
                glossary=False,
                caption="Power-adjusted recall (only truth variants with power >= --detection_power).",
            )
        )


def section_truth(rep, d):
    """Results: truth variants."""
    rep.h(3, "truth", "Truth variants")
    rep.p(
        "Every truth variant that was compared, with whether each query technology found it "
        "and, if not, why (class). Variants found by no query technology deserve a look: they "
        "may be artefacts of the truth technology."
    )
    truth = d["truth"]
    if truth.empty:
        rep.p("<span class='na'>No truth variants.</span>")
        return
    order = truth.sort_values(["found_by", "sample", "chrom", "pos"])
    rep.add(
        table_html(
            order,
            max_rows=50,
            caption="Truth variants, those found by the fewest query technologies first.",
        )
    )


def section_discordance(rep, d):
    """Results: discordant variants per class."""
    rep.h(3, "discordance", "Discordant variants")
    rep.p(
        "A discordant variant is a truth variant missed by the query (FN) or a query call "
        "without truth variant (FP). Each is explained from the point of view of the "
        "technology that misses it: for an FN the query, for an FP the truth technology. The "
        "classes are tested in this order; the first that applies is reported."
    )
    variants = d["variants"]
    if variants.empty:
        return
    disc = variants[variants["status"].isin(["FN", "FP"])]
    counts = disc.groupby(["class", "status"]).size().unstack(fill_value=0)
    for cls in CLASS_ORDER:
        sub = disc[disc["class"] == cls]
        n_fn = int(counts.loc[cls, "FN"]) if cls in counts.index and "FN" in counts.columns else 0
        n_fp = int(counts.loc[cls, "FP"]) if cls in counts.index and "FP" in counts.columns else 0
        rep.h(4, f"class-{cls}", f"{cls}: {n_fn} FN, {n_fp} FP")
        rep.p(esc(CLASS_TEXT[cls]))
        if not sub.empty:
            cols = [
                c
                for c in (
                    "sample",
                    "truth",
                    "query",
                    "status",
                    "chrom",
                    "pos",
                    "ref",
                    "alt",
                    "type",
                    "truth_vaf",
                    "query_vaf",
                    "other_filter",
                    "other_tumor_dp",
                    "other_tumor_alt",
                    "normal_vaf",
                    "gene",
                    "consequence",
                )
                if c in sub.columns
            ]
            rep.add(table_html(sub[cols], max_rows=20, glossary=False))


def section_annotation(rep, d):
    """Results: annotation of truth, query and discordant variants."""
    rep.h(3, "annotation", "Annotation")
    rep.p(
        "Consequence (Ensembl VEP), impact, variant class and repeat context per technology "
        "and variant set: truth variants (all compared truth variants), query calls, missed "
        "(FN) and extra (FP) variants. If discordant variants differ from the others (e.g. "
        "more indels or repeats among the FN), the enrichment tests below quantify it."
    )
    counts = d["annotation_counts"]
    if counts.empty:
        rep.p("<span class='na'>No annotation.</span>")
        return
    rep.add("<div class='figs'>")
    for annotation in ("consequence_group", "impact", "type", "in_repeat"):
        rep.fig(
            plot_annotation(counts, annotation),
            f"<b>{esc(annotation)}</b> per variant set and technology (all samples).",
        )
    rep.add("</div>")
    wide = (
        counts.pivot_table(
            index=["sample", "technology", "set", "annotation"],
            columns="value",
            values="count",
            fill_value=0,
        )
        .rename_axis(columns=None)
        .reset_index()
    )
    rep.add(
        f"<details><summary>Annotation counts per sample</summary>{table_html(wide, glossary=False)}</details>"
    )
    if not d["enrichment"].empty:
        rep.add(
            table_html(
                d["enrichment"],
                max_rows=30,
                glossary=False,
                caption="Enrichment of discordant vs. concordant variants per feature (Fisher's exact test, Benjamini-Hochberg q).",
            )
        )
    if not d["genes"].empty:
        rep.add(
            table_html(
                d["genes"].sort_values(
                    "fn" if "fn" in d["genes"] else d["genes"].columns[0], ascending=False
                ),
                max_rows=20,
                glossary=False,
                caption="Genes with discordant variants.",
            )
        )


def section_known(rep, d):
    """Results: known variants."""
    if d["required_summary"].empty and d["known_benchmark"].empty:
        return
    rep.h(3, "known", "Known variants")
    if not d["required_summary"].empty:
        rep.p(
            "Required variants must be found by every technology: PASS (found and compared) or FAIL (filtered, outside the region or not called)."
        )
        rep.add(table_html(d["required_summary"], glossary=False))
        if not d["required_variants"].empty:
            rep.add(table_html(d["required_variants"], max_rows=50, glossary=False))
    if not d["known_benchmark"].empty:
        rep.p("Every technology compared with an external truth set within its confident regions.")
        rep.add(table_html(overall_rows(d["known_benchmark"]), glossary=False))


def section_methods(rep, d, args):
    """Methods: steps and tools, parameters, coverage, checks."""
    rep.h(2, "methods", "Methods")
    versions = d["versions"]
    tools = {}
    for tool_versions in versions.values():
        for tool, version in tool_versions.items():
            tools.setdefault(tool, set()).add(str(version))
    tool_text = ", ".join(f"{t} {'/'.join(sorted(v))}" for t, v in sorted(tools.items()))
    rep.p(
        f"The analysis was run with core-unit-bioinformatics/workflow-nxf-somatic-variantcomp "
        f"{esc(args.pipeline_version)} (Nextflow). Every input was first checked against the "
        "reference genome (contig names and lengths, REF alleles). The VCFs were normalized with "
        "bcftools (sample selection, PASS split, contig renaming, cleaning, norm: split, "
        "left-align, atomize, deduplicate) and harmonized with harmonize_vcf.py (VAF, depth and "
        "alt reads from FORMAT/AF or AD). The comparison region was built with bedtools from the "
        "target BEDs; callable bases came from mosdepth. Variants were matched with RTG vcfeval "
        "(--squash-ploidy, --decompose), exact and position matches with bcftools isec. "
        "Discordant variants were classified with read counts of bcftools mpileup and annotated "
        "with Ensembl VEP and msisensor-pro. "
        f"Tools: {esc(tool_text) or 'see the software versions'}."
    )
    rep.h(3, "parameters", "Parameters")
    rep.add(
        f"<details><summary>Parameters of the run</summary>{table_html(d['params'], glossary=False)}</details>"
    )
    rep.h(3, "coverage", "Coverage")
    rep.p(
        "Mean depth and percent of the bases with at least each depth (--coverage_thresholds), in the common region and in the target region of the technology (genome for WGS; mosdepth's global distribution, 1 % resolution)."
    )
    rep.fig(
        plot_coverage(d["coverage"]),
        "<b>Mean depth per alignment</b> in the common region and in its target region.",
    )
    rep.add(table_html(d["coverage"], glossary=False))
    rep.h(3, "checks", "Consistency checks")
    rep.p(
        "Checks that the counts of the steps add up (input = removed + kept per step; metrics = variant tables; classes = discordant variants; aggregated = sum over the samples). A failed check stops the pipeline after this report was written."
    )
    rep.add(table_html(d["checks"], glossary=False))
    rep.h(3, "links", "Further reports")
    rep.p(
        f"<a href='{esc(args.datavzrd_link)}'>Table browser (datavzrd)</a>: all variants, truth "
        "variants, every discordance class, removed records, regions and coverage as linked, "
        f"searchable tables. <a href='{esc(args.multiqc_link)}'>MultiQC report</a>: quality "
        "control of every step (bcftools norm, harmonization, mosdepth, VEP, input checks)."
    )


def html_report(args):
    """Write the narrative HTML report."""
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
        "accounting": read_many(args.accounting),
        "steps": read(args.region_steps),
        "targets": read(args.targets),
        "callable": read(args.callable),
        "coverage": read(args.coverage),
        "required_summary": read(args.required_summary),
        "required_variants": read(args.required_variants),
        "known_benchmark": read(args.known_benchmark),
        "checks": read(args.checks),
        "samplesheet": read(args.samplesheet, sep=","),
        "params": read(args.params),
        "versions": read_versions(args.versions),
    }
    rep = Report()
    rep.add("<h1>Comparison of somatic variant calls between sequencing technologies</h1>")
    rep.add(
        f"<p class='meta'>workflow-nxf-somatic-variantcomp {esc(args.pipeline_version)} &middot; run {esc(args.run_name)} &middot; {esc(args.start)} &middot; report written {datetime.now().isoformat(timespec='seconds')}</p>"
    )
    failed = d["checks"][d["checks"]["status"] != "ok"] if not d["checks"].empty else d["checks"]
    if d["checks"].empty:
        rep.add("<div class='banner fail'>No consistency checks available.</div>")
    elif failed.empty:
        rep.add(f"<div class='banner ok'>All {len(d['checks'])} consistency checks passed.</div>")
    else:
        rep.add(
            f"<div class='banner fail'>{len(failed)} of {len(d['checks'])} consistency checks FAILED: the results may be wrong. See Methods &rarr; Consistency checks.</div>"
        )
    section_summary(rep, d)
    section_introduction(rep, d)
    rep.h(2, "results", "Results")
    section_regions(rep, d)
    section_variants(rep, d)
    section_metrics(rep, d)
    section_truth(rep, d)
    section_discordance(rep, d)
    section_annotation(rep, d)
    section_known(rep, d)
    section_methods(rep, d, args)
    with open(args.html, "w") as fh:
        fh.write(rep.render("variantcomp report"))
    print(f"INFO: report written: {args.html}", file=sys.stderr)


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
        "checks",
        "samplesheet",
        "params",
        "versions",
    ):
        parser.add_argument(f"--{name}")
    parser.add_argument("--accounting", nargs="*", default=[])
    parser.add_argument("--removed", nargs="*", default=[])
    parser.add_argument("--html", help="output: narrative report (HTML)")
    parser.add_argument("--datavzrd", help="output: folder of the datavzrd tables and config")
    parser.add_argument("--pipeline-version", default="unknown")
    parser.add_argument("--run-name", default="unknown")
    parser.add_argument("--start", default="unknown")
    parser.add_argument("--datavzrd-link", default="tables/index.html")
    parser.add_argument("--multiqc-link", default="../multiqc/multiqc_report.html")
    args = parser.parse_args()
    if not args.html and not args.datavzrd:
        sys.exit("ERROR: nothing to write: give --html and/or --datavzrd")
    if args.datavzrd:
        datavzrd(args)
    if args.html:
        html_report(args)


if __name__ == "__main__":
    main()
