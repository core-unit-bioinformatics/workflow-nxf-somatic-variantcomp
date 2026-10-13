#!/usr/bin/env python3
"""Statistics, plots and MultiQC files of a variantcomp run.

Inputs (from the benchmark and discordance steps):
  --metrics     metrics.aggregated.tsv
  --per-pair    metrics.per_pair.tsv
  --power       power.aggregated.tsv (optional)
  --variants    variants.annotated.tsv (every compared variant, TP/FN/FP)

Outputs (prefix "variantcomp"):
  <prefix>.enrichment.tsv    per truth/query comparison: are discordant variants
                             (FN, FP, or both) enriched for a feature compared with
                             the concordant variants (TP)? Features: indel, low VAF
                             (< --low-vaf; truth VAF for FN/TP, query VAF for FP),
                             and every stratification column (yes/no).
                             Fisher's exact test, odds ratio, Benjamini-Hochberg q.
  <prefix>.genes.tsv         per comparison and gene: TP/FN/FP counts, discordance
                             rate, Fisher's exact test of discordant vs. concordant
                             in the gene vs. all other genes, BH q
  <prefix>.classes.tsv       per comparison: counts of discordance classes
  *_mqc.*                    MultiQC custom content (tables, bar and line graphs)
  plots/*.png                recall/precision, recall by VAF, classes, enrichment

Checks: the required columns of every input, status values TP/FN/FP and known
discordance classes (an unknown class would be missing from the class summary).

Why a script: the statistics (Fisher's exact test, Benjamini-Hochberg), plots and
MultiQC custom-content files are analysis code; no tool produces them from these
tables.
"""

import argparse
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

BASE_COLUMNS = {
    "sample",
    "truth",
    "query",
    "status",
    "chrom",
    "pos",
    "ref",
    "alt",
    "type",
    "length",
    "truth_vaf",
    "truth_dp",
    "truth_alt",
    "query_vaf",
    "query_dp",
    "query_alt",
    "query_filter",
    "class",
    "other_in_normalized",
    "other_filter",
    "other_vaf_normalized",
    "other_tumor_dp",
    "other_tumor_alt",
    "normal_dp",
    "normal_alt",
    "normal_vaf",
    "query_power",
    "detectable",
    "gene",
}
CLASS_ORDER = [
    "filtered_non_pass",
    "present_but_filtered",
    "germline_in_normal",
    "low_coverage",
    "subthreshold_evidence",
    "no_evidence",
    "absent",
]
METRICS_COLUMNS = [
    "truth",
    "query",
    "stratum_type",
    "stratum",
    "n_samples",
    "tp_baseline",
    "fn",
    "tp_call",
    "fp",
    "micro_recall",
    "micro_precision",
    "micro_f1",
    "macro_recall_mean",
    "macro_precision_mean",
]
PER_PAIR_COLUMNS = [
    "sample",
    "truth",
    "query",
    "stratum_type",
    "stratum",
    "tp_baseline",
    "fn",
    "tp_call",
    "fp",
    "recall",
    "precision",
    "f1",
]
VARIANT_COLUMNS = [
    "sample",
    "truth",
    "query",
    "status",
    "pos",
    "type",
    "truth_vaf",
    "query_vaf",
    "class",
]
COLORS = {"recall": "#2a6f97", "precision": "#c8553d", "f1": "#6a994e"}


def require_columns(df, path, columns):
    """Stop if a table misses one of the columns."""
    missing = [c for c in columns if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: {path}: missing column(s) {', '.join(missing)}")


def check_variants(variants, path):
    """Required columns and known status/class values of the variant table."""
    require_columns(variants, path, VARIANT_COLUMNS)
    # status and class values must be known, otherwise counts would be lost silently
    bad_status = sorted(set(variants["status"]) - {"TP", "FN", "FP"})
    if bad_status:
        sys.exit(f"ERROR: {path}: unknown status value(s) {', '.join(bad_status)}")
    bad_class = sorted(set(variants["class"]) - set(CLASS_ORDER) - {"concordant"})
    if bad_class:
        sys.exit(
            f"ERROR: {path}: unknown discordance class(es) {', '.join(bad_class)} "
            f"(known: concordant, {', '.join(CLASS_ORDER)})"
        )


def bh(pvalues):
    """Benjamini-Hochberg adjusted p-values (NaN kept, not counted as tests)."""
    p = np.asarray(pvalues, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if not ok.any():
        return q
    ps = p[ok]
    # p-values in ascending order, scaled by n / rank
    order = np.argsort(ps)
    ranked = ps[order] * len(ps) / (np.arange(len(ps)) + 1)
    # enforce monotonicity from the largest p-value downwards
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ps)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def fisher(a, b, c, d):
    """Fisher's exact test of [[a, b], [c, d]]; (odds ratio, p).

    NaN if a row or column of the table is empty (the test is not defined).
    """
    if min(a + b, c + d, a + c, b + d) == 0:
        return math.nan, math.nan
    odds, p = fisher_exact([[a, b], [c, d]])
    return odds, p


def feature_flags(df, low_vaf):
    """Per-variant boolean features; VAF from truth (TP, FN) or query (FP)."""
    vaf = np.where(df["status"] == "FP", df["query_vaf"], df["truth_vaf"])
    vaf = pd.to_numeric(pd.Series(vaf, index=df.index), errors="coerce")
    features = {
        "indel": df["type"].isin(["INS", "DEL", "COMPLEX"]),
        f"vaf_lt_{low_vaf:g}": vaf < low_vaf,
    }
    # every other yes/no column is a stratification (from annotate_variants.py)
    for col in df.columns:
        if col not in BASE_COLUMNS and set(df[col].dropna().unique()) <= {"yes", "no"}:
            features[col] = df[col] == "yes"
    return features


def enrichment(df, low_vaf):
    """Fisher's exact test per comparison, discordance status and feature (BH over all tests)."""
    rows = []
    for (truth, query), group in df.groupby(["truth", "query"], sort=True):
        features = feature_flags(group, low_vaf)
        concordant = group["status"] == "TP"
        for label, discordant in (
            ("FN", group["status"] == "FN"),
            ("FP", group["status"] == "FP"),
            ("FN+FP", group["status"] != "TP"),
        ):
            for name, flag in features.items():
                # 2x2 table: discordant / concordant x with / without the feature
                a = int((discordant & flag).sum())
                b = int((discordant & ~flag).sum())
                c = int((concordant & flag).sum())
                d = int((concordant & ~flag).sum())
                odds, p = fisher(a, b, c, d)
                rows.append(
                    {
                        "truth": truth,
                        "query": query,
                        "discordant": label,
                        "feature": name,
                        "discordant_with": a,
                        "discordant_without": b,
                        "concordant_with": c,
                        "concordant_without": d,
                        "fraction_discordant": a / (a + b) if a + b else math.nan,
                        "fraction_concordant": c / (c + d) if c + d else math.nan,
                        "odds_ratio": odds,
                        "p_value": p,
                    }
                )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_value"] = bh(out["p_value"])
    return out


def gene_summary(df):
    """Per comparison and gene: counts and Fisher's test of the gene vs. all other genes."""
    if "gene" not in df.columns:
        return pd.DataFrame()
    # a variant overlapping several genes counts for each; '.' = no gene
    exploded = df.assign(gene=df["gene"].fillna(".").str.split(",")).explode("gene")
    exploded = exploded[exploded["gene"] != "."]
    rows = []
    for (truth, query), group in exploded.groupby(["truth", "query"], sort=True):
        disc_total = int((group["status"] != "TP").sum())
        conc_total = int((group["status"] == "TP").sum())
        for gene, g in group.groupby("gene", sort=True):
            tp = int((g["status"] == "TP").sum())
            fn = int((g["status"] == "FN").sum())
            fp = int((g["status"] == "FP").sum())
            odds, p = fisher(fn + fp, disc_total - fn - fp, tp, conc_total - tp)
            rows.append(
                {
                    "truth": truth,
                    "query": query,
                    "gene": gene,
                    "tp": tp,
                    "fn": fn,
                    "fp": fp,
                    "discordance_rate": (fn + fp) / (tp + fn + fp),
                    "samples": g["sample"].nunique(),
                    "odds_ratio": odds,
                    "p_value": p,
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_value"] = bh(out["p_value"])
        out = out.sort_values(["truth", "query", "p_value", "gene"], na_position="last")
    return out


def class_summary(df):
    """Counts of the discordance classes per comparison and status (FN, FP)."""
    disc = df[df["status"] != "TP"]
    if disc.empty:
        return pd.DataFrame(columns=["truth", "query", "status", *CLASS_ORDER])
    table = disc.pivot_table(
        index=["truth", "query", "status"],
        columns="class",
        values="pos",
        aggfunc="count",
        fill_value=0,
    )
    for c in CLASS_ORDER:
        if c not in table.columns:
            table[c] = 0
    return table[CLASS_ORDER].reset_index()


# ------------------------------------------------------------------ MultiQC files


def mqc_header(section_id, name, description, plot_type, extra=""):
    return (
        f"# id: 'variantcomp_{section_id}'\n# section_name: '{name}'\n"
        f"# description: '{description}'\n# plot_type: '{plot_type}'\n{extra}"
    )


def write_mqc(prefix, metrics, per_pair, power, classes):
    """MultiQC custom-content files: metrics tables, class bar graph, recall-by-VAF line graph."""
    # overall metrics per comparison, with the power-adjusted recall if available
    overall = metrics[(metrics["stratum_type"] == "overall")].copy()
    if power is not None and not power.empty:
        overall = overall.merge(
            power[["truth", "query", "micro_recall"]].rename(
                columns={"micro_recall": "power_adjusted_recall"}
            ),
            on=["truth", "query"],
            how="left",
        )
    overall.insert(0, "comparison", overall["truth"] + " vs " + overall["query"])
    for col in ("n_samples", "tp_baseline", "fn", "tp_call", "fp"):
        overall[col] = pd.to_numeric(overall[col], errors="coerce").astype("Int64")
    cols = [
        "comparison",
        "n_samples",
        "tp_baseline",
        "fn",
        "tp_call",
        "fp",
        "micro_recall",
        "micro_precision",
        "micro_f1",
        "macro_recall_mean",
        "macro_precision_mean",
    ]
    if "power_adjusted_recall" in overall.columns:
        cols.append("power_adjusted_recall")
    with open(f"{prefix}_metrics_mqc.tsv", "w") as fh:
        fh.write(
            mqc_header(
                "metrics",
                "Recall and precision",
                "Per comparison, summed over samples (micro) and mean over samples " "(macro).",
                "table",
            )
        )
        overall[cols].to_csv(fh, sep="\t", index=False, na_rep="NA", float_format="%.4f")

    pairs = per_pair[per_pair["stratum_type"] == "overall"].copy()
    pairs.insert(0, "comparison", pairs["sample"] + ": " + pairs["truth"] + " vs " + pairs["query"])
    with open(f"{prefix}_per_sample_mqc.tsv", "w") as fh:
        fh.write(
            mqc_header(
                "per_sample",
                "Recall and precision per sample",
                "Every sample and technology pair.",
                "table",
            )
        )
        pairs[
            ["comparison", "tp_baseline", "fn", "tp_call", "fp", "recall", "precision", "f1"]
        ].to_csv(fh, sep="\t", index=False)

    if not classes.empty:
        cl = classes.copy()
        cl.insert(0, "comparison", cl["truth"] + " vs " + cl["query"] + " " + cl["status"])
        with open(f"{prefix}_classes_mqc.tsv", "w") as fh:
            fh.write(
                mqc_header(
                    "classes",
                    "Discordance classes",
                    "Why variants were missed (FN) or only called by the query " "(FP).",
                    "bargraph",
                )
            )
            cl[["comparison", *CLASS_ORDER]].to_csv(fh, sep="\t", index=False)

    vaf = metrics[metrics["stratum_type"] == "truth_vaf"]
    lines = [
        "id: 'variantcomp_recall_vaf'",
        "section_name: 'Recall by truth VAF'",
        "description: 'Micro recall per truth VAF bin.'",
        "plot_type: 'linegraph'",
        "pconfig:",
        "  id: 'variantcomp_recall_vaf_plot'",
        "  title: 'Recall by truth VAF'",
        "  xlab: 'VAF bin (lower edge)'",
        "  ylab: 'Recall'",
        "  ymin: 0",
        "  ymax: 1",
        "data:",
    ]
    for (truth, query), group in vaf.groupby(["truth", "query"]):
        lines.append(f"  '{truth} vs {query}':")
        for _, row in group.iterrows():
            if row["stratum"] == "unknown" or pd.isna(row["micro_recall"]):
                continue
            lower = float(row["stratum"].strip("[)").split(",")[0])
            lines.append(f"    {lower}: {row['micro_recall']:.4f}")
    Path(f"{prefix}_recall_vaf_mqc.yaml").write_text("\n".join(lines) + "\n")


# ------------------------------------------------------------------ plots


def plot_metrics(per_pair, metrics, outdir):
    overall = metrics[metrics["stratum_type"] == "overall"]
    if overall.empty:
        return
    labels = (overall["truth"] + " vs " + overall["query"]).tolist()
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(max(4, 1.8 * len(labels) + 2), 4))
    width = 0.25
    for i, metric in enumerate(["recall", "precision", "f1"]):
        ax.bar(
            x + (i - 1) * width,
            overall[f"micro_{metric}"].astype(float),
            width,
            label=metric,
            color=COLORS[metric],
            alpha=0.85,
        )
        samples = per_pair[per_pair["stratum_type"] == "overall"]
        for j, (_, row) in enumerate(overall.iterrows()):
            values = pd.to_numeric(
                samples[(samples["truth"] == row["truth"]) & (samples["query"] == row["query"])][
                    metric
                ],
                errors="coerce",
            ).dropna()
            ax.scatter(
                np.full(len(values), j + (i - 1) * width), values, color="black", s=10, zorder=3
            )
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("value (bars: all samples; dots: per sample)")
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout()
    fig.savefig(outdir / "recall_precision.png", dpi=150)
    plt.close(fig)


def plot_vaf(metrics, outdir):
    vaf = metrics[(metrics["stratum_type"] == "truth_vaf") & (metrics["stratum"] != "unknown")]
    if vaf.empty:
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    for (truth, query), group in vaf.groupby(["truth", "query"]):
        ax.plot(
            group["stratum"],
            group["micro_recall"].astype(float),
            marker="o",
            label=f"{truth} vs {query}",
        )
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("truth VAF")
    ax.set_ylabel("recall")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "recall_by_vaf.png", dpi=150)
    plt.close(fig)


def plot_classes(classes, outdir):
    if classes.empty:
        return
    labels = (classes["truth"] + " vs " + classes["query"] + " " + classes["status"]).tolist()
    fig, ax = plt.subplots(figsize=(9, max(2.5, 0.5 * len(labels) + 1)))
    left = np.zeros(len(labels))
    cmap = plt.get_cmap("tab10")
    for i, c in enumerate(CLASS_ORDER):
        values = classes[c].to_numpy(dtype=float)
        ax.barh(labels, values, left=left, label=c, color=cmap(i))
        left += values
    ax.set_xlabel("variants")
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    fig.savefig(outdir / "discordance_classes.png", dpi=150)
    plt.close(fig)


def plot_enrichment(enrich, outdir, fdr=0.05):
    data = enrich[(enrich["discordant"] == "FN+FP") & enrich["odds_ratio"].notna()].copy()
    data = data[np.isfinite(data["odds_ratio"]) & (data["odds_ratio"] > 0)]
    if data.empty:
        return
    data["label"] = data["truth"] + " vs " + data["query"] + ": " + data["feature"]
    fig, ax = plt.subplots(figsize=(6, max(2.5, 0.35 * len(data) + 1)))
    colors = ["#c8553d" if q < fdr else "#7f7f7f" for q in data["q_value"].fillna(1)]
    ax.scatter(np.log2(data["odds_ratio"]), data["label"], color=colors)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel(f"log2 odds ratio, discordant vs. concordant (red: q < {fdr:g})")
    fig.tight_layout()
    fig.savefig(outdir / "enrichment.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--per-pair", required=True)
    parser.add_argument("--power")
    parser.add_argument("--variants", required=True)
    parser.add_argument("--low-vaf", type=float, default=0.1)
    parser.add_argument(
        "--fdr", type=float, default=0.05, help="q-value below which a test is significant"
    )
    parser.add_argument("--prefix", default="variantcomp")
    args = parser.parse_args()

    def read(path):
        return pd.read_csv(path, sep="\t", na_values=["NA"], keep_default_na=False)

    # inputs with their required columns
    metrics = read(args.metrics)
    require_columns(metrics, args.metrics, METRICS_COLUMNS)
    per_pair = read(args.per_pair)
    require_columns(per_pair, args.per_pair, PER_PAIR_COLUMNS)
    power = read(args.power) if args.power else None
    if power is not None:
        require_columns(power, args.power, ["truth", "query", "micro_recall"])
    variants = read(args.variants)
    check_variants(variants, args.variants)
    # VAFs as numbers (NA -> NaN)
    for col in ("truth_vaf", "query_vaf"):
        variants[col] = pd.to_numeric(variants[col], errors="coerce")
    print(
        f"INFO: {len(variants)} variants, {variants.groupby(['truth', 'query']).ngroups} "
        "comparisons",
        file=sys.stderr,
    )

    enrich = enrichment(variants, args.low_vaf)
    genes = gene_summary(variants)
    classes = class_summary(variants)
    enrich.to_csv(
        f"{args.prefix}.enrichment.tsv", sep="\t", index=False, na_rep="NA", float_format="%.4g"
    )
    genes.to_csv(
        f"{args.prefix}.genes.tsv", sep="\t", index=False, na_rep="NA", float_format="%.4g"
    )
    classes.to_csv(f"{args.prefix}.classes.tsv", sep="\t", index=False)
    write_mqc(args.prefix, metrics, per_pair, power, classes)

    plots = Path("plots")
    plots.mkdir(exist_ok=True)
    plot_metrics(per_pair, metrics, plots)
    plot_vaf(metrics, plots)
    plot_classes(classes, plots)
    if not enrich.empty:
        plot_enrichment(enrich, plots, args.fdr)
    n_sig = int((enrich["q_value"] < args.fdr).sum()) if not enrich.empty else 0
    print(
        f"INFO: {len(enrich)} enrichment tests ({n_sig} with q < {args.fdr:g}), {len(genes)} gene rows; "
        f"plots in {plots}/",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
