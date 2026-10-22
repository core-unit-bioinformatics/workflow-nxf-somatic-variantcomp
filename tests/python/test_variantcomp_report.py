"""Unit tests for bin/variantcomp_report.py."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

BIN = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN))
import variantcomp_report as report  # noqa: E402


def test_bh():
    q = report.bh([0.01, 0.04, np.nan, 0.03])
    assert q[0] == pytest.approx(0.03)
    assert q[1] == pytest.approx(0.04)
    assert np.isnan(q[2])
    assert q[3] == pytest.approx(0.04)


def variants():
    rows = []
    # 20 concordant SNVs at high VAF, 8 FP at low VAF (6 in a homopolymer), 2 FN indels
    for _ in range(20):
        rows.append(("TP", "SNV", 0.3, 0.3, "no", "GENE_A"))
    for i in range(8):
        rows.append(("FP", "SNV", np.nan, 0.03, "yes" if i < 6 else "no", "GENE_B"))
    for _ in range(2):
        rows.append(("FN", "DEL", 0.3, np.nan, "no", "GENE_A"))
    df = pd.DataFrame(
        rows, columns=["status", "type", "truth_vaf", "query_vaf", "homopolymer", "gene"]
    )
    df["sample"], df["truth"], df["query"] = "S1", "wgs", "wes"
    df["pos"] = range(len(df))
    df["class"] = np.where(df["status"] == "TP", "concordant", "absent")
    return df


def test_enrichment_finds_low_vaf_and_homopolymer():
    out = report.enrichment(variants(), 0.1)
    fp = out[out["discordant"] == "FP"].set_index("feature")
    assert set(fp.index) == {"indel", "vaf_lt_0.1", "homopolymer"}
    assert fp.loc["vaf_lt_0.1", "discordant_with"] == 8
    assert fp.loc["vaf_lt_0.1", "p_value"] < 1e-5
    assert fp.loc["homopolymer", "discordant_with"] == 6
    assert fp.loc["homopolymer", "q_value"] < 0.05


def test_gene_summary():
    genes = report.gene_summary(variants()).set_index("gene")
    assert (genes.loc["GENE_B", "fp"], genes.loc["GENE_B", "tp"]) == (8, 0)
    assert genes.loc["GENE_B", "discordance_rate"] == 1.0
    assert genes.loc["GENE_B", "p_value"] < 0.001


def test_class_summary():
    classes = report.class_summary(variants())
    assert classes.set_index("status").loc["FP", "absent"] == 8


def test_check_variants_accepts_valid_table():
    report.check_variants(variants(), "variants.tsv")


@pytest.mark.parametrize(
    ("column", "value", "message"),
    [
        ("class", "outside_callable", "unknown discordance class(es) outside_callable"),
        ("status", "TN", "unknown status value(s) TN"),
    ],
)
def test_check_variants_rejects_unknown_values(capsys, column, value, message):
    df = variants()
    df.loc[0, column] = value
    with pytest.raises(SystemExit) as err:
        report.check_variants(df, "variants.tsv")
    assert message in str(err.value)


def test_check_variants_missing_column():
    with pytest.raises(SystemExit) as err:
        report.check_variants(variants().drop(columns="class"), "variants.tsv")
    assert "missing column(s) class" in str(err.value)


def annotated():
    df = variants()
    df["consequence_group"] = ["non_synonymous"] * 10 + ["synonymous"] * 10 + ["other"] * 10
    df["impact"] = "MODERATE"
    df["in_repeat"] = np.where(df["status"] == "FP", "yes", "no")
    return df


def test_strata_summary_counts_per_value():
    strata = report.strata_summary(annotated())
    rows = strata[(strata["annotation"] == "in_repeat") & (strata["value"] == "yes")]
    assert (rows["tp"].item(), rows["fn"].item(), rows["fp"].item()) == (0, 0, 8)
    assert rows["discordance_rate"].item() == pytest.approx(1.0)
    assert set(strata["annotation"]) == {"consequence_group", "impact", "in_repeat", "homopolymer"}


def test_feature_flags_include_consequence_and_impact():
    flags = report.feature_flags(annotated(), 0.1)
    assert {"non_synonymous", "synonymous", "high_or_moderate_impact", "in_repeat"} <= set(flags)


def comparison_rows():
    """Two queries (wes, panel) against wgs for one sample: v1 found by both, v2 missed by
    panel, one FP of wes."""
    base = {
        "sample": "S1",
        "truth": "wgs",
        "chrom": "chr1",
        "ref": "A",
        "alt": "C",
        "type": "SNV",
        "consequence_group": "non_synonymous",
        "impact": "MODERATE",
        "variant_class": "SNV",
        "in_repeat": "no",
        "truth_vaf": 0.3,
    }
    rows = [
        {**base, "query": "wes", "pos": 100, "status": "TP", "class": "concordant"},
        {**base, "query": "panel", "pos": 100, "status": "TP", "class": "concordant"},
        {**base, "query": "wes", "pos": 200, "status": "TP", "class": "concordant"},
        {**base, "query": "panel", "pos": 200, "status": "FN", "class": "absent"},
        {
            **base,
            "query": "wes",
            "pos": 300,
            "status": "FP",
            "class": "no_evidence",
            "truth_vaf": None,
        },
    ]
    return pd.DataFrame(rows)


def test_truth_variants_found_per_query():
    table = report.truth_variants(comparison_rows())
    assert list(table["pos"]) == [100, 200]
    v2 = table[table["pos"] == 200].iloc[0]
    assert (
        v2["found_by"],
        v2["queries"],
        v2["found_wes"],
        v2["found_panel"],
        v2["class_panel"],
        v2["class_wes"],
    ) == (1, 2, "yes", "no", "absent", ".")


def test_annotation_counts_per_set():
    counts = report.annotation_counts(comparison_rows())

    def get(tech, set_):
        rows = counts[(counts["technology"] == tech) & (counts["set"] == set_)]
        return int(rows[rows["annotation"] == "type"]["count"].sum())

    # truth variants: each once (2); wes calls: 2 TP + 1 FP; panel misses 1; wes has 1 FP
    assert (
        get("wgs", "truth_variants"),
        get("wes", "query_calls"),
        get("panel", "missed_FN"),
        get("wes", "extra_FP"),
    ) == (2, 3, 1, 1)
