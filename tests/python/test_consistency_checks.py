"""Unit tests for bin/consistency_checks.py."""

import csv
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "consistency_checks.py"


def tsv(path, header, rows):
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))
    return str(path)


VARIANT_HEADER = [
    "sample",
    "truth",
    "query",
    "status",
    "class",
    "reasons",
    "type",
    "truth_vaf",
    "query_vaf",
    "other_in_normalized",
    "other_tumor_dp",
    "other_tumor_alt",
    "normal_alt",
    "normal_vaf",
    "in_repeat",
    "giab_homopolymer",
]
TP_ROW = [
    "S1",
    "wgs",
    "wes",
    "TP",
    "concordant",
    "concordant",
    "SNV",
    0.3,
    0.3,
    "yes",
    60,
    20,
    0,
    0,
    "no",
    "no",
]
FN_ROW = [
    "S1",
    "wgs",
    "wes",
    "FN",
    "low_coverage",
    "low_coverage,in_repeat",
    "SNV",
    0.3,
    "NA",
    "no",
    5,
    0,
    0,
    0,
    "yes",
    "no",
]
FP_ROW = [
    "S1",
    "wgs",
    "wes",
    "FP",
    "no_evidence",
    "no_evidence,indel,low_vaf,giab_homopolymer",
    "DEL",
    "NA",
    0.05,
    "no",
    50,
    0,
    0,
    0,
    "no",
    "yes",
]


def inputs(tmp_path, fn=1, fn_row=None, tp_alt=20):
    """A consistent run of one pair (S1 wgs vs wes); fn sets FN of the metrics, fn_row
    replaces the FN variant row, tp_alt the alt reads of the TP rows."""
    metric = [
        "sample",
        "truth",
        "query",
        "stratum_type",
        "stratum",
        "tp_baseline",
        "fn",
        "tp_call",
        "fp",
    ]
    per_pair = tsv(
        tmp_path / "pp.tsv", metric, [["S1", "wgs", "wes", "overall", "all", 9, fn, 9, 1]]
    )
    aggregated = tsv(
        tmp_path / "agg.tsv",
        ["truth", "query", "stratum_type", "stratum", "tp_baseline", "fn", "tp_call", "fp"],
        [["wgs", "wes", "overall", "all", 9, fn, 9, 1]],
    )
    tp_row = TP_ROW[:11] + [tp_alt] + TP_ROW[12:]
    variants = tsv(
        tmp_path / "v.tsv",
        VARIANT_HEADER,
        [tp_row] * 9 + [fn_row or FN_ROW, FP_ROW],
    )
    classes = tsv(
        tmp_path / "c.tsv",
        ["truth", "query", "status", "low_coverage", "no_evidence"],
        [["wgs", "wes", "FN", 1, 0], ["wgs", "wes", "FP", 0, 1]],
    )
    steps = tsv(
        tmp_path / "steps.tsv",
        ["step", "name", "intervals", "bp"],
        [[1, "genome", 1, 1000], [2, "targets_intersection", 1, 500]],
    )
    targets = tsv(
        tmp_path / "t.tsv", ["region", "type", "intervals", "bp"], [["common", "common", 1, 500]]
    )
    accounting = tsv(
        tmp_path / "S1_wgs.accounting.tsv",
        ["sample", "technology", "step", "step_name", "outcome", "reason", "records"],
        [["S1", "wgs", 7, "filters", "kept", ".", 10]],
    )
    return [
        "--per-pair",
        per_pair,
        "--aggregated",
        aggregated,
        "--variants",
        variants,
        "--classes",
        classes,
        "--region-steps",
        steps,
        "--targets",
        targets,
        "--accounting",
        accounting,
    ]


def run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_consistent_run_passes_and_gate_is_silent(tmp_path):
    out = tmp_path / "checks.tsv"
    result = run("check", *inputs(tmp_path), "--output", str(out))
    assert result.returncode == 0, result.stderr
    checks = rows(out)
    assert checks and all(c["status"] == "ok" for c in checks), checks
    # TP_baseline + FN (10) >= compared truth records (10)
    assert any(c["check"].startswith("TP_baseline + FN") for c in checks)
    assert run("gate", "--checks", str(out)).returncode == 0


def test_contradictory_reasons_fail(tmp_path):
    # called in the query (non-PASS) but also "alt reads, not called": exclusive reasons
    bad = FN_ROW[:4] + ["filtered_non_pass", "filtered_non_pass,subthreshold_evidence,in_repeat"]
    bad += ["SNV", 0.3, "NA", "yes", 50, 3, 0, 0, "yes", "no"]
    out = tmp_path / "checks.tsv"
    assert run("check", *inputs(tmp_path, fn_row=bad), "--output", str(out)).returncode == 0
    failed = [c for c in rows(out) if c["status"] != "ok"]
    assert [c["category"] for c in failed] == ["reasons"], failed
    assert "main reason: 1" in failed[0]["check"]


def test_context_reason_missing_fails(tmp_path):
    # in a repeat, but the reason in_repeat is missing
    bad = FN_ROW[:5] + ["low_coverage"] + FN_ROW[6:]
    out = tmp_path / "checks.tsv"
    assert run("check", *inputs(tmp_path, fn_row=bad), "--output", str(out)).returncode == 0
    failed = [c for c in rows(out) if c["status"] != "ok"]
    assert len(failed) == 1 and "in_repeat: 1" in failed[0]["check"]


def test_tp_without_alt_reads_fails_the_evidence_check(tmp_path):
    out = tmp_path / "checks.tsv"
    assert run("check", *inputs(tmp_path, tp_alt=0), "--output", str(out)).returncode == 0
    failed = [c for c in rows(out) if c["status"] != "ok"]
    assert [c["category"] for c in failed] == ["evidence"]
    assert failed[0]["observed"] == "0.000"


def test_inconsistent_counts_fail_and_gate_stops(tmp_path):
    out = tmp_path / "checks.tsv"
    # metrics say 2 FN, but the variant table and the classes have 1
    result = run("check", *inputs(tmp_path, fn=2), "--output", str(out))
    assert result.returncode == 0  # check never stops
    failed = {c["check"] for c in rows(out) if c["status"] == "FAIL"}
    assert "FN rows of the variant table = FN of the metrics" in failed
    assert "discordance classes add up to the variants" in failed
    gate = run("gate", "--checks", str(out))
    assert gate.returncode == 1 and "consistency checks failed" in gate.stderr
