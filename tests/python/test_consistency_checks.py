"""Unit tests for bin/consistency_checks.py."""

import csv
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "consistency_checks.py"


def tsv(path, header, rows):
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))
    return str(path)


def inputs(tmp_path, fn=1):
    """A consistent run of one pair (S1 wgs vs wes); fn sets FN of the metrics."""
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
    variants = tsv(
        tmp_path / "v.tsv",
        ["sample", "truth", "query", "status", "class"],
        [["S1", "wgs", "wes", "TP", "."]] * 9
        + [["S1", "wgs", "wes", "FN", "absent"], ["S1", "wgs", "wes", "FP", "no_evidence"]],
    )
    classes = tsv(
        tmp_path / "c.tsv",
        ["truth", "query", "status", "absent", "no_evidence"],
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
