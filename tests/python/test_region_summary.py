"""Unit tests for bin/region_summary.py (target and callable region summaries)."""

import csv
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "region_summary.py"


def write(path, text):
    path.write_text(text)
    return str(path)


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)


def test_targets_summary(tmp_path):
    common = write(tmp_path / "common.bed", "chr1\t100\t200\n")
    wes = write(tmp_path / "wes.bed", "chr1\t0\t200\nchr1\t300\t400\n")
    genome = write(tmp_path / "genome.bed", "chr1\t0\t1000\n")
    out = tmp_path / "targets.tsv"
    proc = run(
        "targets",
        "--common",
        common,
        "--genome",
        genome,
        "--target",
        f"wes={wes}",
        "--genome-wide",
        "wgs",
        "--output",
        str(out),
    )
    assert proc.returncode == 0, proc.stderr
    table = {r["region"]: r for r in rows(out)}
    assert table["wes"]["bp"] == "300" and table["wes"]["bp_in_common"] == "100"
    assert table["wes"]["fraction_in_common"] == "0.3333"
    assert table["wgs"]["type"] == "genome" and table["wgs"]["fraction_in_common"] == "0.1000"
    assert table["common"]["bp"] == "100"


def test_targets_empty_common_region_is_an_error(tmp_path):
    common = write(tmp_path / "common.bed", "")
    genome = write(tmp_path / "genome.bed", "chr1\t0\t1000\n")
    proc = run(
        "targets", "--common", common, "--genome", genome, "--output", str(tmp_path / "t.tsv")
    )
    assert proc.returncode == 1
    assert "common comparison region is empty" in proc.stderr


def test_overlapping_intervals_are_an_error(tmp_path):
    common = write(tmp_path / "common.bed", "chr1\t100\t200\nchr1\t150\t300\n")
    genome = write(tmp_path / "genome.bed", "chr1\t0\t1000\n")
    proc = run(
        "targets", "--common", common, "--genome", genome, "--output", str(tmp_path / "t.tsv")
    )
    assert proc.returncode == 1
    assert "overlap or are unsorted" in proc.stderr


def test_callable_summary(tmp_path):
    common = write(tmp_path / "common.bed", "chr1\t0\t100\nchr2\t0\t100\n")
    callable_bed = write(tmp_path / "callable.bed", "chr1\t10\t60\n")
    out = tmp_path / "callable.tsv"
    proc = run(
        "callable",
        "--sample",
        "S1",
        "--common",
        common,
        "--callable-bed",
        callable_bed,
        "--alignments",
        "4",
        "--mode",
        "base",
        "--min-coverage",
        "10",
        "--output",
        str(out),
    )
    assert proc.returncode == 0, proc.stderr
    assert rows(out) == [
        {
            "sample": "S1",
            "alignments": "4",
            "mode": "base",
            "min_coverage": "10",
            "common_bp": "200",
            "callable_bp": "50",
            "fraction_callable": "0.2500",
        }
    ]


def test_callable_outside_common_region_is_an_error(tmp_path):
    common = write(tmp_path / "common.bed", "chr1\t0\t100\n")
    callable_bed = write(tmp_path / "callable.bed", "chr1\t50\t150\n")
    proc = run(
        "callable",
        "--sample",
        "S1",
        "--common",
        common,
        "--callable-bed",
        callable_bed,
        "--alignments",
        "2",
        "--mode",
        "base",
        "--min-coverage",
        "10",
        "--output",
        str(tmp_path / "c.tsv"),
    )
    assert proc.returncode == 1
    assert "extends beyond the common region" in proc.stderr


def test_targets_below_min_bp_fails(tmp_path):
    common = write(tmp_path / "common.bed", "chr1\t100\t200\n")
    genome = write(tmp_path / "genome.bed", "chr1\t0\t1000\n")
    proc = run(
        "targets",
        "--common",
        common,
        "--genome",
        genome,
        "--output",
        str(tmp_path / "t.tsv"),
        "--min-bp",
        "500",
    )
    assert proc.returncode == 1
    assert "100 bp, less than --min_region_bp 500" in proc.stderr


def test_callable_below_min_bp_fails(tmp_path):
    common = write(tmp_path / "common.bed", "chr1\t0\t100\n")
    callable_bed = write(tmp_path / "callable.bed", "chr1\t10\t60\n")
    proc = run(
        "callable",
        "--sample",
        "S1",
        "--common",
        common,
        "--callable-bed",
        callable_bed,
        "--alignments",
        "4",
        "--mode",
        "base",
        "--min-coverage",
        "10",
        "--output",
        str(tmp_path / "c.tsv"),
        "--min-bp",
        "51",
    )
    assert proc.returncode == 1
    assert "S1: only 50 bp of the common region are callable" in proc.stderr
