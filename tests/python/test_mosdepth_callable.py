"""Unit tests for bin/mosdepth_callable.py (callable intervals from mosdepth output)."""

import gzip
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "mosdepth_callable.py"


def run(tmp_path, mode, quantized, regions, rename="", extra=()):
    """Write the inputs, run the script and return (process, output lines)."""
    q = tmp_path / "x.quantized.bed.gz"
    r = tmp_path / "x.regions.bed.gz"
    with gzip.open(q, "wt") as fh:
        fh.write(quantized)
    with gzip.open(r, "wt") as fh:
        fh.write(regions)
    (tmp_path / "rename.tsv").write_text(rename)
    (tmp_path / "ref.fa.fai").write_text("chr1\t1000\t6\t60\t61\nchr2\t1000\t1100\t60\t61\n")
    out = tmp_path / "callable.bed"
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--mode",
            mode,
            "--quantized",
            str(q),
            "--regions",
            str(r),
            "--rename",
            str(tmp_path / "rename.tsv"),
            "--fai",
            str(tmp_path / "ref.fa.fai"),
            "--output",
            str(out),
            *extra,
        ],
        capture_output=True,
        text=True,
    )
    lines = out.read_text().splitlines() if out.exists() else []
    return proc, lines


QUANTIZED = (
    "1\t0\t100\tLOW_COVERAGE\n1\t100\t400\tCALLABLE\n2\t0\t50\tCALLABLE\nGL1\t0\t9\tCALLABLE\n"
)
REGIONS = "1\t0\t200\t12.5\n1\t300\t400\t3.0\n"


def test_base_mode_keeps_callable_label_renames_and_drops_unknown_contigs(tmp_path):
    proc, lines = run(tmp_path, "base", QUANTIZED, REGIONS, rename="1\tchr1\n2\tchr2\n")
    assert proc.returncode == 0, proc.stderr
    assert lines == ["chr1\t100\t400", "chr2\t0\t50"]
    assert "1 callable intervals on contigs not in the reference dropped" in proc.stderr


def test_region_mode_uses_mean_depth_threshold(tmp_path):
    proc, lines = run(
        tmp_path, "region", QUANTIZED, REGIONS, rename="1\tchr1\n", extra=["--min-coverage", "10"]
    )
    assert proc.returncode == 0, proc.stderr
    assert lines == ["chr1\t0\t200"]


def test_region_mode_requires_min_coverage(tmp_path):
    proc, _ = run(tmp_path, "region", QUANTIZED, REGIONS)
    assert proc.returncode != 0
    assert "--min-coverage is required" in proc.stderr


def test_unexpected_quantize_label_is_an_error(tmp_path):
    proc, _ = run(tmp_path, "base", "chr1\t0\t10\tMOSDEPTH_Q1\n", REGIONS)
    assert proc.returncode == 1
    assert "unexpected quantize label(s) MOSDEPTH_Q1" in proc.stderr


def test_non_numeric_depth_is_an_error(tmp_path):
    proc, _ = run(
        tmp_path, "region", QUANTIZED, "chr1\t0\t10\tCALLABLE\n", extra=["--min-coverage", "10"]
    )
    assert proc.returncode == 1
    assert "is not a number" in proc.stderr


def test_empty_mosdepth_output_is_an_error(tmp_path):
    proc, _ = run(tmp_path, "base", "", REGIONS)
    assert proc.returncode == 1
    assert "no intervals" in proc.stderr
