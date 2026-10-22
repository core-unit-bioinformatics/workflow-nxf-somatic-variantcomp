"""Unit tests for bin/coverage_summary.py."""

import csv
import gzip
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "coverage_summary.py"


def gz(path, text):
    with gzip.open(path, "wt") as fh:
        fh.write(text)


def run(tmp_path, thresholds="10,30"):
    out = tmp_path / "cov"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--alignments",
            str(tmp_path / "alignments.tsv"),
            "--dir",
            str(tmp_path),
            "--thresholds",
            thresholds,
            "--prefix",
            str(out),
        ],
        capture_output=True,
        text=True,
    )
    return result, out


def test_common_target_and_genome(tmp_path):
    (tmp_path / "alignments.tsv").write_text(
        "alignment\tsample\ttechnology\trole\nS1_wes_tumor\tS1\twes\ttumor\nS1_wgs_tumor\tS1\twgs\ttumor\n"
    )
    for aid in ("S1_wes_tumor", "S1_wgs_tumor"):
        # common region: 100 bp at depth 20 and 300 bp at depth 40 -> mean 35;
        # >= 10X: 400 bp (100 %), >= 30X: 300 bp (75 %)
        gz(tmp_path / f"{aid}.regions.bed.gz", "chr1\t0\t100\t20\nchr1\t200\t500\t40\n")
        gz(
            tmp_path / f"{aid}.thresholds.bed.gz",
            "#chrom\tstart\tend\tregion\t10X\t30X\nchr1\t0\t100\tunknown\t100\t0\nchr1\t200\t500\tunknown\t300\t300\n",
        )
    # WES target: 1000 bp at depth 50, all >= 30X
    gz(tmp_path / "S1_wes_tumor.target.regions.bed.gz", "chr1\t0\t1000\tT1\t50\n")
    gz(
        tmp_path / "S1_wes_tumor.target.thresholds.bed.gz",
        "#chrom\tstart\tend\tregion\t10X\t30X\nchr1\t0\t1000\tT1\t1000\t1000\n",
    )
    # WGS: genome from the global distribution
    (tmp_path / "S1_wgs_tumor.mosdepth.global.dist.txt").write_text(
        "total\t30\t0.80\ntotal\t10\t0.95\nchr1\t10\t0.95\n"
    )
    (tmp_path / "S1_wgs_tumor.mosdepth.summary.txt").write_text(
        "chrom\tlength\tbases\tmean\tmin\tmax\nchr1\t5000\t150000\t30.00\t0\t90\ntotal\t5000\t150000\t30.00\t0\t90\n"
    )
    result, out = run(tmp_path)
    assert result.returncode == 0, result.stderr
    with open(f"{out}.coverage.tsv") as fh:
        rows = {(r["alignment"], r["region"]): r for r in csv.DictReader(fh, delimiter="\t")}
    common = rows[("S1_wes_tumor", "common")]
    assert (common["bp"], common["mean_depth"], common["pct_ge_10X"], common["pct_ge_30X"]) == (
        "400",
        "35.00",
        "100.00",
        "75.00",
    )
    target = rows[("S1_wes_tumor", "target")]
    assert (target["bp"], target["mean_depth"], target["pct_ge_30X"]) == ("1000", "50.00", "100.00")
    genome = rows[("S1_wgs_tumor", "genome")]
    assert (genome["bp"], genome["mean_depth"], genome["pct_ge_10X"], genome["pct_ge_30X"]) == (
        "5000",
        "30.00",
        "95.00",
        "80.00",
    )
    stats = Path(f"{out}_generalstats_mqc.tsv").read_text()
    assert "plot_type: 'generalstats'" in stats
    assert "S1_wes_tumor\t35.00\t100.00\t75.00\t50.00\t100.00\t100.00" in stats


def test_threshold_mismatch_fails(tmp_path):
    (tmp_path / "alignments.tsv").write_text(
        "alignment\tsample\ttechnology\trole\nA\tS1\twes\ttumor\n"
    )
    gz(tmp_path / "A.regions.bed.gz", "chr1\t0\t100\t20\n")
    gz(
        tmp_path / "A.thresholds.bed.gz",
        "#chrom\tstart\tend\tregion\t10X\nchr1\t0\t100\tunknown\t100\n",
    )
    result, _out = run(tmp_path, "10,30")
    assert result.returncode == 1 and "must match --thresholds" in result.stderr
