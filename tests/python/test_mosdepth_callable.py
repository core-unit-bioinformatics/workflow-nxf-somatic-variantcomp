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


def test_clip_keeps_only_the_parts_inside_the_clip_region(tmp_path):
    # clip BED unsorted, with overlapping and adjacent intervals (merged by the script):
    # chr1 50-150 + 140-160 -> 50-160, 300-350 + 350-380 -> 300-380; chr2 is not in it
    (tmp_path / "clip.bed").write_text(
        "chr1\t300\t350\nchr1\t50\t150\nchr1\t140\t160\nchr1\t350\t380\n"
    )
    proc, lines = run(
        tmp_path,
        "base",
        QUANTIZED,
        REGIONS,
        rename="1\tchr1\n2\tchr2\n",
        extra=["--clip", str(tmp_path / "clip.bed")],
    )
    assert proc.returncode == 0, proc.stderr
    # callable chr1 100-400 is cut into its parts inside the clip region; chr2 0-50 dropped
    assert lines == ["chr1\t100\t160", "chr1\t300\t380"]
    assert (
        "1 callable intervals outside the clip region dropped; 2 intervals written" in proc.stderr
    )


def test_clip_gives_the_same_intersection_as_clipping_afterwards(tmp_path):
    # (C1 n C2) n K == (C1 n K) n (C2 n K), checked base by base on random intervals
    import random

    sys.path.insert(0, str(SCRIPT.parent))
    from mosdepth_callable import clip_interval, read_clip_region

    rng = random.Random(7)

    def random_bed(n):
        """n random intervals on chr1 (0-1000), as a set of covered bases and a list."""
        intervals = []
        for _ in range(n):
            start = rng.randrange(0, 990)
            intervals.append((start, start + rng.randrange(1, 60)))
        return {b for s, e in intervals for b in range(s, e)}, intervals

    k_bases, k_intervals = random_bed(15)
    clip_bed = tmp_path / "k.bed"
    clip_bed.write_text("".join(f"chr1\t{s}\t{e}\n" for s, e in k_intervals))

    region = read_clip_region(clip_bed)
    for _ in range(20):
        c1_bases, c1 = random_bed(30)
        c2_bases, c2 = random_bed(30)
        clipped = [
            {
                b
                for s, e in bed
                for ps, pe in clip_interval("chr1", s, e, region)
                for b in range(ps, pe)
            }
            for bed in (c1, c2)
        ]
        assert clipped[0] & clipped[1] == (c1_bases & c2_bases) & k_bases


def test_region_mode_with_clip(tmp_path):
    # region mode: the regions of mosdepth --by already lie inside the common region; the
    # clip keeps them unchanged
    (tmp_path / "clip.bed").write_text("chr1\t0\t200\nchr1\t300\t400\n")
    proc, lines = run(
        tmp_path,
        "region",
        QUANTIZED,
        REGIONS,
        rename="1\tchr1\n",
        extra=["--min-coverage", "10", "--clip", str(tmp_path / "clip.bed")],
    )
    assert proc.returncode == 0, proc.stderr
    assert lines == ["chr1\t0\t200"]
