"""Unit tests for bin/annotate_variants.py (variant table <-> bedtools results)."""

import csv
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "annotate_variants.py"
TABLE = "chrom\tpos\tref\talt\nchr1\t100\tA\tC\nchr1\t200\tGTT\tG\n"


def run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)


def test_to_bed_covers_the_ref_allele(tmp_path):
    (tmp_path / "v.tsv").write_text(TABLE)
    out = tmp_path / "v.bed"
    proc = run("to-bed", "--table", str(tmp_path / "v.tsv"), "--output", str(out))
    assert proc.returncode == 0, proc.stderr
    assert out.read_text() == "chr1\t99\t100\t1\nchr1\t199\t202\t2\n"


def test_merge_adds_gene_and_stratum_columns(tmp_path):
    (tmp_path / "v.tsv").write_text(TABLE)
    # bedtools output is in sorted order, not table order: rows are matched by name
    (tmp_path / "genes.bed").write_text("chr1\t199\t202\t2\tTP53\nchr1\t99\t100\t1\t.\n")
    (tmp_path / "hp.bed").write_text("chr1\t99\t100\t1\t0\nchr1\t199\t202\t2\t2\n")
    out = tmp_path / "annotated.tsv"
    proc = run(
        "merge",
        "--table",
        str(tmp_path / "v.tsv"),
        "--genes",
        str(tmp_path / "genes.bed"),
        "--stratum",
        f"homopolymer={tmp_path / 'hp.bed'}",
        "--output",
        str(out),
    )
    assert proc.returncode == 0, proc.stderr
    with open(out, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert [(r["gene"], r["homopolymer"]) for r in rows] == [(".", "no"), ("TP53", "yes")]


def test_merge_missing_variant_is_an_error(tmp_path):
    (tmp_path / "v.tsv").write_text(TABLE)
    (tmp_path / "genes.bed").write_text("chr1\t99\t100\t1\tKRAS\n")
    proc = run(
        "merge",
        "--table",
        str(tmp_path / "v.tsv"),
        "--genes",
        str(tmp_path / "genes.bed"),
        "--output",
        str(tmp_path / "a.tsv"),
    )
    assert proc.returncode == 1
    assert "1 variant(s) missing (rows 2)" in proc.stderr


def test_merge_invalid_stratum_name_is_an_error(tmp_path):
    (tmp_path / "v.tsv").write_text(TABLE)
    (tmp_path / "s.bed").write_text("chr1\t99\t100\t1\t0\nchr1\t199\t202\t2\t0\n")
    proc = run(
        "merge",
        "--table",
        str(tmp_path / "v.tsv"),
        "--stratum",
        f"ref={tmp_path / 's.bed'}",
        "--output",
        str(tmp_path / "a.tsv"),
    )
    assert proc.returncode == 1
    assert "stratum name 'ref' is invalid or already a column" in proc.stderr


def test_missing_table_column_is_an_error(tmp_path):
    (tmp_path / "v.tsv").write_text("chrom\tposition\nchr1\t1\n")
    proc = run("to-bed", "--table", str(tmp_path / "v.tsv"), "--output", str(tmp_path / "v.bed"))
    assert proc.returncode == 1
    assert "missing column(s) pos, ref" in proc.stderr
