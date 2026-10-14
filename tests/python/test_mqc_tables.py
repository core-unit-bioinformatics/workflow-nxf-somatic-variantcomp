"""Unit tests for bin/mqc_tables.py (MultiQC tables of the step summaries)."""

import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "mqc_tables.py"


def run(tmp_path, *args):
    out = tmp_path / "x"
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            *args,
            "--section",
            "test",
            "--title",
            "Test",
            "--prefix",
            str(out),
        ],
        capture_output=True,
        text=True,
    )
    lines = Path(f"{out}_mqc.tsv").read_text().splitlines() if proc.returncode == 0 else []
    return proc, lines


def test_key_value_counts_messages(tmp_path):
    a = tmp_path / "S1_wes.vcf.reference_check.tsv"
    a.write_text("key\tvalue\nstatus\tok\nrecords\t10\nwarning\tw1\nwarning\tw2\n")
    b = tmp_path / "S1_wgs.vcf.reference_check.tsv"
    b.write_text("key\tvalue\nstatus\terror\nrecords\t5\nerror\te1\n")
    proc, lines = run(tmp_path, "key-value", str(a), str(b), "--strip", ".reference_check.tsv")
    assert proc.returncode == 0, proc.stderr
    assert lines[4] == "id\twarnings\terrors\tstatus\trecords"
    assert lines[5:] == ["S1_wes.vcf\t2\t0\tok\t10", "S1_wgs.vcf\t0\t1\terror\t5"]


def test_norm_logs(tmp_path):
    log = tmp_path / "S1_wes.norm.log"
    log.write_text(
        "Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:\t11/2/0/1/0/0/0\n"
    )
    proc, lines = run(tmp_path, "norm-logs", str(log), "--strip", ".norm.log")
    assert proc.returncode == 0, proc.stderr
    assert lines[-1] == "S1_wes\t11\t2\t0\t1\t0\t0\t0"


def test_norm_log_without_summary_fails(tmp_path):
    log = tmp_path / "bad.log"
    log.write_text("nothing\n")
    proc, _ = run(tmp_path, "norm-logs", str(log))
    assert proc.returncode == 1 and "no summary line of bcftools norm" in proc.stderr


def test_table_rows_and_mismatching_header(tmp_path):
    a = tmp_path / "a.tsv"
    a.write_text("region\ttype\tbp\nwes\ttarget\t100\n")
    proc, lines = run(tmp_path, "table", str(a))
    assert proc.returncode == 0, proc.stderr
    assert lines[-2:] == ["id\ttype\tbp", "wes\ttarget\t100"]
    b = tmp_path / "b.tsv"
    b.write_text("other\tcolumns\n")
    proc, _ = run(tmp_path, "table", str(a), str(b))
    assert proc.returncode == 1 and "header differs" in proc.stderr


def test_rows_sorted_by_file_name_not_staging_folder(tmp_path):
    # Nextflow stages inputs as in1/, in2/, ... in arrival order: the row order must
    # only depend on the file names
    for folder, name in (
        ("in1", "b.vcf.reference_check.tsv"),
        ("in2", "a.vcf.reference_check.tsv"),
    ):
        (tmp_path / folder).mkdir()
        (tmp_path / folder / name).write_text("key\tvalue\nstatus\tok\n")
    paths = [
        str(tmp_path / "in1" / "b.vcf.reference_check.tsv"),
        str(tmp_path / "in2" / "a.vcf.reference_check.tsv"),
    ]
    proc, lines = run(tmp_path, "key-value", *paths, "--strip", ".reference_check.tsv")
    assert proc.returncode == 0, proc.stderr
    assert [line.split("\t")[0] for line in lines[5:]] == ["a.vcf", "b.vcf"]
