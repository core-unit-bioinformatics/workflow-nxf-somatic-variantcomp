"""Unit tests for bin/known_summary.py (summary of the required known variants)."""

import csv
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "known_summary.py"
VARIANTS = (
    "sample\ttechnology\tlabel\tgene\tchrom\tpos\tref\talt\tstatus\n"
    "S1\twes\tV1\tKRAS\tchr1\t100\tA\tC\tPASS\n"
    "S1\twes\tV2\tTP53\tchr1\t50\tG\tT\tabsent\n"
)
SUMMARY = "sample\ttechnology\trequired\tpassed\tfailed\tall_passed\nS1\twes\t2\t1\t1\tno\n"


def run(tmp_path, *extra):
    (tmp_path / "v.tsv").write_text(VARIANTS)
    (tmp_path / "s.tsv").write_text(SUMMARY)
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--variants",
            f"clinical={tmp_path / 'v.tsv'}",
            "--summary",
            f"clinical={tmp_path / 's.tsv'}",
            "--prefix",
            str(tmp_path / "known"),
            *extra,
        ],
        capture_output=True,
        text=True,
    )


def test_summary_tables_and_warnings(tmp_path):
    proc = run(tmp_path)
    assert proc.returncode == 0, proc.stderr
    with open(tmp_path / "known.required_variants.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    # sorted by position; the list name is the first column
    assert [(r["list"], r["label"]) for r in rows] == [("clinical", "V2"), ("clinical", "V1")]
    assert (
        "required variant V2 chr1:50 G>T (TP53) of list 'clinical' not found in S1 wes: absent"
        in proc.stderr
    )
    mqc = (tmp_path / "known_required_mqc.tsv").read_text().splitlines()
    assert mqc[-1] == "clinical S1 wes\tclinical\tS1\twes\t2\t1\t1\tno"


def test_fail_on_missing_exits_after_writing_outputs(tmp_path):
    proc = run(tmp_path, "--fail-on-missing")
    assert proc.returncode == 1
    assert "1 required variant check(s) failed" in proc.stderr
    assert (tmp_path / "known.required_summary.tsv").exists()


def test_missing_column_is_an_error(tmp_path):
    (tmp_path / "bad.tsv").write_text("sample\ttechnology\nS1\twes\n")
    (tmp_path / "s.tsv").write_text(SUMMARY)
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--variants",
            f"x={tmp_path / 'bad.tsv'}",
            "--summary",
            f"x={tmp_path / 's.tsv'}",
            "--prefix",
            str(tmp_path / "k"),
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 1
    assert "missing column(s) label" in proc.stderr
