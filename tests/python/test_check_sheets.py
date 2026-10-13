"""Unit tests for bin/check_sheets.py."""

import csv
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "check_sheets.py"
HEADER = (
    "sample,technology,library_type,tumor_id,normal_id,caller,vcf,tumor_bam,normal_bam,target_bed\n"
)


def row(sample, tech, lib="wgs", bed="n/a", tumor="T", normal="N", bam="n/a"):
    return f"{sample},{tech},{lib},{tumor},{normal},dragen,vcf/{sample}_{tech}.vcf.gz,{bam},n/a,{bed}\n"


def run(tmp_path, text, command="samplesheet", *extra):
    sheet = tmp_path / "sheet.csv"
    sheet.write_text(text)
    out = tmp_path / "out.csv"
    cmd = [
        sys.executable,
        str(SCRIPT),
        command,
        "--sheet",
        str(sheet),
        "--base-dir",
        str(tmp_path),
        "--output",
        str(out),
        *extra,
    ]
    if command == "samplesheet":
        cmd += ["--truth-technology", "wgs"]
    return subprocess.run(cmd, capture_output=True, text=True), out


VALID = HEADER + row("S1", "wgs") + row("S1", "wes", "capture", "bed/wes.bed")


def test_valid_samplesheet_resolves_paths(tmp_path):
    result, out = run(tmp_path, VALID)
    assert result.returncode == 0, result.stderr
    rows = list(csv.DictReader(out.open()))
    assert [r["is_truth"] for r in rows] == ["true", "false"]
    assert rows[1]["target_bed"] == str(tmp_path / "bed" / "wes.bed")
    assert rows[0]["target_bed"] == "" and rows[0]["tumor_bam"] == ""


@pytest.mark.parametrize(
    ("text", "message"),
    [
        (
            HEADER + row("S1", "wgs") + row("S1", "wgs") + row("S1", "wes", "capture", "x.bed"),
            "occur twice",
        ),
        (HEADER + row("S1", "wgs") + row("S2", "wgs"), "at least two technologies"),
        (
            HEADER
            + row("S1", "wgs")
            + row("S1", "wes", "capture", "x.bed")
            + row("S2", "wes", "capture", "x.bed"),
            "sample 'S2' has no row for the truth technology 'wgs'",
        ),
        (
            HEADER + row("S1", "wgs") + row("S1", "wes", "capture", "x.bed") + row("S2", "wgs"),
            "sample 'S2' has no technology to compare",
        ),
        (
            HEADER
            + row("S1", "wgs")
            + row("S1", "wes", "capture", "a.bed")
            + row("S2", "wgs")
            + row("S2", "wes", "capture", "b.bed"),
            "must use one target_bed",
        ),
        (HEADER + row("S1", "wgs") + row("S1", "panel", "amplicon"), "needs a target_bed"),
        (
            HEADER + row("S1", "wgs") + row("S1", "wes", "capture", "x.bed", tumor="T", normal="T"),
            "tumor_id and normal_id are identical",
        ),
        (HEADER + row("S1", "wgs") + row("S1", "wes", "exome", "x.bed"), "library_type 'exome'"),
        (
            HEADER.replace("target_bed", "target_bed,vcf_index")
            + row("S1", "wgs").rstrip()
            + ",x.tbi\n",
            "column 'vcf_index' is no longer supported",
        ),
        (VALID.replace(",", "\t"), "tab-separated"),
    ],
)
def test_invalid_samplesheets(tmp_path, text, message):
    result, _ = run(tmp_path, text)
    assert result.returncode == 1
    assert message in result.stderr


def test_truth_technology_missing(tmp_path):
    sheet = HEADER + row("S1", "panel", "amplicon", "x.bed") + row("S1", "wes", "capture", "x.bed")
    result, _ = run(tmp_path, sheet)
    assert result.returncode == 1
    assert "truth technology 'wgs' (--truth_technology) is not in the samplesheet" in result.stderr


def test_all_errors_reported_at_once(tmp_path):
    result, _ = run(tmp_path, HEADER + row("S1 x", "WGS", "foo"))
    assert result.returncode == 1
    assert result.stderr.count("ERROR:") >= 4  # sample, technology, library_type + summary


def test_stratification(tmp_path):
    result, out = run(
        tmp_path,
        "name,bed\nhomopolymer,giab/hp.bed\ntandem_repeat,/abs/tr.bed.gz\n",
        "stratification",
    )
    assert result.returncode == 0, result.stderr
    rows = list(csv.DictReader(out.open()))
    assert rows[0]["bed"] == str(tmp_path / "giab" / "hp.bed")
    assert rows[1]["bed"] == "/abs/tr.bed.gz"
    result, _ = run(
        tmp_path, "name,bed\nhomo-polymer,x.bed\nhomo-polymer,y.bed\n", "stratification"
    )
    assert result.returncode == 1
    assert "may only contain" in result.stderr and "used twice" in result.stderr


def test_known(tmp_path):
    samplesheet = tmp_path / "validated.csv"
    samplesheet.write_text("sample,technology\nS1,wgs\nS1,wes\n")
    text = (
        "name,sample,mode,variants,regions\n"
        "clinical,S1,required,known/a.csv,n/a\n"
        "truth,all,benchmark,t/snv.vcf.gz;t/indel.vcf.gz,t/hc.bed\n"
    )
    result, out = run(tmp_path, text, "known", "--samplesheet", str(samplesheet))
    assert result.returncode == 0, result.stderr
    rows = list(csv.DictReader(out.open()))
    assert rows[1]["variants"].split(";") == [
        str(tmp_path / "t" / "snv.vcf.gz"),
        str(tmp_path / "t" / "indel.vcf.gz"),
    ]
    bad = "name,sample,mode,variants,regions\nx,S9,check,a.txt,n/a\n"
    result, _ = run(tmp_path, bad, "known", "--samplesheet", str(samplesheet))
    assert result.returncode == 1
    for message in ("sample 'S9' is not in the samplesheet", "mode 'check'", "must end with"):
        assert message in result.stderr
