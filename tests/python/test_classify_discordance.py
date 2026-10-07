"""Unit tests for bin/classify_discordance.py."""

import csv
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN))
import classify_discordance  # noqa: E402

CALL_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=weak_evidence,Description="weak">
##INFO=<ID=VC_T_VAF,Number=1,Type=Float,Description="vaf">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""
PILEUP_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tbam
"""
VARIANT_COLUMNS = ["sample", "truth", "query", "status", "chrom", "pos", "ref", "alt",
                   "truth_vaf", "query_vaf"]  # fmt: skip


def write(path, text):
    path.write_text(text)
    return str(path)


@pytest.fixture()
def run(tmp_path):
    variants = [
        # status, pos, ref, alt, truth_vaf
        ("TP", 100, "A", "C", "0.30"),
        ("FN", 200, "A", "C", "0.30"),  # filtered in query
        ("FN", 300, "A", "C", "0.30"),  # alt reads in normal -> germline
        ("FN", 400, "A", "C", "0.30"),  # low coverage in query
        ("FN", 500, "A", "C", "0.02"),  # 1 alt read in query -> subthreshold; low power
        ("FN", 600, "A", "C", "0.30"),  # no alt read in query
        ("FP", 700, "A", "C", "NA"),  # no alt read in truth
    ]
    with open(tmp_path / "variants.tsv", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(VARIANT_COLUMNS)
        for status, pos, ref, alt, vaf in variants:
            writer.writerow(["S1", "wgs", "wes", status, "chr1", pos, ref, alt, vaf, "NA"])
    truth_norm = write(tmp_path / "truth.vcf", CALL_HEADER)
    query_norm = write(
        tmp_path / "query.vcf",
        CALL_HEADER + "chr1\t200\t.\tA\tC\t.\tweak_evidence\tVC_T_VAF=0.05\n",
    )
    pile = "chr1\t{pos}\t.\tA\tC,<*>\t.\t.\t.\tAD:DP\t{ref},{alt},0:{dp}\n"
    query_tumor = write(
        tmp_path / "query_tumor.vcf",
        PILEUP_HEADER
        + "".join(
            pile.format(pos=p, ref=r, alt=a, dp=r + a)
            for p, r, a in [
                (100, 70, 30),
                (200, 95, 5),
                (300, 100, 0),
                (400, 5, 0),
                (500, 99, 1),
                (600, 200, 0),
            ]  # fmt: skip
        ),
    )
    truth_tumor = write(tmp_path / "truth_tumor.vcf", PILEUP_HEADER + pile.format(
        pos=700, ref=50, alt=0, dp=50))  # fmt: skip
    normal = write(tmp_path / "normal.vcf", PILEUP_HEADER + pile.format(
        pos=300, ref=50, alt=48, dp=98))  # fmt: skip
    prefix = tmp_path / "out"
    subprocess.run(
        [sys.executable, str(BIN / "classify_discordance.py"), "--variants",
         str(tmp_path / "variants.tsv"), "--truth-normalized", truth_norm, "--query-normalized",
         query_norm, "--truth-tumor-evidence", truth_tumor, "--query-tumor-evidence", query_tumor,
         "--normal-evidence", normal, "--prefix", str(prefix)],
        check=True,
    )  # fmt: skip
    with open(f"{prefix}.discordance.tsv") as fh:
        rows = {int(r["pos"]): r for r in csv.DictReader(fh, delimiter="\t")}
    with open(f"{prefix}.power.tsv") as fh:
        power = next(csv.DictReader(fh, delimiter="\t"))
    return rows, power


def test_classes(run):
    rows, _ = run
    assert {pos: r["class"] for pos, r in rows.items()} == {
        100: "concordant",
        200: "present_but_filtered",
        300: "germline_in_normal",
        400: "low_coverage",
        500: "subthreshold_evidence",
        600: "no_evidence",
        700: "no_evidence",
    }
    assert rows[200]["other_filter"] == "weak_evidence"
    assert rows[300]["normal_alt"] == "48"


def test_power(run):
    rows, power = run
    assert rows[100]["detectable"] == "yes"
    assert rows[400]["detectable"] == "no"  # depth 5
    assert rows[500]["detectable"] == "no"  # VAF 2 % at depth 100
    assert rows[700]["query_power"] == "NA"  # FP: not a truth variant
    # detectable truth variants: 100 (TP), 200, 300, 600 (FN) -> recall 1/4
    assert (power["tp_baseline"], power["fn"], power["recall"]) == ("1", "3", "0.2500")


def test_detection_power_function():
    assert classify_discordance.detection_power(0, 0.3, 3) == 0.0
    assert classify_discordance.detection_power(100, 0.3, 3) == pytest.approx(1.0)
    # P(X >= 1) for n=10, p=0.1
    assert classify_discordance.detection_power(10, 0.1, 1) == pytest.approx(1 - 0.9**10)
