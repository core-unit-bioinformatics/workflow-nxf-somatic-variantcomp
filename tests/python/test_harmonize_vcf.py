"""Unit tests for bin/harmonize_vcf.py."""

import subprocess
import sys
from pathlib import Path

import pysam
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "harmonize_vcf.py"
sys.path.insert(0, str(SCRIPT.parent))
import harmonize_vcf  # noqa: E402

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=weak_evidence,Description="weak">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT\tN
"""


@pytest.mark.parametrize(
    ("ref", "alt", "expected"),
    [
        ("A", "C", "SNV"),
        ("AC", "GT", "MNV"),
        ("A", "AT", "INS"),
        ("AT", "A", "DEL"),
        ("AT", "GC", "MNV"),
        ("ATG", "C", "COMPLEX"),
    ],
)
def test_variant_type(ref, alt, expected):
    assert harmonize_vcf.variant_type(ref, alt) == expected


def run(tmp_path, body, *extra):
    vcf = tmp_path / "in.vcf"
    vcf.write_text(HEADER + body)
    out = tmp_path / "out.vcf.gz"
    cmd = [sys.executable, str(SCRIPT), "--input", str(vcf), "--output", str(out), *extra]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result, out


def test_tags_gt_and_symbolic(tmp_path):
    body = (
        "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AD:DP\t0/0:70,30:110\t0/0:50,1:51\n"  # GT forced
        "chr1\t20\t.\tA\tAT\t.\tweak_evidence\t.\tGT:AD\t0/1:9,1\t0/0:10,0\n"  # DP from AD
        "chr1\t30\t.\tG\tT\t.\tPASS\t.\tGT:AF:DP\t0/1:0.25:40\t0/0:0:30\n"  # AF fallback
        "chr1\t40\t.\tG\t<NON_REF>\t.\tPASS\t.\tGT:AD\t0/0:10,0\t0/0:10,0\n"  # dropped
    )
    result, out = run(tmp_path, body, "--tumor", "T", "--normal", "N", "--index")
    assert result.returncode == 0, result.stderr
    assert Path(f"{out}.tbi").exists()
    recs = list(pysam.VariantFile(str(out)))
    assert [r.pos for r in recs] == [10, 20, 30]
    first, second, third = recs
    assert first.samples["T"]["GT"] == (0, 1)
    assert first.samples["N"]["GT"] == (0, 0)
    assert first.info["VC_T_VAF"] == pytest.approx(0.3)
    assert first.info["VC_T_DP"] == 110
    assert first.info["VC_T_ALT"] == 30
    assert first.info["VC_N_ALT"] == 1
    assert first.info["VC_TYPE"] == "SNV"
    assert second.info["VC_T_DP"] == 10
    assert second.info["VC_TYPE"] == "INS"
    assert second.info["VC_LEN"] == 1
    assert list(second.filter.keys()) == ["weak_evidence"]
    assert third.info["VC_T_VAF"] == pytest.approx(0.25)
    assert third.info["VC_T_ALT"] == 10
    assert "3 written, 1 dropped (symbolic ALT)" in result.stderr
    # per-record source of the tumor values
    assert [r.info.get("VC_T_SRC") for r in recs] == ["AD", "AD", "AF_DP"]
    assert "1 records without FORMAT/AD" in result.stderr


def test_summary_counts_sources(tmp_path):
    body = (
        "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AD:DP\t0/1:70,30:110\t0/0:50,1:51\n"
        "chr1\t20\t.\tA\tAT\t.\tPASS\t.\tGT:AD\t0/1:9,1\t0/0:.,.\n"
    )
    summary = tmp_path / "summary.tsv"
    result, _ = run(tmp_path, body, "--tumor", "T", "--normal", "N", "--summary", str(summary))
    assert result.returncode == 0, result.stderr
    values = dict(line.split("\t") for line in summary.read_text().splitlines()[1:])
    assert values["records_written"] == "2"
    assert values["T_source_AD"] == "2" and values["T_depth_DP"] == "1" and values["T_depth_AD_sum"] == "1"
    # a normal without reads is only counted
    assert values["N_source_none"] == "1"


@pytest.mark.parametrize(
    ("record", "message"),
    [
        ("chr1\t10\t.\tA\tC,G\t.\tPASS\t.\tGT:AD\t0/1:5,3,2\t0/0:9,0,0\n", "2 ALT alleles; split"),
        ("chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AD\t0/1:5,-3\t0/0:9,0\n", "negative value"),
        ("chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AF:DP\t0/1:1.5:20\t0/0:0:20\n", "outside [0, 1]"),
        ("chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n", "1 of 1 tumor records have no VAF"),
    ],
)
def test_invalid_input_fails_with_message(tmp_path, record, message):
    result, _ = run(tmp_path, record, "--tumor", "T", "--normal", "N")
    assert result.returncode == 1
    assert message in result.stderr


def test_missing_vaf_allowed_up_to_threshold(tmp_path):
    body = (
        "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AD\t0/1:5,5\t0/0:9,0\n"
        "chr1\t20\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n"
    )
    result, _ = run(tmp_path, body, "--tumor", "T", "--normal", "N", "--max-missing-vaf", "0.5")
    assert result.returncode == 0, result.stderr


def test_wrong_sample_order_fails(tmp_path):
    result, _ = run(tmp_path, "", "--tumor", "N", "--normal", "T")
    assert result.returncode != 0
    assert "expected the samples ['N', 'T'] in this order" in result.stderr


def test_nan_af_is_missing_not_an_error(tmp_path):
    # GATK MuTect2 3.8 writes AF=nan for a normal without alt reads; AD is used
    body = "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AD:AF\t0/1:3,5:1\t0/0:2,0:nan\n"
    summary = tmp_path / "summary.tsv"
    result, out = run(tmp_path, body, "--tumor", "T", "--normal", "N", "--summary", str(summary))
    assert result.returncode == 0, result.stderr
    rec = next(iter(pysam.VariantFile(str(out))))
    assert rec.info["VC_N_VAF"] == pytest.approx(0.0)
    assert "N_af_nan\t1" in summary.read_text()


def test_zero_depth_tumor_is_undefined_not_missing(tmp_path):
    # filtered MuTect2 3.8 call without tumor reads: AD 0,0, AF 1, no DP
    body = "chr1\t10\t.\tC\tCA\t.\tweak_evidence\t.\tGT:AD:AF\t0/1:0,0:1\t0/0:5,0:0\n"
    summary = tmp_path / "summary.tsv"
    result, out = run(tmp_path, body, "--tumor", "T", "--normal", "N", "--summary", str(summary))
    assert result.returncode == 0, result.stderr
    rec = next(iter(pysam.VariantFile(str(out))))
    # AF=1 is not used as VAF
    assert "VC_T_VAF" not in rec.info and rec.info["VC_T_ALT"] == 0 and rec.info["VC_T_DP"] == 0
    assert "T_source_zero_depth\t1" in summary.read_text()
    assert "1 records with no tumor reads" in result.stderr
