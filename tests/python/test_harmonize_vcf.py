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
    assert "dropped 1 records with a symbolic ALT" in result.stderr


def test_wrong_sample_order_fails(tmp_path):
    result, _ = run(tmp_path, "", "--tumor", "N", "--normal", "T")
    assert result.returncode != 0
    assert "expected samples ['N', 'T']" in result.stderr
