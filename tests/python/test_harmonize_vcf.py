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
    assert "3 written, 1 dropped (symbolic ALT), 0 dropped (overlapping allele *)" in result.stderr
    # per-record source of the tumor values
    # auto: FORMAT/AF where present (third record), else the allelic depths
    assert [r.info.get("VC_T_SRC") for r in recs] == ["AD", "AD", "AF"]


def test_overlapping_allele_is_dropped_and_counted(tmp_path):
    # 'A>*' halves of split records (a spanning deletion of the caller, or of an atomized
    # complex call): no sequence change of their own, AD of '*' may be missing
    body = (
        "chr1\t10\t.\tACGTA\tA\t.\tPASS\t.\tGT:AD\t0/1:20,10\t0/0:30,0\n"
        "chr1\t13\t.\tT\tTGTATAAGAGG\t.\tPASS\t.\tGT:AD\t0/1:20,5\t0/0:30,0\n"
        "chr1\t13\t.\tT\t*\t.\tPASS\t.\tGT:AD\t0/1:20,10\t0/0:30,0\n"
        "chr1\t20\t.\tG\t*\t.\tPASS\t.\tGT:AD\t0/1:5,.\t0/0:9,.\n"
    )
    summary = tmp_path / "summary.tsv"
    result, out = run(tmp_path, body, "--tumor", "T", "--normal", "N", "--summary", str(summary))
    assert result.returncode == 0, result.stderr
    recs = list(pysam.VariantFile(str(out)))
    assert [(r.pos, r.alts[0]) for r in recs] == [(10, "A"), (13, "TGTATAAGAGG")]
    assert "2 written, 0 dropped (symbolic ALT), 2 dropped (overlapping allele *)" in result.stderr
    values = dict(line.split("\t") for line in summary.read_text().splitlines()[1:])
    assert values["records_dropped_overlap"] == "2"
    assert values["records_dropped_symbolic"] == "0"
    # the '*' records don't count as records without VAF
    assert values["vaf_available"] == "yes"


def test_unsplit_overlapping_allele_fails(tmp_path):
    # bcftools norm --atomize with the default --atom-overlaps '*' writes such records
    body = "chr1\t13\t.\tT\tTGTATAAGAGG,*\t.\tPASS\t.\tGT:AD\t1/2:20,5,10\t0/0:30,0,0\n"
    result, _ = run(tmp_path, body, "--tumor", "T", "--normal", "N")
    assert result.returncode == 1
    assert "chr1:13 T>TGTATAAGAGG,*: 2 ALT alleles; split multiallelic" in result.stderr


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
    assert (
        values["T_source_AD"] == "2"
        and values["T_depth_DP"] == "1"
        and values["T_depth_AD_sum"] == "1"
    )
    # a normal without reads is only counted
    assert values["N_source_none"] == "1"


@pytest.mark.parametrize(
    ("record", "message"),
    [
        ("chr1\t10\t.\tA\tC,G\t.\tPASS\t.\tGT:AD\t0/1:5,3,2\t0/0:9,0,0\n", "2 ALT alleles; split"),
        ("chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AD\t0/1:5,-3\t0/0:9,0\n", "negative value"),
        ("chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AF:DP\t0/1:1.5:20\t0/0:0:20\n", "outside [0, 1]"),
    ],
)
def test_invalid_input_fails_with_message(tmp_path, record, message):
    result, _ = run(tmp_path, record, "--tumor", "T", "--normal", "N")
    assert result.returncode == 1
    assert message in result.stderr


def test_missing_vaf_above_threshold_marks_vaf_unavailable(tmp_path):
    body = (
        "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AD\t0/1:5,5\t0/0:9,0\n"
        "chr1\t20\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n"
    )
    # 1 of 2 without VAF: allowed with 0.5, VAF available
    result, out = run(tmp_path, body, "--tumor", "T", "--normal", "N", "--max-missing-vaf", "0.5")
    assert result.returncode == 0, result.stderr
    assert "variantcomp_vaf" not in str(pysam.VariantFile(str(out)).header)
    # above the default 0.05: a warning and the header flag, no error
    result, out = run(tmp_path, body, "--tumor", "T", "--normal", "N")
    assert result.returncode == 0, result.stderr
    assert "1 of 2 tumor records have no VAF" in result.stderr
    assert "##variantcomp_vaf=unavailable" in str(pysam.VariantFile(str(out)).header)


@pytest.mark.parametrize(
    ("mode", "expected"),
    [("auto", (0.25, "AF")), ("af", (0.25, "AF")), ("ad", (0.3, "AD"))],
)
def test_vaf_source_modes(tmp_path, mode, expected):
    # AF 0.25 and AD 70,30 disagree: the mode decides
    body = "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AD:AF:DP\t0/1:70,30:0.25:100\t0/0:50,0:0:50\n"
    result, out = run(tmp_path, body, "--tumor", "T", "--normal", "N", "--vaf-source", mode)
    assert result.returncode == 0, result.stderr
    rec = next(iter(pysam.VariantFile(str(out))))
    assert (rec.info["VC_T_VAF"], rec.info["VC_T_SRC"]) == (pytest.approx(expected[0]), expected[1])
    # alt reads always from AD, depth from DP
    assert (rec.info["VC_T_ALT"], rec.info["VC_T_DP"]) == (30, 100)


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
