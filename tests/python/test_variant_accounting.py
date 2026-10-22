"""Unit tests for bin/variant_accounting.py."""

import subprocess
import sys
from pathlib import Path

import pysam

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "variant_accounting.py"
HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID={c1},length=1000>\n##contig=<ID={mt},length=100>\n"
    '##FILTER=<ID=weak,Description="w">\n'
    '##INFO=<ID=VC_T_VAF,Number=1,Type=Float,Description="v">\n'
    '##INFO=<ID=VC_T_DP,Number=1,Type=Integer,Description="d">\n'
    '##INFO=<ID=VC_T_ALT,Number=1,Type=Integer,Description="a">\n'
    '##INFO=<ID=VC_TYPE,Number=1,Type=String,Description="t">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def vcf(path, records, c1="chr1", mt="chrM"):
    """bgzipped VCF with the given record lines."""
    plain = path.with_suffix("")
    plain.write_text(HEADER.format(c1=c1, mt=mt) + "".join(r + "\n" for r in records))
    pysam.tabix_compress(str(plain), str(path), force=True)
    return path


def run_accounting(tmp_path, compare_records):
    """Run the script on a small case; return (exit code, accounting, removed, checks)."""
    # input: 1 PASS SNV, 1 non-PASS, 1 on MT, 1 gVCF block, 1 low-VAF SNV, 1 outside region
    vcf(
        tmp_path / "in.vcf.gz",
        [
            "1\t100\t.\tA\tC\t.\tPASS\t.",
            "1\t110\t.\tA\tC\t.\tweak\t.",
            "MT\t10\t.\tA\tC\t.\tPASS\t.",
            "1\t120\t.\tA\t<NON_REF>\t.\t.\t.",
            "1\t130\t.\tA\tC\t.\tPASS\t.",
            "1\t900\t.\tAT\tA\t.\tPASS\t.",
        ],
        c1="1",
        mt="MT",
    )
    (tmp_path / "rename.tsv").write_text("1\tchr1\nMT\tchrM\n")
    (tmp_path / "genome.bed").write_text("chr1\t0\t1000\n")
    (tmp_path / "norm.log").write_text(
        "Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:\t3/0/0/0/0/0/0\n"
    )
    (tmp_path / "summary.tsv").write_text(
        "key\tvalue\nrecords_read\t3\nrecords_written\t3\nrecords_dropped_symbolic\t0\n"
    )
    harmonized = [
        "chr1\t100\t.\tA\tC\t.\tPASS\tVC_T_VAF=0.3;VC_T_DP=50;VC_T_ALT=15;VC_TYPE=SNV",
        "chr1\t130\t.\tA\tC\t.\tPASS\tVC_T_VAF=0.01;VC_T_DP=50;VC_T_ALT=1;VC_TYPE=SNV",
        "chr1\t900\t.\tAT\tA\t.\tPASS\tVC_T_VAF=0.2;VC_T_DP=50;VC_T_ALT=10;VC_TYPE=DEL",
    ]
    vcf(tmp_path / "harmonized.vcf.gz", harmonized)
    vcf(tmp_path / "compare.vcf.gz", compare_records)
    (tmp_path / "regions.bed").write_text("chr1\t50\t500\n")
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--sample",
            "S1",
            "--technology",
            "wes",
            "--input-vcf",
            str(tmp_path / "in.vcf.gz"),
            "--rename",
            str(tmp_path / "rename.tsv"),
            "--genome-bed",
            str(tmp_path / "genome.bed"),
            "--norm-log",
            str(tmp_path / "norm.log"),
            "--harmonize-summary",
            str(tmp_path / "summary.tsv"),
            "--harmonized",
            str(tmp_path / "harmonized.vcf.gz"),
            "--compare",
            str(tmp_path / "compare.vcf.gz"),
            "--regions",
            str(tmp_path / "regions.bed"),
            "--pass-only",
            "--min-vaf",
            "0.05",
            "--prefix",
            str(tmp_path / "out"),
        ],
        capture_output=True,
        text=True,
    )

    def read(name):
        """Rows (without header) of an output table."""
        lines = (tmp_path / f"out.{name}.tsv").read_text().splitlines()[1:]
        return [line.split("\t") for line in lines]

    return result, read("accounting"), read("removed"), read("checks")


def test_steps_reasons_and_checks(tmp_path):
    compare = ["chr1\t100\t.\tA\tC\t.\tPASS\tVC_T_VAF=0.3;VC_T_DP=50;VC_T_ALT=15;VC_TYPE=SNV"]
    result, accounting, removed, checks = run_accounting(tmp_path, compare)
    assert result.returncode == 0, result.stderr
    reasons = {(r[3], r[5]): int(r[6]) for r in accounting}
    assert reasons[("input", ".")] == 6
    assert reasons[("pass_split", "filter_not_pass")] == 1
    assert reasons[("contigs", "mitochondrial")] == 1
    assert reasons[("gvcf_blocks", "gvcf_block")] == 1
    assert reasons[("filters", "vaf_below_min")] == 1
    assert reasons[("filters", "outside_region")] == 1
    assert reasons[("filters", ".")] == 1
    assert sorted((r[2], r[3]) for r in removed) == [
        ("contigs", "mitochondrial"),
        ("filters", "outside_region"),
        ("filters", "vaf_below_min"),
        ("gvcf_blocks", "gvcf_block"),
        ("pass_split", "filter_not_pass"),
    ]
    assert all(c[5] == "ok" for c in checks), checks


def test_compare_vcf_that_does_not_match_the_filters_fails_a_check(tmp_path):
    # the low-VAF record is in the compared VCF although --min-vaf removes it
    compare = [
        "chr1\t100\t.\tA\tC\t.\tPASS\tVC_T_VAF=0.3;VC_T_DP=50;VC_T_ALT=15;VC_TYPE=SNV",
        "chr1\t130\t.\tA\tC\t.\tPASS\tVC_T_VAF=0.01;VC_T_DP=50;VC_T_ALT=1;VC_TYPE=SNV",
    ]
    result, _accounting, _removed, checks = run_accounting(tmp_path, compare)
    # the script itself doesn't stop: the checks are evaluated at the end of the pipeline
    assert result.returncode == 0
    failed = [c[2] for c in checks if c[5] == "FAIL"]
    assert "compared VCF = harmonized records passing the filters and region" in failed
