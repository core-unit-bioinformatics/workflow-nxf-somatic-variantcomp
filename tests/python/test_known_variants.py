"""Unit tests for bin/known_variants.py."""

import csv
import random
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "known_variants.py"
CALL_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=2000>
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=weak_evidence,Description="weak">
##INFO=<ID=VC_T_VAF,Number=1,Type=Float,Description="vaf">
##INFO=<ID=VC_T_DP,Number=1,Type=Integer,Description="dp">
##INFO=<ID=VC_T_ALT,Number=1,Type=Integer,Description="alt">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


@pytest.fixture()
def ref(tmp_path):
    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(2000))
    fasta = tmp_path / "ref.fa"
    fasta.write_text(f">chr1\n{seq}\n")
    pysam.faidx(str(fasta))
    return fasta, seq


def alt_of(base):
    return "A" if base != "A" else "C"


def run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)


def test_prepare_tsv_renames_and_sorts(tmp_path, ref):
    fasta, seq = ref
    tsv = tmp_path / "known.tsv"
    tsv.write_text(
        "chrom\tpos\tref\talt\tgene\tlabel\n"
        f"1\t500\t{seq[499]}\t{alt_of(seq[499])}\tTP53\tR175H\n"
        f"chr1\t100\t{seq[99]}\t{alt_of(seq[99])}\tKRAS\tG12D\n"
    )
    out = tmp_path / "known.vcf.gz"
    result = run("prepare", "--input", str(tsv), "--fai", f"{fasta}.fai", "--output", str(out))
    assert result.returncode == 0, result.stderr
    recs = list(pysam.VariantFile(str(out)))
    assert [(r.chrom, r.pos, r.id, r.info["KV_GENE"]) for r in recs] == [
        ("chr1", 100, "G12D", "KRAS"),
        ("chr1", 500, "R175H", "TP53"),
    ]


def test_prepare_unknown_contig_and_position_fail(tmp_path, ref):
    fasta, seq = ref
    tsv = tmp_path / "known.tsv"
    tsv.write_text(
        f"chrom\tpos\tref\talt\tlabel\nchrX\t100\tA\tC\tv1\nchr1\t2500\tA\tC\tv2\nchr1\t10\t{seq[9]}\t<DEL>\tv3\n"
    )
    result = run(
        "prepare",
        "--input",
        str(tsv),
        "--fai",
        f"{fasta}.fai",
        "--output",
        str(tmp_path / "k.vcf.gz"),
    )
    assert result.returncode != 0
    for message in (
        "contig 'chrX' is not in the reference",
        "position beyond the end of chr1",
        "REF/ALT must be",
    ):
        assert message in result.stderr


def test_check_statuses(tmp_path, ref):
    fasta, seq = ref
    variants = {pos: (seq[pos - 1], alt_of(seq[pos - 1])) for pos in (100, 200, 300, 400)}
    tsv = tmp_path / "known.tsv"
    tsv.write_text(
        "chrom\tpos\tref\talt\n"
        + "".join(f"chr1\t{p}\t{r}\t{a}\n" for p, (r, a) in variants.items())
    )
    known = tmp_path / "known.vcf.gz"
    assert (
        run(
            "prepare", "--input", str(tsv), "--fai", f"{fasta}.fai", "--output", str(known)
        ).returncode
        == 0
    )

    def call(pos, filt="PASS"):
        r, a = variants[pos]
        return f"chr1\t{pos}\t.\t{r}\t{a}\t.\t{filt}\tVC_T_VAF=0.3;VC_T_DP=100;VC_T_ALT=30\n"

    # 100: compared (PASS); 200: filtered; 300: not called; 400: outside the region
    normalized = tmp_path / "norm.vcf"
    normalized.write_text(CALL_HEADER + call(100) + call(200, "weak_evidence") + call(400))
    compare = tmp_path / "compare.vcf"
    compare.write_text(CALL_HEADER + call(100))
    regions = tmp_path / "regions.bed"
    regions.write_text("chr1\t0\t350\n")
    prefix = tmp_path / "S1"
    result = run(
        "check",
        "--known",
        str(known),
        "--sample",
        "S1",
        "--technology",
        "wes",
        "--normalized",
        str(normalized),
        "--compare",
        str(compare),
        "--regions",
        str(regions),
        "--prefix",
        str(prefix),
    )
    assert result.returncode == 0, result.stderr
    with open(f"{prefix}.known_variants.tsv") as fh:
        rows = {int(r["pos"]): r for r in csv.DictReader(fh, delimiter="\t")}
    assert {p: r["status"] for p, r in rows.items()} == {
        100: "PASS",
        200: "FAIL_filtered",
        300: "FAIL_not_called",
        400: "FAIL_outside_region",
    }
    assert rows[200]["filter"] == "weak_evidence"
    with open(f"{prefix}.known_summary.tsv") as fh:
        summary = next(csv.DictReader(fh, delimiter="\t"))
    assert (summary["required"], summary["passed"], summary["all_passed"]) == ("4", "1", "no")


@pytest.mark.parametrize(
    ("row", "message"),
    [
        ("chr1\tx10\tA\tC\t0.3", "line 2: pos 'x10' is not a positive integer"),
        ("chr1\t10\tA\tC\thigh", "line 2: vaf 'high' is not a number"),
        ("chr1\t10\tA\tC\t1.4", "line 2: vaf 1.4 is outside [0, 1]"),
    ],
)
def test_prepare_invalid_tsv_values_fail(tmp_path, ref, row, message):
    tsv = tmp_path / "known.tsv"
    tsv.write_text("chrom\tpos\tref\talt\tvaf\n" + row + "\n")
    result = run(
        "prepare",
        "--input",
        str(tsv),
        "--fai",
        f"{ref[0]}.fai",
        "--output",
        str(tmp_path / "k.vcf.gz"),
    )
    assert result.returncode == 1
    assert message in result.stderr


def test_check_evidence_of_unknown_technology_fails(tmp_path, ref):
    known = tmp_path / "known.vcf"
    known.write_text(
        CALL_HEADER.replace(
            "#CHROM",
            '##INFO=<ID=KV_GENE,Number=1,Type=String,Description="g">\n'
            '##INFO=<ID=KV_LABEL,Number=1,Type=String,Description="l">\n#CHROM',
        )
        + "chr1\t10\tv1\tA\tC\t.\t.\tKV_GENE=G\n"
    )
    calls = tmp_path / "calls.vcf"
    calls.write_text(CALL_HEADER)
    bed = tmp_path / "r.bed"
    bed.write_text("chr1\t0\t100\n")
    result = run(
        "check",
        "--known",
        str(known),
        "--sample",
        "S1",
        "--technology",
        "wes",
        "--normalized",
        str(calls),
        "--compare",
        str(calls),
        "--regions",
        str(bed),
        "--evidence",
        f"wgs={calls}",
        "--prefix",
        str(tmp_path / "x"),
    )
    assert result.returncode == 1
    assert (
        "--evidence 'wgs=" in result.stderr and "technology of --technology (wes)" in result.stderr
    )


def test_check_output_order_does_not_depend_on_input_order(tmp_path, ref):
    known = tmp_path / "known.vcf"
    known.write_text(
        CALL_HEADER.replace(
            "#CHROM",
            '##INFO=<ID=KV_GENE,Number=1,Type=String,Description="g">\n'
            '##INFO=<ID=KV_LABEL,Number=1,Type=String,Description="l">\n#CHROM',
        )
        + "chr1\t10\tv1\tA\tC\t.\t.\tKV_GENE=G;KV_LABEL=v1\n"
    )
    calls = tmp_path / "calls.vcf"
    calls.write_text(CALL_HEADER)
    bed = tmp_path / "r.bed"
    bed.write_text("chr1\t0\t100\n")
    # technologies given as wgs, wes: the tables list wes first
    result = run(
        "check",
        "--known",
        str(known),
        "--sample",
        "S1",
        "--technology",
        "wgs",
        "wes",
        "--normalized",
        str(calls),
        str(calls),
        "--compare",
        str(calls),
        str(calls),
        "--regions",
        str(bed),
        str(bed),
        "--prefix",
        str(tmp_path / "x"),
    )
    assert result.returncode == 0, result.stderr
    summary = (tmp_path / "x.known_summary.tsv").read_text().splitlines()
    assert [line.split("\t")[1] for line in summary[1:]] == ["wes", "wgs"]


def test_prepare_ref_mismatch_fails_unless_warn_only(tmp_path, ref):
    fasta, seq = ref
    wrong = alt_of(seq[99])
    tsv = tmp_path / "known.tsv"
    tsv.write_text(
        "chrom\tpos\tref\talt\tlabel\n"
        f"chr1\t100\t{wrong}\t{alt_of(wrong)}\tbad\n"
        f"chr1\t200\t{seq[199]}\t{alt_of(seq[199])}\tgood\n"
    )
    args = ["prepare", "--input", str(tsv), "--fai", f"{fasta}.fai", "--fasta", str(fasta)]
    out = ["--output", str(tmp_path / "k.vcf.gz"), "--report", str(tmp_path / "r.tsv")]
    result = run(*args, *out)
    assert result.returncode == 1
    assert f"REF {wrong} differs from the reference ({seq[99]})" in result.stderr
    assert "ref_mismatches\t1" in (tmp_path / "r.tsv").read_text()
    # check switched off: the mismatching variant is dropped, the run goes on
    result = run(*args, *out, "--warn-only")
    assert result.returncode == 0, result.stderr
    assert [r.id for r in pysam.VariantFile(str(tmp_path / "k.vcf.gz"))] == ["good"]
