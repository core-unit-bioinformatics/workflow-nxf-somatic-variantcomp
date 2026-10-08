"""Unit tests for bin/check_reference.py (decisions on tool outputs)."""

import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "check_reference.py"
FAI = "chr1\t2000\t6\t60\t61\nchr2\t1000\t2046\t60\t61\nchrM\t500\t3070\t60\t61\n"


@pytest.fixture()
def fai(tmp_path):
    path = tmp_path / "ref.fa.fai"
    path.write_text(FAI)
    return path


def vcf_header(contigs, samples=("TUMOR", "NORMAL")):
    lines = ["##fileformat=VCFv4.2"]
    lines += [f"##contig=<ID={c},length={n}>" if n else f"##contig=<ID={c}>" for c, n in contigs]
    lines.append(
        "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]
        )
    )
    return "\n".join(lines) + "\n"


def run(tmp_path, *args):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, args)], capture_output=True, text=True
    )
    report_file = next(tmp_path.glob("out.*.tsv"), None)
    report = {}
    if report_file:
        for line in report_file.read_text().splitlines()[1:]:
            key, value = line.split("\t", 1)
            report.setdefault(key, value)
    return result, report


def check_vcf(tmp_path, fai, header, records, *extra):
    (tmp_path / "header.vcf").write_text(header)
    (tmp_path / "contigs.txt").write_text("".join(f"{c}\n" for c in records))
    return run(tmp_path, "contigs", "--type", "vcf", "--name", "test.vcf", "--fai", fai,
               "--prefix", tmp_path / "out", "--vcf-header", tmp_path / "header.vcf",
               "--record-contigs", tmp_path / "contigs.txt", *extra)  # fmt: skip


def test_vcf_ok(tmp_path, fai):
    result, report = check_vcf(tmp_path, fai, vcf_header([("chr1", 2000), ("chr2", 1000)]), ["chr1", "chr2"],
                               "--samples", "TUMOR,NORMAL")  # fmt: skip
    assert result.returncode == 0, result.stderr
    assert report["status"] == "ok"
    assert (tmp_path / "out.rename_chrs.tsv").read_text() == ""
    assert (tmp_path / "out.keep_contigs.txt").read_text().split() == ["chr1", "chr2"]


def test_vcf_nochr_renamed(tmp_path, fai):
    result, report = check_vcf(tmp_path, fai, vcf_header([]), ["1", "2", "MT"])
    assert result.returncode == 0, result.stderr
    assert report["status"] == "fixed"
    assert (tmp_path / "out.rename_chrs.tsv").read_text().splitlines() == [
        "1\tchr1",
        "2\tchr2",
        "MT\tchrM",
    ]


def test_vcf_mito_length_fixed(tmp_path, fai):
    result, report = check_vcf(tmp_path, fai, vcf_header([("chr1", 2000), ("chrM", 501)]), ["chrM"])
    assert result.returncode == 0, result.stderr
    assert report["header_length_fixes"] == "chrM:501->500"


def test_vcf_nuclear_length_is_build_error(tmp_path, fai):
    result, _ = check_vcf(tmp_path, fai, vcf_header([("chr1", 2500)]), ["chr1"])
    assert result.returncode == 1
    assert "another genome build" in result.stderr


def test_vcf_missing_sample(tmp_path, fai):
    result, _ = check_vcf(
        tmp_path, fai, vcf_header([("chr1", 2000)]), ["chr1"], "--samples", "T1,NORMAL"
    )
    assert result.returncode == 1
    assert "sample(s) T1 not in the VCF (samples in the VCF: TUMOR, NORMAL)" in result.stderr


def test_vcf_unknown_contig_records_dropped(tmp_path, fai):
    result, report = check_vcf(
        tmp_path, fai, vcf_header([("chr1", 2000), ("decoy", 10)]), ["chr1", "decoy", "decoy"]
    )
    assert result.returncode == 0, result.stderr
    assert report["records_on_unknown_contigs"] == "2"
    assert "decoy (2)" in result.stderr


def test_vcf_no_record_on_reference(tmp_path, fai):
    result, _ = check_vcf(tmp_path, fai, vcf_header([]), ["chrUn1", "chrUn2"])
    assert result.returncode == 1
    assert "no record is on a contig of the reference" in result.stderr


def test_bed(tmp_path, fai):
    bed = tmp_path / "in.bed"
    bed.write_text("track name=x\n1\t0\t100\tGENE1\nchr2\t900\t1000\tGENE2\ndecoy\t0\t5\n")
    result, report = run(tmp_path, "contigs", "--type", "bed", "--name", "in.bed", "--fai", fai,
                         "--prefix", tmp_path / "out", "--bed", bed)  # fmt: skip
    assert result.returncode == 0, result.stderr
    assert (
        tmp_path / "out.harmonized.bed"
    ).read_text() == "chr1\t0\t100\tGENE1\nchr2\t900\t1000\tGENE2\n"
    assert report["intervals_on_unknown_contigs"] == "1"


def test_bed_beyond_contig_end(tmp_path, fai):
    bed = tmp_path / "in.bed"
    bed.write_text("chr2\t900\t1200\n")
    result, _ = run(tmp_path, "contigs", "--type", "bed", "--name", "in.bed", "--fai", fai,
                    "--prefix", tmp_path / "out", "--bed", bed)  # fmt: skip
    assert result.returncode == 1
    assert "end beyond their contig (e.g. chr2 at 1200, reference length 1000)" in result.stderr


def test_bam_idxstats(tmp_path, fai):
    idx = tmp_path / "x.idxstats"
    idx.write_text("1\t2000\t10\t0\n2\t1000\t5\t0\nMT\t500\t1\t0\n*\t0\t0\t3\n")
    regions = tmp_path / "regions.bed"
    regions.write_text("chr1\t10\t20\nchrM\t0\t5\n")
    result, report = run(tmp_path, "contigs", "--type", "bam", "--name", "x.bam", "--fai", fai,
                         "--prefix", tmp_path / "out", "--idxstats", idx, "--regions", regions)  # fmt: skip
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "out.regions.bed").read_text() == "1\t10\t20\nMT\t0\t5\n"
    assert report["renamed_contigs"] == "1->chr1,2->chr2,MT->chrM"


def test_bam_other_build(tmp_path, fai):
    idx = tmp_path / "x.idxstats"
    idx.write_text("chr1\t2100\t10\t0\n")
    result, _ = run(tmp_path, "contigs", "--type", "bam", "--name", "x.bam", "--fai", fai,
                    "--prefix", tmp_path / "out", "--idxstats", idx)  # fmt: skip
    assert result.returncode == 1
    assert "another genome build" in result.stderr


@pytest.mark.parametrize(("line", "code"), [
    ("Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:\t1000/2/0/5/3/1/0", 0),
    ("Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:\t100/0/0/0/50/0/0", 1),
])  # fmt: skip
def test_ref_alleles(tmp_path, line, code):
    log = tmp_path / "x.norm.log"
    log.write_text("REF_MISMATCH\tchr1\t10\tA\tC\n" + line + "\n")
    result, report = run(tmp_path, "ref-alleles", "--name", "x.vcf", "--norm-log", log,
                         "--max-ref-mismatch", "0.01", "--prefix", tmp_path / "out")  # fmt: skip
    assert result.returncode == code, result.stderr
    if code:
        assert "50 of 100 REF alleles (50.00%) differ from the reference" in result.stderr
    else:
        assert report["ref_mismatches_removed"] == "3"


def test_ref_alleles_without_summary(tmp_path):
    log = tmp_path / "x.norm.log"
    log.write_text("something else\n")
    result, _ = run(
        tmp_path, "ref-alleles", "--name", "x.vcf", "--norm-log", log, "--prefix", tmp_path / "out"
    )
    assert result.returncode == 1
    assert "no summary line of bcftools norm" in result.stderr


def test_mixed_naming_styles_fail(tmp_path, fai):
    bed = tmp_path / "mixed.bed"
    bed.write_text("1\t10\t20\nchr1\t30\t40\n")
    result = subprocess.run([sys.executable, str(SCRIPT), "contigs", "--type", "bed", "--name", "mixed",
                             "--fai", str(fai), "--prefix", str(tmp_path / "x"), "--bed", str(bed)],
                            capture_output=True, text=True)  # fmt: skip
    assert result.returncode == 1
    assert "contigs 1, chr1 all correspond to the reference contig chr1" in result.stderr


def test_no_contig_on_reference_fails(tmp_path, fai):
    bed = tmp_path / "other.bed"
    bed.write_text("scaffold_1\t10\t20\n")
    result = subprocess.run([sys.executable, str(SCRIPT), "contigs", "--type", "bed", "--name", "other",
                             "--fai", str(fai), "--prefix", str(tmp_path / "x"), "--bed", str(bed)],
                            capture_output=True, text=True)  # fmt: skip
    assert result.returncode == 1
    assert "none of the 1 contigs (scaffold_1) is in the reference" in result.stderr
