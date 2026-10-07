"""Unit tests for bin/check_reference.py."""

import random
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "check_reference.py"
CONTIGS = {"chr1": 2000, "chr2": 1000, "chrM": 500}


@pytest.fixture()
def reference(tmp_path):
    rng = random.Random(1)
    seqs = {c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in CONTIGS.items()}
    fasta = tmp_path / "ref.fa"
    fasta.write_text("".join(f">{c}\n{s}\n" for c, s in seqs.items()))
    pysam.faidx(str(fasta))
    return fasta, seqs


def write_vcf(path, records, contig_lines):
    header = "##fileformat=VCFv4.2\n"
    header += "".join(
        f"##contig=<ID={c},length={n}>\n" if n else f"##contig=<ID={c}>\n" for c, n in contig_lines
    )
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    body = "".join(f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\n" for c, p, r, a in records)
    path.write_text(header + body)
    return path


def run_check(tmp_path, input_path, kind, reference, *extra):
    fasta, _ = reference
    prefix = tmp_path / "out"
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--input",
        str(input_path),
        "--type",
        kind,
        "--fai",
        f"{fasta}.fai",
        "--fasta",
        str(fasta),
        "--prefix",
        str(prefix),
        *extra,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    report = dict(
        line.split("\t", 1)
        for line in (tmp_path / "out.reference_check.tsv").read_text().splitlines()[1:]
        if not line.startswith(("warning", "error"))
    )
    renames = (tmp_path / "out.rename_chrs.tsv").read_text()
    keep = (tmp_path / "out.keep_contigs.txt").read_text().split()
    return result, report, renames, keep


def snv(seqs, chrom, pos, name=None):
    ref = seqs[chrom][pos - 1]
    return (name or chrom, pos, ref, "A" if ref != "A" else "C")


def test_matching_vcf_is_ok(tmp_path, reference):
    _, seqs = reference
    vcf = write_vcf(
        tmp_path / "in.vcf", [snv(seqs, "chr1", 10), snv(seqs, "chr2", 20)], list(CONTIGS.items())
    )
    result, report, renames, keep = run_check(tmp_path, vcf, "vcf", reference)
    assert result.returncode == 0
    assert report["status"] == "ok"
    assert report["ref_mismatches"] == "0"
    assert renames == ""
    # all reference contigs named in header or records
    assert keep == ["chr1", "chr2", "chrM"]


def test_nochr_names_are_renamed(tmp_path, reference):
    _, seqs = reference
    records = [snv(seqs, "chr1", 10, "1"), snv(seqs, "chrM", 5, "MT")]
    vcf = write_vcf(tmp_path / "in.vcf", records, [("1", 2000), ("2", 1000), ("MT", 500)])
    result, report, renames, keep = run_check(tmp_path, vcf, "vcf", reference)
    assert result.returncode == 0
    assert report["status"] == "fixed"
    assert renames.splitlines() == ["1\tchr1", "2\tchr2", "MT\tchrM"]
    assert keep == ["chr1", "chr2", "chrM"]


def test_mito_header_length_is_fixed(tmp_path, reference):
    _, seqs = reference
    vcf = write_vcf(
        tmp_path / "in.vcf", [snv(seqs, "chrM", 5)], [("chr1", 2000), ("chr2", 1000), ("chrM", 501)]
    )
    result, report, _, _ = run_check(tmp_path, vcf, "vcf", reference)
    assert result.returncode == 0
    assert report["status"] == "fixed"
    assert report["header_length_fixes"] == "chrM:501->500"


def test_nuclear_length_mismatch_is_build_error(tmp_path, reference):
    _, seqs = reference
    vcf = write_vcf(tmp_path / "in.vcf", [snv(seqs, "chr1", 10)], [("chr1", 2500), ("chr2", 1000)])
    result, report, _, _ = run_check(tmp_path, vcf, "vcf", reference)
    assert result.returncode == 1
    assert report["status"] == "error"
    assert "different genome build" in result.stderr


def test_ref_mismatch_is_error(tmp_path, reference):
    _, seqs = reference
    records = []
    for pos in range(10, 30):
        c, p, ref, alt = snv(seqs, "chr1", pos)
        records.append((c, p, alt, ref))  # REF and ALT swapped -> REF mismatch
    vcf = write_vcf(tmp_path / "in.vcf", records, list(CONTIGS.items()))
    result, report, _, _ = run_check(tmp_path, vcf, "vcf", reference)
    assert result.returncode == 1
    assert report["ref_mismatches"] == "20"
    assert "do not match the reference" in result.stderr


def test_unknown_contig_records_are_dropped(tmp_path, reference):
    _, seqs = reference
    records = [snv(seqs, "chr1", 10), ("chrUn_decoy", 5, "A", "C")]
    vcf = write_vcf(tmp_path / "in.vcf", records, [*CONTIGS.items(), ("chrUn_decoy", 100)])
    result, report, _, keep = run_check(tmp_path, vcf, "vcf", reference)
    assert result.returncode == 0
    assert report["records_on_unknown_contigs"] == "1"
    assert "chrUn_decoy" not in keep


def test_bed_beyond_contig_end_is_error(tmp_path, reference):
    bed = tmp_path / "in.bed"
    bed.write_text("1\t0\t100\nchr2\t900\t1200\n")
    result, report, renames, _ = run_check(tmp_path, bed, "bed", reference)
    assert result.returncode == 1
    assert renames == "1\tchr1\n"
    assert report["intervals_beyond_contig_end"] == "1"
    assert (tmp_path / "out.harmonized.bed").read_text() == "chr1\t0\t100\nchr2\t900\t1200\n"


def test_bed_ok(tmp_path, reference):
    bed = tmp_path / "in.bed"
    bed.write_text("track name=x\nchr1\t0\t100\tGENE1\nchr2\t900\t1000\tGENE2\n")
    result, report, _, keep = run_check(tmp_path, bed, "bed", reference)
    assert (
        tmp_path / "out.harmonized.bed"
    ).read_text() == "chr1\t0\t100\tGENE1\nchr2\t900\t1000\tGENE2\n"
    assert result.returncode == 0
    assert report["status"] == "ok"
    assert keep == ["chr1", "chr2"]


def test_missing_sample_is_error(tmp_path, reference):
    _, seqs = reference
    vcf = tmp_path / "in.vcf"
    write_vcf(vcf, [], list(CONTIGS.items()))
    text = vcf.read_text().replace("\tINFO\n", "\tINFO\tFORMAT\tTUMOR\tNORMAL\n")
    vcf.write_text(text)
    result, report, _, _ = run_check(tmp_path, vcf, "vcf", reference, "--samples", "T1,NORMAL")
    assert result.returncode == 1
    assert report["samples"] == "TUMOR,NORMAL"
    assert "sample(s) T1 not in the VCF" in result.stderr


def write_bam(path, contigs):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": c, "LN": n} for c, n in contigs]}
    with pysam.AlignmentFile(str(path), "wb", header=header):
        pass
    return path


def test_bam_nochr_is_renamed(tmp_path, reference):
    bam = write_bam(tmp_path / "in.bam", [("1", 2000), ("2", 1000), ("MT", 500), ("decoy", 50)])
    regions = tmp_path / "regions.bed"
    regions.write_text("chr1\t10\t20\nchrM\t0\t5\nchrX\t0\t5\n")
    result, report, renames, keep = run_check(
        tmp_path, bam, "bam", reference, "--regions", str(regions)
    )
    assert (tmp_path / "out.regions.bed").read_text() == "1\t10\t20\nMT\t0\t5\n"
    assert result.returncode == 0, result.stderr
    assert report["status"] == "fixed"
    assert renames.splitlines() == ["1\tchr1", "2\tchr2", "MT\tchrM"]
    assert keep == ["chr1", "chr2", "chrM"]


def test_bam_other_build_is_error(tmp_path, reference):
    bam = write_bam(tmp_path / "in.bam", [("chr1", 2100), ("chr2", 1000)])
    result, report, _, _ = run_check(tmp_path, bam, "bam", reference)
    assert result.returncode == 1
    assert "different genome build" in result.stderr
