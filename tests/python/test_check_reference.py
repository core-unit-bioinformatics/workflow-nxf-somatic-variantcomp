"""Unit tests for bin/check_reference.py (input checks against the reference)."""

import random
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "check_reference.py"
sys.path.insert(0, str(SCRIPT.parent))
import check_reference  # noqa: E402

LENGTHS = {"chr1": 2000, "chr2": 1000, "chrM": 500}


@pytest.fixture()
def ref(tmp_path):
    """Reference FASTA with index; returns (fasta, fai, {contig: sequence})."""
    rng = random.Random(7)
    seqs = {name: "".join(rng.choice("ACGT") for _ in range(n)) for name, n in LENGTHS.items()}
    fasta = tmp_path / "ref.fa"
    fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in seqs.items()))
    pysam.faidx(str(fasta))
    return fasta, Path(f"{fasta}.fai"), seqs


def vcf_text(contigs, records, samples=("TUMOR", "NORMAL")):
    """VCF with the given contig lines [(name, length or None)] and records [(chrom, pos, ref)]."""
    lines = ["##fileformat=VCFv4.2"]
    lines += [f"##contig=<ID={c},length={n}>" if n else f"##contig=<ID={c}>" for c, n in contigs]
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    lines.append(
        "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]
        )
    )
    for chrom, pos, ref_allele in records:
        alt = "A" if ref_allele[0] != "A" else "C"
        lines.append(
            "\t".join(
                [chrom, str(pos), ".", ref_allele, alt, ".", "PASS", ".", "GT"]
                + ["0/1"] * len(samples)
            )
        )
    return "\n".join(lines) + "\n"


def run(tmp_path, *args):
    """Run the script; return (process, report dict or None)."""
    prefix = tmp_path / "out"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), *args, "--name", "x", "--prefix", str(prefix)],
        capture_output=True,
        text=True,
    )
    report_path = Path(f"{prefix}.reference_check.tsv")
    report = None
    if report_path.exists():
        rows = [line.split("\t", 1) for line in report_path.read_text().splitlines()[1:]]
        report = {}
        for key, value in rows:
            report.setdefault(key, value)
    return proc, report


def check_vcf(tmp_path, ref, contigs, records, *extra):
    fasta, fai, _ = ref
    vcf = tmp_path / "in.vcf"
    vcf.write_text(vcf_text(contigs, records))
    return run(tmp_path, "vcf", "--vcf", str(vcf), "--fai", str(fai), "--fasta", str(fasta), *extra)


def real_ref(ref, chrom, pos, n=1):
    """REF allele of the test reference at a 1-based position."""
    return ref[2][chrom][pos - 1 : pos - 1 + n]


def test_vcf_ok(tmp_path, ref):
    proc, report = check_vcf(
        tmp_path,
        ref,
        [("chr1", 2000), ("chr2", 1000)],
        [("chr1", 10, real_ref(ref, "chr1", 10)), ("chr2", 5, real_ref(ref, "chr2", 5, 2))],
        "--samples",
        "TUMOR,NORMAL",
    )
    assert proc.returncode == 0, proc.stderr
    assert report["status"] == "ok"
    assert report["records"] == "2" and report["ref_mismatches"] == "0"


def test_vcf_nochr_renamed(tmp_path, ref):
    proc, report = check_vcf(tmp_path, ref, [("1", 2000)], [("1", 10, real_ref(ref, "chr1", 10))])
    assert proc.returncode == 0, proc.stderr
    assert report["status"] == "fixed" and report["renamed_contigs"] == "1->chr1"
    assert (tmp_path / "out.rename_chrs.tsv").read_text() == "1\tchr1\n"


def test_vcf_mito_dropped_even_with_other_length(tmp_path, ref):
    # chrM with another length is not checked: the contig is not compared
    proc, report = check_vcf(
        tmp_path,
        ref,
        [("chr1", 2000), ("MT", 16571)],
        [("chr1", 10, real_ref(ref, "chr1", 10)), ("MT", 5, "A")],
    )
    assert proc.returncode == 0, proc.stderr
    assert report["dropped_mito"] == "MT" and report["records_on_mito"] == "1"
    assert report["length_mismatches"] == "none"


def test_vcf_mito_kept_is_length_checked(tmp_path, ref):
    proc, report = check_vcf(
        tmp_path,
        ref,
        [("chr1", 2000), ("chrM", 16571)],
        [("chr1", 10, real_ref(ref, "chr1", 10))],
        "--keep-mito",
    )
    assert proc.returncode == 1
    assert "another genome build" in proc.stderr and report["status"] == "error"


def test_vcf_nuclear_length_is_build_error(tmp_path, ref):
    proc, _ = check_vcf(tmp_path, ref, [("chr1", 2001)], [("chr1", 10, real_ref(ref, "chr1", 10))])
    assert proc.returncode == 1
    assert "chr1:2001 (reference chr1:2000)" in proc.stderr


def test_warn_only_turns_errors_into_warnings(tmp_path, ref):
    proc, report = check_vcf(
        tmp_path, ref, [("chr1", 2001)], [("chr1", 10, real_ref(ref, "chr1", 10))], "--warn-only"
    )
    assert proc.returncode == 0, proc.stderr
    assert "(check switched off)" in proc.stderr and report["status"] == "fixed"


def test_vcf_missing_sample(tmp_path, ref):
    proc, _ = check_vcf(
        tmp_path,
        ref,
        [("chr1", 2000)],
        [("chr1", 10, real_ref(ref, "chr1", 10))],
        "--samples",
        "TUMOR1",
    )
    assert proc.returncode == 1
    assert "sample(s) TUMOR1 not in the VCF (samples in the VCF: TUMOR, NORMAL)" in proc.stderr


def test_vcf_unknown_contig_records_dropped(tmp_path, ref):
    proc, report = check_vcf(
        tmp_path,
        ref,
        [("chr1", 2000)],
        [("chr1", 10, real_ref(ref, "chr1", 10)), ("chrUn_x", 3, "A")],
    )
    assert proc.returncode == 0, proc.stderr
    assert report["records_on_unknown_contigs"] == "1" and report["unknown_contigs"] == "chrUn_x"


def test_vcf_no_record_on_reference(tmp_path, ref):
    proc, _ = check_vcf(tmp_path, ref, [], [("scaffold9", 3, "A")])
    assert proc.returncode == 1
    assert "none of the 1 contigs (scaffold9) is in the reference" in proc.stderr


def test_vcf_ref_mismatch_rate(tmp_path, ref):
    wrong = "A" if real_ref(ref, "chr1", 20) != "A" else "C"
    records = [("chr1", 10, real_ref(ref, "chr1", 10)), ("chr1", 20, wrong)]
    # 1 of 2 mismatching: error with the default threshold of 1 %
    proc, report = check_vcf(tmp_path, ref, [("chr1", 2000)], records)
    assert proc.returncode == 1
    assert report["ref_mismatches"] == "1" and "differ from the reference" in proc.stderr
    # accepted with a threshold of 60 %
    proc, _ = check_vcf(tmp_path, ref, [("chr1", 2000)], records, "--max-ref-mismatch", "0.6")
    assert proc.returncode == 0, proc.stderr


def test_bed_harmonized(tmp_path, ref):
    bed = tmp_path / "t.bed"
    bed.write_text("1\t100\t200\tA\nchrUn_x\t1\t5\tB\nMT\t10\t20\tC\n")
    proc, report = run(tmp_path, "bed", "--bed", str(bed), "--fai", str(ref[1]))
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "out.harmonized.bed").read_text() == "chr1\t100\t200\tA\n"
    assert report["intervals_dropped"] == "2" and report["dropped_mito"] == "MT"


def test_bed_beyond_contig_end(tmp_path, ref):
    bed = tmp_path / "t.bed"
    bed.write_text("chr1\t100\t200\nchr1\t1900\t2100\n")
    proc, _ = run(tmp_path, "bed", "--bed", str(bed), "--fai", str(ref[1]))
    assert proc.returncode == 1
    assert "line 2: chr1 at 2100" in proc.stderr and "another genome build" in proc.stderr


@pytest.mark.parametrize(
    ("line", "message"),
    [
        ("chr1 0 10", "space-separated? BED needs tabs"),
        ("chr1\ta\t10", "are not integers"),
        ("chr1\t20\t10", "invalid interval 20-10"),
    ],
)
def test_read_bed_errors(tmp_path, line, message):
    bed = tmp_path / "x.bed"
    bed.write_text(line + "\n")
    with pytest.raises(SystemExit) as err:
        check_reference.read_bed(bed)
    assert "line 1" in str(err.value) and message in str(err.value)


def test_mixed_naming_styles_fail(tmp_path, ref):
    bed = tmp_path / "mixed.bed"
    bed.write_text("1\t10\t20\nchr1\t30\t40\n")
    proc, _ = run(tmp_path, "bed", "--bed", str(bed), "--fai", str(ref[1]))
    assert proc.returncode == 1
    assert "contigs 1, chr1 all correspond to the reference contig chr1" in proc.stderr


def make_bam(path, contigs):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": c, "LN": n} for c, n in contigs]}
    with pysam.AlignmentFile(str(path), "wb", header=header):
        pass
    return path


def test_bam_renamed(tmp_path, ref):
    bam = make_bam(tmp_path / "x.bam", [("1", 2000), ("2", 1000), ("MT", 16569)])
    proc, report = run(tmp_path, "bam", "--alignment", str(bam), "--fai", str(ref[1]))
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "out.rename_chrs.tsv").read_text() == "1\tchr1\n2\tchr2\n"
    assert report["dropped_mito"] == "MT"


def test_bam_other_build(tmp_path, ref):
    bam = make_bam(tmp_path / "x.bam", [("chr1", 2500)])
    proc, _ = run(tmp_path, "bam", "--alignment", str(bam), "--fai", str(ref[1]))
    assert proc.returncode == 1 and "another genome build" in proc.stderr


def test_rename_bed(tmp_path):
    rename = tmp_path / "r.tsv"
    rename.write_text("1\tchr1\n")
    bed = tmp_path / "common.bed"
    bed.write_text("chr1\t0\t10\nchr2\t5\t9\n")
    out = tmp_path / "out.bed"
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "rename-bed",
            "--bed",
            str(bed),
            "--rename",
            str(rename),
            "--output",
            str(out),
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert out.read_text() == "1\t0\t10\nchr2\t5\t9\n"
