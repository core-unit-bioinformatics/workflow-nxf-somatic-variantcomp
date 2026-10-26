"""Unit tests for bin/annotate_variants.py (variant table <-> tool inputs/results)."""

import csv
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "annotate_variants.py"
sys.path.insert(0, str(SCRIPT.parent))
import annotate_variants  # noqa: E402

TABLE = "chrom\tpos\tref\talt\nchr1\t100\tA\tC\nchr1\t200\tGTT\tG\nchr1\t100\tA\tC\n"
VEP_HEADER = (
    "# [1]CHROM\t[2]POS\t[3]REF\t[4]ALT\t[5]Consequence\t[6]SYMBOL\t[7]IMPACT\t"
    "[8]VARIANT_CLASS\t[9]BIOTYPE\n"
)


def run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)


def write(path, text):
    path.write_text(text)
    return str(path)


def rows_of(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.mark.parametrize(
    ("consequence", "group"),
    [
        ("missense_variant", "non_synonymous"),
        ("synonymous_variant", "synonymous"),
        ("splice_region_variant&synonymous_variant", "synonymous"),
        ("frameshift_variant&splice_region_variant", "non_synonymous"),
        ("intron_variant", "other"),
        ("", "other"),
    ],
)
def test_consequence_group(consequence, group):
    assert annotate_variants.consequence_group(consequence) == group


def test_to_bed_covers_the_ref_allele(tmp_path):
    out = tmp_path / "v.bed"
    proc = run("to-bed", "--table", write(tmp_path / "v.tsv", TABLE), "--output", str(out))
    assert proc.returncode == 0, proc.stderr
    assert out.read_text() == "chr1\t99\t100\t1\nchr1\t199\t202\t2\nchr1\t99\t100\t3\n"


def test_to_vcf_distinct_and_sorted(tmp_path):
    fai = write(tmp_path / "ref.fa.fai", "chr1\t1000\t6\t60\t61\n")
    out = tmp_path / "v.vcf.gz"
    proc = run(
        "to-vcf", "--table", write(tmp_path / "v.tsv", TABLE), "--fai", fai, "--output", str(out)
    )
    assert proc.returncode == 0, proc.stderr
    recs = [(r.pos, r.ref, r.alts[0]) for r in pysam.VariantFile(str(out))]
    assert recs == [(100, "A", "C"), (200, "GTT", "G")]


def test_repeats_to_bed(tmp_path):
    scan = write(
        tmp_path / "scan.list",
        "chromosome\tlocation\trepeat_unit_length\trepeat_times\trepeat_unit_bases\tleft_flank_bases\n"
        "chr1\t4999\t1\t12\tA\tTTAGG\nchr2\t2999\t2\t8\tCA\tTGAGG\n",
    )
    out = tmp_path / "r.bed"
    proc = run("repeats-to-bed", "--list", scan, "--output", str(out))
    assert proc.returncode == 0, proc.stderr
    assert out.read_text() == "chr1\t4999\t5011\tA\t12\nchr2\t2999\t3015\tCA\t8\n"


def test_merge_adds_vep_repeat_and_stratum_columns(tmp_path):
    table = write(tmp_path / "v.tsv", TABLE)
    vep = write(
        tmp_path / "vep.tsv",
        VEP_HEADER
        + "chr1\t100\tA\tC\tmissense_variant\tKRAS\tMODERATE\tSNV\tprotein_coding\n"
        + "chr1\t200\tGTT\tG\tintron_variant\t\tMODIFIER\tdeletion\t\n",
    )
    # bedtools output is in sorted order, not table order: rows are matched by name
    repeats = write(
        tmp_path / "rep.bed",
        "chr1\t99\t100\t1\t.\t.\nchr1\t199\t202\t2\tT\t9\nchr1\t99\t100\t3\t.\t.\n",
    )
    hp = write(
        tmp_path / "hp.bed", "chr1\t99\t100\t1\t0\nchr1\t199\t202\t2\t2\nchr1\t99\t100\t3\t0\n"
    )
    out = tmp_path / "annotated.tsv"
    proc = run(
        "merge",
        "--table",
        table,
        "--vep",
        vep,
        "--repeats",
        repeats,
        "--stratum",
        f"homopolymer={hp}",
        "--output",
        str(out),
    )
    assert proc.returncode == 0, proc.stderr
    rows = rows_of(out)
    assert [
        (r["gene"], r["consequence_group"], r["in_repeat"], r["repeat_unit"], r["homopolymer"])
        for r in rows
    ] == [
        ("KRAS", "non_synonymous", "no", ".", "no"),
        (".", "other", "yes", "T", "yes"),
        ("KRAS", "non_synonymous", "no", ".", "no"),
    ]


def test_merge_optional_vep_fields(tmp_path):
    table = write(tmp_path / "v.tsv", TABLE)
    header = VEP_HEADER.rstrip("\n") + "\t[10]Feature\t[11]HGVSc\t[12]HGVSp\t"
    vep = write(
        tmp_path / "vep.tsv",
        header
        + "[13]Existing_variation\t[14]CLIN_SIG\n"
        + "chr1\t100\tA\tC\tmissense_variant\tKRAS\tMODERATE\tSNV\tprotein_coding\t"
        + "ENST1.4\tENST1.4:c.35G>T\tENSP1.1:p.Gly12Val\trs121913529&COSV55497369\t"
        + "pathogenic&likely_pathogenic\n"
        + "chr1\t200\tGTT\tG\tintron_variant\t\tMODIFIER\tdeletion\t\t.\t.\t.\t.\t.\n",
    )
    out = tmp_path / "annotated.tsv"
    proc = run("merge", "--table", table, "--vep", vep, "--output", str(out))
    assert proc.returncode == 0, proc.stderr
    first, second, _ = rows_of(out)
    assert [first[c] for c in annotate_variants.OPTIONAL_COLUMNS] == [
        "ENST1.4",
        "c.35G>T",
        "p.Gly12Val",
        "rs121913529",
        "COSV55497369",
        "pathogenic,likely_pathogenic",
    ]
    assert {second[c] for c in annotate_variants.OPTIONAL_COLUMNS} == {"."}


def test_merge_without_optional_fields_writes_dots(tmp_path):
    table = write(tmp_path / "v.tsv", "chrom\tpos\tref\talt\nchr1\t100\tA\tC\n")
    vep = write(
        tmp_path / "vep.tsv",
        VEP_HEADER + "chr1\t100\tA\tC\tmissense_variant\tKRAS\tMODERATE\tSNV\tx\n",
    )
    out = tmp_path / "annotated.tsv"
    assert run("merge", "--table", table, "--vep", vep, "--output", str(out)).returncode == 0
    assert {rows_of(out)[0][c] for c in annotate_variants.OPTIONAL_COLUMNS} == {"."}


def test_merge_optional_fields_out_of_order_is_an_error(tmp_path):
    table = write(tmp_path / "v.tsv", TABLE)
    vep = write(tmp_path / "vep.tsv", VEP_HEADER.rstrip("\n") + "\t[10]HGVSp\t[11]HGVSc\n")
    proc = run("merge", "--table", table, "--vep", vep, "--output", str(tmp_path / "a.tsv"))
    assert proc.returncode == 1 and "in this order" in proc.stderr


def test_merge_missing_vep_annotation_is_an_error(tmp_path):
    table = write(tmp_path / "v.tsv", TABLE)
    vep = write(
        tmp_path / "vep.tsv",
        VEP_HEADER + "chr1\t100\tA\tC\tmissense_variant\tKRAS\tMODERATE\tSNV\tx\n",
    )
    proc = run("merge", "--table", table, "--vep", vep, "--output", str(tmp_path / "a.tsv"))
    assert proc.returncode == 1
    assert "no VEP annotation for line 3" in proc.stderr


def test_merge_wrong_vep_columns_is_an_error(tmp_path):
    table = write(tmp_path / "v.tsv", TABLE)
    vep = write(tmp_path / "vep.tsv", "# [1]CHROM\t[2]POS\n")
    proc = run("merge", "--table", table, "--vep", vep, "--output", str(tmp_path / "a.tsv"))
    assert proc.returncode == 1 and "expected CHROM, POS, REF, ALT" in proc.stderr


def test_merge_missing_bedtools_row_is_an_error(tmp_path):
    table = write(tmp_path / "v.tsv", TABLE)
    hp = write(tmp_path / "hp.bed", "chr1\t99\t100\t1\t0\n")
    proc = run(
        "merge", "--table", table, "--stratum", f"hp={hp}", "--output", str(tmp_path / "a.tsv")
    )
    assert proc.returncode == 1
    assert "2 variant(s) missing (rows 2, 3)" in proc.stderr


def test_missing_table_column_is_an_error(tmp_path):
    proc = run(
        "to-bed",
        "--table",
        write(tmp_path / "v.tsv", "chrom\tposition\nchr1\t1\n"),
        "--output",
        str(tmp_path / "v.bed"),
    )
    assert proc.returncode == 1
    assert "missing column(s) pos, ref, alt" in proc.stderr
