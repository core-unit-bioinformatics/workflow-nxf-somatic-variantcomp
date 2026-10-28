"""Unit tests for bin/read_evidence.py: allele-specific read counts with known answers.

The reference has an AG repeat; the variants are a 2-bp deletion of AG (left-aligned),
an insertion of AG and an SNV. The reads are built so that every count is known: the
deletion placed at the left or the right end of the repeat by the aligner, a read with
the deletion in its soft-clipped part, a read with another indel (4-bp deletion), reads
that don't cover the window, low mapping quality, a duplicate, a low base quality.
"""

import subprocess
import sys
from pathlib import Path

import pysam

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "read_evidence.py"
sys.path.insert(0, str(SCRIPT.parent))
import read_evidence  # noqa: E402

# 0-based layout: 40 unique bases, AGAGAGAG at 40-47, 40 unique bases
LEFT = "GATCCTAGCTTGACCGTATGCCAGTTCAGGCTAACGTTCC"
RIGHT = "TCGGACTTACCGATGCATTGCCAAGTCTGATCGGTACCAT"
REPEAT = "AGAGAGAG"
REF = LEFT + REPEAT + RIGHT
# deletion of one AG unit, left-aligned: anchor = last base before the repeat (pos 40, 1-based)
DEL = ("chr1", 40, REF[39:42], REF[39])  # C AG -> C
INS = ("chr1", 40, REF[39], REF[39] + "AG")  # C -> C AG
SNV = ("chr1", 70, REF[69], "A" if REF[69] != "A" else "C")


def write_reference(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + REF + "\n")
    pysam.faidx(str(fasta))
    return str(fasta)


def read(name, seq, start, cigar, mapq=60, dup=False, qual=30):
    """An aligned read (0-based start, CIGAR as a pysam tuple list)."""
    r = pysam.AlignedSegment()
    r.query_name = name
    r.query_sequence = seq
    r.reference_id = 0
    r.reference_start = start
    r.cigartuples = cigar
    r.mapping_quality = mapq
    r.query_qualities = pysam.qualitystring_to_array(chr(qual + 33) * len(seq))
    r.is_duplicate = dup
    return r


def write_bam(tmp_path, reads):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(REF)}]}
    path = tmp_path / "unsorted.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for r in reads:
            out.write(r)
    sorted_path = tmp_path / "reads.bam"
    pysam.sort("-o", str(sorted_path), str(path))
    pysam.index(str(sorted_path))
    return str(sorted_path)


M, C_INS, DEL_OP, S = 0, 1, 2, 4  # CIGAR operations
DEL_HAP = LEFT + REPEAT[2:] + RIGHT  # haplotype with one AG unit less
INS_HAP = LEFT + "AG" + REPEAT + RIGHT


def deletion_reads():
    reads = [
        # 1. deletion placed at the left end of the repeat (positions 40-41)
        read("del_left", DEL_HAP, 0, [(M, 40), (DEL_OP, 2), (M, 46)]),
        # 2. the same deletion placed by the aligner at the right end (positions 46-47)
        read("del_right", DEL_HAP, 0, [(M, 46), (DEL_OP, 2), (M, 40)]),
        # 3. aligned part ends before the repeat; the deletion haplotype is in the soft clip
        read("del_clipped", DEL_HAP[:70], 0, [(M, 38), (S, 32)]),
        # 4. reference reads
        read("ref_1", REF, 0, [(M, len(REF))]),
        read("ref_2", REF[10:80], 10, [(M, 70)]),
        # 5. another indel: 4-bp deletion (two AG units)
        read("other_del4", LEFT + REPEAT[4:] + RIGHT, 0, [(M, 40), (DEL_OP, 4), (M, 44)]),
        # 6. doesn't cover the window (ends inside the repeat): ignored
        read("short", REF[:44], 0, [(M, 44)]),
        # 7. low mapping quality: ignored
        read("low_mq", DEL_HAP, 0, [(M, 40), (DEL_OP, 2), (M, 46)], mapq=5),
        # 8. duplicate: ignored unless --keep-duplicates
        read("dup", DEL_HAP, 0, [(M, 40), (DEL_OP, 2), (M, 46)], dup=True),
    ]
    return reads


def run(tmp_path, reads, sites, *extra):
    fasta = write_reference(tmp_path)
    bam = write_bam(tmp_path, reads)
    sites_file = tmp_path / "sites.tsv"
    sites_file.write_text("".join(f"{c}\t{p}\t{r}\t{a}\n" for c, p, r, a in sites))
    out = tmp_path / "ev.vcf.gz"
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--alignment",
            bam,
            "--fasta",
            fasta,
            "--sites",
            str(sites_file),
            "--sample",
            "S1_tumor",
            "--output",
            str(out),
            *extra,
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr
    counts = {}
    with pysam.VariantFile(str(out)) as vcf:
        for rec in vcf:
            fmt = rec.samples["S1_tumor"]
            counts[(rec.chrom, rec.pos, rec.ref, rec.alts[0])] = (
                fmt["DP"],
                tuple(fmt["AD"]),
                fmt["OT"],
            )
    return counts


def test_deletion_counted_wherever_the_aligner_placed_it(tmp_path):
    counts = run(tmp_path, deletion_reads(), [DEL])
    # alt: left, right placement and soft-clipped; ref: 2; other: the 4-bp deletion
    assert counts[DEL] == (6, (2, 3), 1)


def test_duplicates_kept_on_request(tmp_path):
    counts = run(tmp_path, deletion_reads(), [DEL], "--keep-duplicates")
    assert counts[DEL][1] == (2, 4)


def test_insertion_and_snv(tmp_path):
    snv_pos0 = SNV[1] - 1
    snv_hap = REF[:snv_pos0] + SNV[3] + REF[snv_pos0 + 1 :]
    reads = [
        # insertion of AG placed at the right end of the repeat by the aligner
        read("ins_right", INS_HAP, 0, [(M, 48), (C_INS, 2), (M, 40)]),
        read("ins_left", INS_HAP, 0, [(M, 40), (C_INS, 2), (M, 48)]),
        read("ref", REF, 0, [(M, len(REF))]),
        # SNV: one alt read, one alt read with low base quality (not counted)
        read("snv_alt", snv_hap, 0, [(M, len(REF))]),
        read("snv_lowq", snv_hap, 0, [(M, len(REF))], qual=5),
    ]
    counts = run(tmp_path, reads, [INS, SNV])
    assert counts[INS] == (5, (3, 2), 0)
    # SNV: ins reads carry the reference base at the SNV (their alignment is shifted by
    # the insertion before it, but the SNV lies after it in the read too)
    assert counts[SNV][1][1] == 1
    assert counts[SNV][0] == 4


def test_ref_mismatch_is_an_error(tmp_path):
    fasta = write_reference(tmp_path)
    bam = write_bam(tmp_path, [read("ref", REF, 0, [(M, len(REF))])])
    sites = tmp_path / "s.tsv"
    sites.write_text("chr1\t40\tGG\tG\n")
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--alignment",
            bam,
            "--fasta",
            fasta,
            "--sites",
            str(sites),
            "--sample",
            "x",
            "--output",
            str(tmp_path / "o.vcf.gz"),
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 1 and "REF differs from the reference" in proc.stderr


def test_shift_range_of_a_deletion_in_a_repeat():
    # deleting AG from AGAGAGAG at 40-41 can be placed anywhere up to 46-47
    assert read_evidence.shift_range(REF, 40, 42, "") == (40, 48)
