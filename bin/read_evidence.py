#!/usr/bin/env python3
"""Allele-specific read counts of one alignment at given variants (read evidence).

For every variant (CHROM, POS, REF, ALT; normalized: split, left-aligned, as in the
compared VCFs) the script counts the reads of one tumor or normal alignment that support
the ALT allele, the REF allele or something else. classify_discordance.py and
known_variants.py read these counts by (CHROM, POS, REF, ALT).

INPUT
  --alignment   BAM or CRAM (indexed); a CRAM needs --fasta
  --fasta       reference FASTA (indexed, .fai)
  --sites       one or more tables CHROM<TAB>POS<TAB>REF<TAB>ALT (1-based POS, no header;
                duplicates and several ALTs at one position are allowed)
  --sample      sample name of the output (e.g. S1_wes_tumor)
OUTPUT
  --output      bgzipped VCF (+ .tbi), one record per distinct variant, sorted by the
                contig order of the FASTA index: FORMAT/DP (int: reads counted at the
                variant = REF + ALT + other), FORMAT/AD (int, int: REF reads, ALT reads),
                FORMAT/OT (int: reads with another allele)

How a read is counted (after the read filters: unmapped, secondary, supplementary,
QC-failed and, unless --keep-duplicates, duplicate reads are skipped, and reads with
mapping quality < --min-mapq; at most --max-depth reads per variant, in file order):
  SNV / MNV (REF and ALT of equal length): the read must align a base to every REF
      position with base quality >= --min-baseq; its bases are compared with ALT and REF
      (ALT, REF or other). Reads with a deletion or a low-quality base there are not
      counted.
  indel / complex (lengths differ): the window around the variant contains every
      position the indel can be shifted to in a repeat, plus --flank bases on each side.
      The ALT haplotype is this window with the variant applied, the REF haplotype the
      reference window. A read supports ALT if its full sequence (soft-clipped bases
      included) contains the ALT haplotype, REF if it contains the REF haplotype; a read
      whose alignment spans the window but contains neither (another indel, or a base
      error inside the window) counts as other. The alignment itself (CIGAR, placement of
      the indel, soft clips) doesn't matter, so the count doesn't depend on where or
      whether the aligner placed the indel. Reads are taken within --max-read-length bases
      of the window (soft-clipped reads align next to the event).

Why a script: bcftools mpileup + call, used before, discovers indel alleles itself and
reports them in its own representation. In repeats it reports other lengths than the
caller (e.g. CAG>C and CAGAGAG>C for a called CAGAG>C), and it misses long deletions (a
37-bp deletion called in both technologies had no indel allele in the pileup), so the
alt reads of called indels were counted as 0. No tool counts reads for given indel
alleles independently of their placement (bam-readcount and samtools mpileup report
the indels as aligned).
"""

import argparse
import sys
from pathlib import Path

import pysam

# readers shared with the reference check (clear errors for unreadable files)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_reference import open_alignment, read_fai  # noqa: E402

# an indel is not shifted further than this many bases (very long repeats)
MAX_SHIFT = 500


def read_sites(paths):
    """Distinct variants of the site tables: {(chrom, pos, ref, alt)}; strict format."""
    sites = set()
    for path in paths:
        with open(path) as fh:
            for line_no, line in enumerate(fh, start=1):
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                # four columns: CHROM, POS (1-based integer), REF, ALT (bases only)
                if len(fields) != 4 or not fields[1].isdigit():
                    sys.exit(
                        f"ERROR: {path}: line {line_no}: expected CHROM<TAB>POS<TAB>REF<TAB>ALT "
                        f"with an integer POS, got '{line.rstrip()}'"
                    )
                chrom, pos, ref, alt = (
                    fields[0],
                    int(fields[1]),
                    fields[2].upper(),
                    fields[3].upper(),
                )
                # symbolic or overlapping alleles have no read evidence: skipped, counted
                if not set(ref + alt) <= set("ACGTN") or alt in ("*", "."):
                    continue
                sites.add((chrom, pos, ref, alt))
    return sites


def shift_range(seq, start, end, inserted):
    """(leftmost, rightmost) 0-based reference positions spanned by all placements of an
    indel that is equivalent in a repeat.

    seq: reference sequence of the contig; for a deletion of seq[start:end] (inserted ''),
    for an insertion of `inserted` before position `start` (start == end). The event is
    shifted left and right while the shifted event gives the same haplotype.
    """
    if inserted:
        # insertion: shift right while the next reference base equals the first inserted
        # base (the inserted sequence rotates), left likewise with the last base
        ins, right = inserted, start
        while right < len(seq) and seq[right] == ins[0] and right - start < MAX_SHIFT:
            ins, right = ins[1:] + ins[0], right + 1
        ins, left = inserted, start
        while left > 0 and seq[left - 1] == ins[-1] and start - left < MAX_SHIFT:
            ins, left = ins[-1] + ins[:-1], left - 1
        return left, right
    # deletion of seq[start:end]: shift right while the base after the deletion equals
    # its first base, left while the base before it equals its last base
    length, left, right = end - start, start, end
    while right < len(seq) and seq[right] == seq[right - length] and right - end < MAX_SHIFT:
        right += 1
    while left > 0 and seq[left - 1] == seq[left - 1 + length] and start - left < MAX_SHIFT:
        left -= 1
    return left, right


def indel_window(seq, start0, ref, alt, flank):
    """Window (ws, we) (0-based, end exclusive) and the REF and ALT haplotype sequences in
    it, for a variant with REF/ALT of different length at 0-based start0, with `flank`
    bases on each side of all placements of the event."""
    # common prefix (anchor base) and suffix: the event itself
    pre = 0
    while pre < min(len(ref), len(alt)) and ref[pre] == alt[pre]:
        pre += 1
    suf = 0
    while (
        suf < min(len(ref), len(alt)) - pre and ref[len(ref) - 1 - suf] == alt[len(alt) - 1 - suf]
    ):
        suf += 1
    ev_start, ev_end = start0 + pre, start0 + len(ref) - suf
    ins = alt[pre : len(alt) - suf]
    # a pure insertion or deletion can shift in a repeat; a complex event stays in place
    if ev_start == ev_end or not ins:
        left, right = shift_range(seq, ev_start, ev_end, ins if ev_start == ev_end else "")
    else:
        left, right = ev_start, ev_end
    # flanking bases on each side of every placement and of the REF span
    ws = max(0, min(left, start0) - flank)
    we = min(len(seq), max(right, start0 + len(ref)) + flank)
    ref_hap = seq[ws:we]
    alt_hap = seq[ws:start0] + alt + seq[start0 + len(ref) : we]
    return ws, we, ref_hap, alt_hap


def passes_filters(read, args):
    """True if the read is counted at all (flags and mapping quality)."""
    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
        return False
    if read.is_duplicate and not args.keep_duplicates:
        return False
    return read.mapping_quality >= args.min_mapq


def count_variant(bam, fasta, chrom, pos, ref, alt, args):
    """(ref reads, alt reads, other reads) of one variant."""
    start0 = pos - 1
    counts = [0, 0, 0]
    n = 0
    if len(ref) == len(alt):
        # SNV / MNV: compare the read's bases at the REF positions
        positions = range(start0, start0 + len(ref))
        for read in bam.fetch(chrom, start0, start0 + len(ref)):
            if not passes_filters(read, args):
                continue
            pairs = {r: q for q, r in read.get_aligned_pairs(matches_only=True)}
            if not all(p in pairs for p in positions):
                continue
            qposs = [pairs[p] for p in positions]
            quals = read.query_qualities
            if quals is not None and any(quals[q] < args.min_baseq for q in qposs):
                continue
            bases = "".join(read.query_sequence[q] for q in qposs).upper()
            counts[1 if bases == alt else 0 if bases == ref else 2] += 1
            n += 1
            if n >= args.max_depth:
                break
        return tuple(counts)
    # indel / complex: search the REF and ALT haplotypes of the window in the reads;
    # the reference around the variant (far enough for the largest shift and the flanks)
    margin = MAX_SHIFT + args.flank + 2
    offset = max(0, start0 - margin)
    seq = fasta.fetch(chrom, offset, start0 + len(ref) + len(alt) + margin).upper()
    ws, we, ref_hap, alt_hap = indel_window(seq, start0 - offset, ref, alt, args.flank)
    # window in contig coordinates
    ws, we = ws + offset, we + offset
    # reads near the window: a soft-clipped read with the event aligns next to it
    for read in bam.fetch(chrom, max(0, ws - args.max_read_length), we + args.max_read_length):
        if not passes_filters(read, args) or read.query_sequence is None:
            continue
        read_seq = read.query_sequence.upper()
        if alt_hap in read_seq:
            counts[1] += 1
        elif ref_hap in read_seq:
            counts[0] += 1
        elif read.reference_start <= ws and (read.reference_end or 0) >= we:
            # aligned across the window, but neither haplotype: another allele or errors
            counts[2] += 1
        else:
            # the read doesn't cover the window: no information about this variant
            continue
        n += 1
        if n >= args.max_depth:
            break
    return tuple(counts)


def write_vcf(path, sample, contigs, results):
    """Bgzipped, indexed VCF with DP, AD and OT per variant (sorted by contig and position)."""
    header = pysam.VariantHeader()
    for name, length in contigs.items():
        header.contigs.add(name, length=length)
    header.formats.add("DP", 1, "Integer", "Reads counted at the variant (REF + ALT + other)")
    header.formats.add("AD", "R", "Integer", "Reads supporting REF and ALT")
    header.formats.add("OT", 1, "Integer", "Reads with another allele")
    header.add_line("##source=read_evidence.py")
    header.add_sample(sample)
    order = {name: i for i, name in enumerate(contigs)}
    with pysam.VariantFile(path, "wz", header=header) as out:
        for (chrom, pos, ref, alt), (n_ref, n_alt, n_other) in sorted(
            results.items(), key=lambda kv: (order[kv[0][0]], kv[0][1], kv[0][2], kv[0][3])
        ):
            rec = out.new_record(contig=chrom, start=pos - 1, alleles=(ref, alt))
            rec.samples[sample]["DP"] = n_ref + n_alt + n_other
            rec.samples[sample]["AD"] = (n_ref, n_alt)
            rec.samples[sample]["OT"] = n_other
            out.write(rec)
    pysam.tabix_index(path, preset="vcf", force=True)


def main():
    """Parse the command line, count the reads of every variant, write the VCF."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--alignment", required=True, help="BAM or CRAM (indexed)")
    parser.add_argument("--fasta", required=True, help="reference FASTA (indexed)")
    parser.add_argument("--sites", nargs="+", required=True, help="CHROM POS REF ALT tables")
    parser.add_argument("--sample", required=True, help="sample name of the output VCF")
    parser.add_argument("--output", required=True, help="output VCF (.vcf.gz)")
    parser.add_argument("--min-mapq", type=int, default=20, help="minimum mapping quality")
    parser.add_argument("--min-baseq", type=int, default=13, help="minimum base quality")
    parser.add_argument("--max-depth", type=int, default=100000, help="reads per variant")
    parser.add_argument(
        "--keep-duplicates", action="store_true", help="count duplicate-flagged reads"
    )
    parser.add_argument(
        "--flank", type=int, default=10, help="reference bases on each side of an indel window"
    )
    parser.add_argument(
        "--max-read-length",
        type=int,
        default=300,
        help="reads are taken up to this distance from an indel window",
    )
    args = parser.parse_args()
    if args.min_mapq < 0 or args.min_baseq < 0 or args.max_depth < 1:
        parser.error("--min-mapq/--min-baseq must be >= 0 and --max-depth >= 1")
    if args.flank < 1 or args.max_read_length < 1:
        parser.error("--flank and --max-read-length must be >= 1")

    contigs = read_fai(f"{args.fasta}.fai")
    sites = read_sites(args.sites)
    unknown = sorted({s[0] for s in sites} - set(contigs))
    if unknown:
        sys.exit(
            f"ERROR: {', '.join(args.sites)}: contig(s) {', '.join(unknown)} not in the "
            f"reference {args.fasta}; the sites must be on the reference contigs"
        )
    fasta = pysam.FastaFile(args.fasta)
    bam = open_alignment(args.alignment, args.fasta)
    results = {}
    for site in sorted(sites):
        chrom, pos, ref, alt = site
        # REF must match the reference (the compared VCFs were checked and normalized)
        genome_ref = fasta.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper()
        if genome_ref != ref:
            sys.exit(
                f"ERROR: {chrom}:{pos} {ref}>{alt}: REF differs from the reference "
                f"({genome_ref}); the sites must be normalized against --fasta"
            )
        # contigs missing in the alignment (e.g. a panel BAM without chrY): no reads
        if chrom not in bam.references:
            results[site] = (0, 0, 0)
            continue
        results[site] = count_variant(bam, fasta, chrom, pos, ref, alt, args)
    write_vcf(args.output, args.sample, contigs, results)
    n_alt = sum(1 for r in results.values() if r[1] > 0)
    print(
        f"INFO: {args.alignment}: {len(results)} variants counted ({n_alt} with alt reads); "
        f"filters MQ >= {args.min_mapq}, BQ >= {args.min_baseq}, max depth {args.max_depth}"
        f"{', duplicates kept' if args.keep_duplicates else ''}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
