#!/usr/bin/env python3
"""Summaries of the comparison regions.

targets   size of every technology's target BED and how much of it lies in the
          common comparison region (the intersection of all target BEDs);
          technologies without BED (WGS) are reported as genome-wide.
callable  size of a sample's callable region relative to the common region.

All BEDs are expected to be sorted and merged (bedtools sort | bedtools merge),
so that bp counts are not inflated by overlapping intervals. The script checks
this and stops if intervals overlap.

Why a script: bedtools computes the intervals; the summary table (bp, fractions,
one row per technology) is bookkeeping that would otherwise be awk in a module.
"""

import argparse
import sys


def read_bed(path):
    """Intervals {chrom: [(start, end), ...]} of a sorted, merged BED; checks the format."""
    intervals = {}
    last = {}
    with open(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                sys.exit(f"ERROR: {path}: line {line_no}: expected at least 3 columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                sys.exit(f"ERROR: {path}: line {line_no}: start/end are not integers")
            if end < start:
                sys.exit(f"ERROR: {path}: line {line_no}: end {end} < start {start}")
            # merged input: an interval must start at or after the end of the previous one
            if fields[0] in last and start < last[fields[0]]:
                sys.exit(f"ERROR: {path}: line {line_no}: intervals overlap or are unsorted; "
                         "sort and merge the BED first")  # fmt: skip
            last[fields[0]] = end
            intervals.setdefault(fields[0], []).append((start, end))
    return intervals


def total_bp(intervals):
    return sum(end - start for ivs in intervals.values() for start, end in ivs)


def overlap_bp(a, b):
    """bp of the intersection of two sorted, merged interval sets (two-pointer sweep)."""
    total = 0
    for chrom, ivs_a in a.items():
        ivs_b = b.get(chrom, [])
        i = j = 0
        while i < len(ivs_a) and j < len(ivs_b):
            # overlap of the current pair
            start = max(ivs_a[i][0], ivs_b[j][0])
            end = min(ivs_a[i][1], ivs_b[j][1])
            if end > start:
                total += end - start
            # advance the interval that ends first
            if ivs_a[i][1] < ivs_b[j][1]:
                i += 1
            else:
                j += 1
    return total


def fraction(part, whole):
    return f"{part / whole:.4f}" if whole else "NA"


def targets(args):
    common = read_bed(args.common)
    common_bp = total_bp(common)
    rows = []
    # one row per targeted technology: its BED and how much of it is in the common region
    for item in args.target:
        name, path = item.split("=", 1)
        bed = read_bed(path)
        bp = total_bp(bed)
        inside = overlap_bp(bed, common)
        rows.append(
            (name, "target", sum(len(v) for v in bed.values()), bp, inside, fraction(inside, bp))
        )
    # one row per genome-wide technology (WGS): the genome and the share in the common region
    genome = read_bed(args.genome)
    genome_bp = total_bp(genome)
    for name in args.genome_wide:
        rows.append((name, "genome", sum(len(v) for v in genome.values()), genome_bp, common_bp,
                     fraction(common_bp, genome_bp)))  # fmt: skip
    rows.append(("common", "common", sum(len(v) for v in common.values()), common_bp, common_bp,
                 fraction(common_bp, common_bp)))  # fmt: skip
    if common_bp == 0:
        sys.exit("ERROR: the common comparison region is empty: the target BEDs (and "
                 "--regions_bed/--exclude_bed) don't overlap. Check the BEDs and the genome build.")  # fmt: skip
    with open(args.output, "w") as fh:
        fh.write("region\ttype\tintervals\tbp\tbp_in_common\tfraction_in_common\n")
        fh.writelines("\t".join(str(x) for x in row) + "\n" for row in rows)
    print(f"INFO: common comparison region: {common_bp} bp", file=sys.stderr)


def callable_(args):
    common = read_bed(args.common)
    callable_bed = read_bed(args.callable_bed)
    common_bp = total_bp(common)
    # the callable region was intersected with the common region, so it must lie inside it
    callable_bp = overlap_bp(callable_bed, common)
    if callable_bp != total_bp(callable_bed):
        sys.exit(f"ERROR: {args.callable_bed}: the callable region extends beyond the common "
                 "region; it must be intersected with the common region first")  # fmt: skip
    with open(args.output, "w") as fh:
        fh.write(
            "sample\talignments\tmode\tmin_coverage\tcommon_bp\tcallable_bp\tfraction_callable\n"
        )
        fh.write(f"{args.sample}\t{args.alignments}\t{args.mode}\t{args.min_coverage}\t"
                 f"{common_bp}\t{callable_bp}\t{fraction(callable_bp, common_bp)}\n")  # fmt: skip
    print(f"INFO: {args.sample}: {callable_bp} of {common_bp} bp of the common region are "
          f"callable in all {args.alignments} alignments", file=sys.stderr)  # fmt: skip
    if callable_bp == 0:
        print(f"WARNING: {args.sample}: no callable base; no variant of this sample will be "
              "compared. Check the coverage (coverage/) and --min_coverage.", file=sys.stderr)  # fmt: skip


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    t = sub.add_parser("targets")
    t.add_argument("--common", required=True, help="common comparison region (BED)")
    t.add_argument("--genome", required=True, help="BED of all reference contigs")
    t.add_argument(
        "--target", nargs="*", default=[], help="technology=BED of targeted technologies"
    )
    t.add_argument("--genome-wide", nargs="*", default=[], help="technologies without target BED")
    t.add_argument("--output", required=True)
    c = sub.add_parser("callable")
    c.add_argument("--sample", required=True)
    c.add_argument("--common", required=True)
    c.add_argument("--callable-bed", required=True)
    c.add_argument("--alignments", type=int, required=True, help="number of alignments")
    c.add_argument("--mode", required=True)
    c.add_argument("--min-coverage", required=True)
    c.add_argument("--output", required=True)
    args = parser.parse_args()
    targets(args) if args.command == "targets" else callable_(args)


if __name__ == "__main__":
    main()
