#!/usr/bin/env python3
"""Summaries of the comparison regions.

targets   size of every technology's target BED and how much of it lies in the
          common comparison region (the intersection of all target BEDs);
          technologies without BED (WGS) are reported as genome-wide.
          INPUT   --common BED, --genome BED (all compared contigs), --target
                  TECHNOLOGY=BED ..., --genome-wide TECHNOLOGY ...
          OUTPUT  TSV with columns region (technology or 'common'), type (target,
                  genome, common), intervals (count), bp, bp_in_common,
                  fraction_in_common (bp_in_common / bp, 4 digits)
callable  size of a sample's callable region relative to the common region.
          INPUT   --common BED, --callable-bed BED (already intersected with it)
          OUTPUT  TSV with columns sample, alignments (number), mode (base/region),
                  min_coverage (tumor/normal), common_bp, callable_bp,
                  fraction_callable
steps     size of the comparison region after every step of PREPARE_REGIONS.
          INPUT   --step NAME=BED ... in the order of the steps (genome, intersection of
                  the target BEDs, after --regions_bed, after --exclude_bed)
          OUTPUT  TSV with columns step (1, 2, ...), name, intervals, bp,
                  fraction_of_genome (bp / bp of step 1), fraction_of_previous
                  (bp / bp of the step before), removed_bp (bp of the step before - bp)
--min-bp: the common resp. callable region must have at least this many bp (error).

Overlapping or adjacent intervals of a BED (common in target BEDs of capture kits) are
merged before counting, so that bp are not counted twice; the number of merged
intervals is logged.

Why a script: bedtools computes the intervals; the summary table (bp, fractions,
one row per technology) is bookkeeping that would otherwise be awk in a module.
"""

import argparse
import sys
from pathlib import Path

# BED reader shared with the reference check (format checks with line numbers)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_reference import read_bed  # noqa: E402


def read_merged_bed(path):
    """Intervals {chrom: [(start, end), ...]} of a BED, sorted and merged per contig.

    The BED format is checked by check_reference.read_bed (columns, integer coordinates,
    0 <= start <= end). Overlapping or adjacent intervals are merged (target BEDs of
    capture kits often overlap), so that no base is counted twice.
    """
    by_chrom = {}
    for iv in read_bed(path):
        by_chrom.setdefault(iv.chrom, []).append((iv.start, iv.end))
    intervals, n_merged = {}, 0
    for chrom, ivs in by_chrom.items():
        merged = []
        for start, end in sorted(ivs):
            # overlaps or touches the previous interval: extend it
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                n_merged += 1
            else:
                merged.append((start, end))
        intervals[chrom] = merged
    if n_merged:
        print(f"INFO: {path}: {n_merged} overlapping or adjacent intervals merged", file=sys.stderr)
    return intervals


def total_bp(intervals):
    """Total size (bp) of sorted, merged intervals."""
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
    """part / whole with 4 digits, NA if whole is 0."""
    return f"{part / whole:.4f}" if whole else "NA"


def targets(args):
    """Summary of the target BEDs and the common region (subcommand targets)."""
    common = read_merged_bed(args.common)
    common_bp = total_bp(common)
    rows = []
    # one row per targeted technology: its BED and how much of it is in the common region
    for item in args.target:
        name, path = item.split("=", 1)
        bed = read_merged_bed(path)
        bp = total_bp(bed)
        inside = overlap_bp(bed, common)
        rows.append(
            (name, "target", sum(len(v) for v in bed.values()), bp, inside, fraction(inside, bp))
        )
    # one row per genome-wide technology (WGS): the genome and the share in the common region
    genome = read_merged_bed(args.genome)
    genome_bp = total_bp(genome)
    for name in args.genome_wide:
        rows.append(
            (
                name,
                "genome",
                sum(len(v) for v in genome.values()),
                genome_bp,
                common_bp,
                fraction(common_bp, genome_bp),
            )
        )
    rows.append(
        (
            "common",
            "common",
            sum(len(v) for v in common.values()),
            common_bp,
            common_bp,
            fraction(common_bp, common_bp),
        )
    )
    # the comparison region must be large enough for a meaningful benchmark
    if 0 < common_bp < args.min_bp:
        sys.exit(
            f"ERROR: the common comparison region has {common_bp} bp, less than --min_region_bp "
            f"{args.min_bp}: the target BEDs (and --regions_bed/--exclude_bed) overlap too little "
            "for a valid benchmark. Check the BEDs or lower --min_region_bp."
        )
    if common_bp == 0:
        sys.exit(
            "ERROR: the common comparison region is empty: the target BEDs (and "
            "--regions_bed/--exclude_bed) don't overlap. Check the BEDs and the genome build."
        )
    with open(args.output, "w") as fh:
        fh.write("region\ttype\tintervals\tbp\tbp_in_common\tfraction_in_common\n")
        fh.writelines("\t".join(str(x) for x in row) + "\n" for row in rows)
    print(f"INFO: common comparison region: {common_bp} bp", file=sys.stderr)


def callable_(args):
    """Summary of the callable region of one sample (subcommand callable)."""
    common = read_merged_bed(args.common)
    callable_bed = read_merged_bed(args.callable_bed)
    common_bp = total_bp(common)
    # the callable region was intersected with the common region, so it must lie inside it
    callable_bp = overlap_bp(callable_bed, common)
    if callable_bp != total_bp(callable_bed):
        sys.exit(
            f"ERROR: {args.callable_bed}: the callable region extends beyond the common "
            "region; it must be intersected with the common region first"
        )
    with open(args.output, "w") as fh:
        fh.write(
            "sample\talignments\tmode\tmin_coverage\tcommon_bp\tcallable_bp\tfraction_callable\n"
        )
        fh.write(
            f"{args.sample}\t{args.alignments}\t{args.mode}\t{args.min_coverage}\t"
            f"{common_bp}\t{callable_bp}\t{fraction(callable_bp, common_bp)}\n"
        )
    print(
        f"INFO: {args.sample}: {callable_bp} of {common_bp} bp of the common region are "
        f"callable in all {args.alignments} alignments",
        file=sys.stderr,
    )
    # the callable region of the sample must be large enough for a meaningful benchmark
    if callable_bp < args.min_bp:
        sys.exit(
            f"ERROR: {args.sample}: only {callable_bp} bp of the common region are callable in all "
            f"{args.alignments} alignments, less than --min_region_bp {args.min_bp}. Check the "
            "coverage (coverage/), --min_coverage_tumor/--min_coverage_normal, or lower "
            "--min_region_bp."
        )
    if callable_bp == 0:
        print(
            f"WARNING: {args.sample}: no callable base; no variant of this sample will be "
            "compared. Check the coverage (coverage/) and --min_coverage.",
            file=sys.stderr,
        )


def steps(args):
    """Size of the comparison region after every step (subcommand steps)."""
    rows, genome_bp, previous_bp = [], None, None
    for number, item in enumerate(args.step, start=1):
        name, path = item.split("=", 1)
        bed = read_merged_bed(path)
        bp = total_bp(bed)
        # step 1 is the genome (all compared contigs): the reference for the fractions
        genome_bp = bp if genome_bp is None else genome_bp
        # a step can only remove bases (intersection, subtraction); more bp is an error
        if previous_bp is not None and bp > previous_bp:
            sys.exit(
                f"ERROR: step {number} ({name}): {bp} bp, more than the step before "
                f"({previous_bp} bp); every region step must be a subset of the step before"
            )
        removed = "NA" if previous_bp is None else str(previous_bp - bp)
        rows.append(
            (
                number,
                name,
                sum(len(v) for v in bed.values()),
                bp,
                fraction(bp, genome_bp),
                "NA" if previous_bp is None else fraction(bp, previous_bp),
                removed,
            )
        )
        previous_bp = bp
    with open(args.output, "w") as fh:
        fh.write(
            "step\tname\tintervals\tbp\tfraction_of_genome\tfraction_of_previous\tremoved_bp\n"
        )
        fh.writelines("\t".join(str(x) for x in row) + "\n" for row in rows)
    print(f"INFO: region steps: {' -> '.join(str(r[3]) for r in rows)} bp", file=sys.stderr)


def main():
    """Parse the command line and run the (sub)command."""
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
    t.add_argument("--min-bp", type=int, default=0, help="minimum size of the common region")
    st = sub.add_parser("steps")
    st.add_argument("--step", nargs="+", required=True, help="NAME=BED in the order of the steps")
    st.add_argument("--output", required=True)
    c = sub.add_parser("callable")
    c.add_argument("--sample", required=True)
    c.add_argument("--common", required=True)
    c.add_argument("--callable-bed", required=True)
    c.add_argument("--alignments", type=int, required=True, help="number of alignments")
    c.add_argument("--mode", required=True)
    c.add_argument("--min-coverage", required=True)
    c.add_argument("--output", required=True)
    c.add_argument("--min-bp", type=int, default=0, help="minimum size of the callable region")
    args = parser.parse_args()
    {"targets": targets, "callable": callable_, "steps": steps}[args.command](args)


if __name__ == "__main__":
    main()
