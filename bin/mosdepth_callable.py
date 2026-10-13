#!/usr/bin/env python3
"""Callable intervals of one alignment from mosdepth output.

INPUT (outputs of mosdepth for one alignment; --mode selects which one is used)
  --quantized  <prefix>.quantized.bed.gz of `mosdepth --quantize 0:<min>:` with the labels
               LOW_COVERAGE / CALLABLE (env MOSDEPTH_Q0/MOSDEPTH_Q1): columns contig,
               start (0-based), end, label. Used with --mode base: every interval labelled
               --label (default CALLABLE) is callable, i.e. every base in it has a depth
               >= the threshold given to --quantize.
  --regions    <prefix>.regions.bed.gz of `mosdepth --by <regions>`: columns contig, start,
               end, mean depth. Used with --mode region: every region with a mean depth
               >= --min-coverage is callable (a mean can hide low-coverage parts).
  --rename     rename map of the alignment's reference check (check_reference.py): two
               columns, alignment contig name -> reference contig name; may be empty
  --fai        reference .fai: intervals on other contigs are dropped

OUTPUT
  --output     BED3 (contig in reference names, start, end) of the callable intervals,
               unsorted (the pipeline sorts and merges it with bedtools)

Checks: base mode accepts only the two expected labels (other labels mean other
mosdepth settings); region mode needs numeric depths; empty mosdepth output is an error.

Why a script: mosdepth writes the depth classes, but selecting a label, applying the
threshold in region mode and renaming contigs is a few lines of logic that would
otherwise be an awk program inside a module.
"""

import argparse
import gzip
import sys


def open_mosdepth_bed(path):
    """Text handle of a mosdepth BED (bgzipped as written by mosdepth, or plain)."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_rename_map(path):
    """{alignment contig name: reference contig name}; names that are equal are not listed."""
    rename_map = {}
    with open(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                sys.exit(f"ERROR: {path}: line {line_no}: expected two tab-separated columns")
            rename_map[fields[0]] = fields[1]
    return rename_map


def read_reference_contigs(fai_path):
    """Set of the reference contig names (first column of the .fai)."""
    with open(fai_path) as fh:
        return {line.split("\t", 1)[0] for line in fh if line.strip()}


def is_callable(value, mode, label, min_coverage, where):
    """Is a mosdepth interval callable? value = label (base mode) or mean depth (region)."""
    # base mode: the 4th column is the quantize label of the interval
    if mode == "base":
        return value == label
    # region mode: the 4th column is the mean depth of the region
    try:
        return float(value) >= min_coverage
    except ValueError:
        sys.exit(
            f"ERROR: {where}: mean depth '{value}' is not a number; is this a mosdepth regions.bed.gz file?"
        )


def write_callable(args, rename_map, reference_contigs):
    """Write the callable intervals; return counts and the labels seen (base mode)."""
    counts = {"intervals": 0, "callable": 0, "unknown_contig": 0}
    labels_seen = set()
    with open_mosdepth_bed(args.input) as mosdepth_bed, open(args.output, "w") as callable_bed:
        for line_no, line in enumerate(mosdepth_bed, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                sys.exit(
                    f"ERROR: {args.input}: line {line_no}: expected at least 4 columns "
                    f"(contig, start, end, label/depth), found {len(fields)}"
                )
            counts["intervals"] += 1
            contig, start, end, value = fields[:4]
            if args.mode == "base":
                labels_seen.add(value)
            if not is_callable(
                value, args.mode, args.label, args.min_coverage, f"{args.input}: line {line_no}"
            ):
                continue
            # contig name of the alignment -> reference name; other contigs are dropped
            reference_contig = rename_map.get(contig, contig)
            if reference_contig not in reference_contigs:
                counts["unknown_contig"] += 1
                continue
            callable_bed.write(f"{reference_contig}\t{start}\t{end}\n")
            counts["callable"] += 1
    return counts, labels_seen


def main():
    """Parse the command line and run the (sub)command."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--quantized", required=True, help="mosdepth <prefix>.quantized.bed.gz")
    parser.add_argument("--regions", required=True, help="mosdepth <prefix>.regions.bed.gz")
    parser.add_argument("--mode", choices=["base", "region"], required=True)
    parser.add_argument("--label", default="CALLABLE", help="quantize label of callable bins")
    parser.add_argument(
        "--low-label", default="LOW_COVERAGE", help="quantize label of bins below the threshold"
    )
    parser.add_argument("--min-coverage", type=float, help="minimum mean depth (region mode)")
    parser.add_argument("--rename", required=True, help="contig rename map (may be empty)")
    parser.add_argument("--fai", required=True, help="reference .fai")
    parser.add_argument("--output", required=True, help="output BED")
    args = parser.parse_args()

    # region mode needs the threshold
    if args.mode == "region" and args.min_coverage is None:
        parser.error("--min-coverage is required with --mode region")
    # the mosdepth output that belongs to the mode
    args.input = args.quantized if args.mode == "base" else args.regions
    print(f"INFO: callable mode '{args.mode}': using {args.input}", file=sys.stderr)

    rename_map = read_rename_map(args.rename)
    reference_contigs = read_reference_contigs(args.fai)
    counts, labels_seen = write_callable(args, rename_map, reference_contigs)

    # base mode: only the two expected labels may occur, otherwise the quantize settings
    # (bins or MOSDEPTH_Q0/Q1 labels) are not the ones this script expects
    unexpected = labels_seen - {args.label, args.low_label}
    if args.mode == "base" and unexpected:
        sys.exit(
            f"ERROR: {args.input}: unexpected quantize label(s) {', '.join(sorted(unexpected))} "
            f"(expected {args.low_label}/{args.label}); run mosdepth with --quantize 0:<min>: and "
            "the labels MOSDEPTH_Q0/MOSDEPTH_Q1"
        )
    if counts["intervals"] == 0:
        sys.exit(f"ERROR: {args.input}: no intervals; mosdepth output is empty")
    print(
        f"INFO: {args.input}: {counts['intervals']} intervals read, {counts['callable']} callable "
        f"written, {counts['unknown_contig']} callable intervals on contigs not in the reference "
        "dropped",
        file=sys.stderr,
    )
    if counts["callable"] == 0:
        print(
            f"WARNING: {args.input}: no callable interval (all below the threshold)",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
