#!/usr/bin/env python3
"""Callable intervals of one alignment from mosdepth output.

Input are the mosdepth outputs of an alignment; --mode selects which one is used:

  --mode base    --quantized <prefix>.quantized.bed.gz (mosdepth --quantize 0:<min>: with the
                 labels LOW_COVERAGE / CALLABLE): every interval labelled
                 --label (default CALLABLE) is callable, i.e. every base in it has
                 a depth >= the threshold used for --quantize.
  --mode region  --regions <prefix>.regions.bed.gz (mosdepth --by <regions>): every region
                 with a mean depth >= --min-coverage is callable. Mean depth can
                 hide low-coverage parts of a region.

The contig names of the alignment are translated to the reference names with the
rename map of check_reference.py (two columns: alignment name, reference name).
Intervals on contigs that are not in the reference .fai are dropped.

Output: BED3 of the callable intervals (unsorted; sort and merge with bedtools).

Why a script: mosdepth writes the depth classes, but selecting a label, applying
the threshold in region mode and renaming contigs is a few lines of logic that
would otherwise be an awk program inside a module.
"""

import argparse
import gzip
import sys


def read_two_columns(path):
    """Mapping first column -> second column of a tab-separated file (may be empty)."""
    mapping = {}
    with open(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                sys.exit(f"ERROR: {path}: line {line_no}: expected two tab-separated columns")
            mapping[fields[0]] = fields[1]
    return mapping


def main():
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

    # alignment contig name -> reference contig name; names already equal are not listed
    rename = read_two_columns(args.rename)
    # reference contigs: first column of the .fai
    with open(args.fai) as fh:
        reference = {line.split("\t", 1)[0] for line in fh if line.strip()}

    opener = gzip.open if args.input.endswith(".gz") else open
    n_in = n_callable = n_unknown_contig = 0
    labels_seen = set()
    with opener(args.input, "rt") as fin, open(args.output, "w") as fout:
        for line_no, line in enumerate(fin, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                sys.exit(
                    f"ERROR: {args.input}: line {line_no}: expected at least 4 columns "
                    f"(chrom, start, end, label/depth), found {len(fields)}"
                )
            n_in += 1
            chrom, start, end, value = fields[0], fields[1], fields[2], fields[3]
            # base mode: the 4th column is the quantize label
            if args.mode == "base":
                labels_seen.add(value)
                callable_ = value == args.label
            # region mode: the 4th column is the mean depth of the region
            else:
                try:
                    callable_ = float(value) >= args.min_coverage
                except ValueError:
                    sys.exit(
                        f"ERROR: {args.input}: line {line_no}: mean depth '{value}' is not "
                        "a number; is this a mosdepth regions.bed.gz file?"
                    )
            if not callable_:
                continue
            # translate the contig name to the reference name
            ref_chrom = rename.get(chrom, chrom)
            if ref_chrom not in reference:
                n_unknown_contig += 1
                continue
            fout.write(f"{ref_chrom}\t{start}\t{end}\n")
            n_callable += 1

    # in base mode, only the two expected labels may occur, otherwise the quantize
    # settings (bins or MOSDEPTH_Q0/Q1 labels) are not the ones this script expects
    unexpected = labels_seen - {args.label, args.low_label}
    if args.mode == "base" and unexpected:
        found = ", ".join(sorted(unexpected))
        sys.exit(
            f"ERROR: {args.input}: unexpected quantize label(s) {found} (expected "
            f"{args.low_label}/{args.label}); run mosdepth with --quantize 0:<min>: and "
            "the labels MOSDEPTH_Q0/MOSDEPTH_Q1"
        )
    if n_in == 0:
        sys.exit(f"ERROR: {args.input}: no intervals; mosdepth output is empty")
    print(
        f"INFO: {args.input}: {n_in} intervals read, {n_callable} callable written, "
        f"{n_unknown_contig} callable intervals on contigs not in the reference dropped",
        file=sys.stderr,
    )
    if n_callable == 0:
        print(
            f"WARNING: {args.input}: no callable interval (all below the threshold)",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
