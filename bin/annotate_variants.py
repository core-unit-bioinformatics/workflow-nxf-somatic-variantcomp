#!/usr/bin/env python3
"""Annotate the compared variants with genes and stratifications.

The interval work is done by bedtools (nf-core modules); this script only
converts the variant table to BED and merges the bedtools results back:

  to-bed   variant table (TSV with columns chrom, pos, ref) -> BED4 with the REF
           span of every variant and its row number as name
           (start = pos - 1, end = pos - 1 + len(ref))
  merge    variant table + the bedtools results -> annotated table:
           --genes    output of `bedtools map -c 4 -o distinct` on the variant BED
                      and a gene BED: column 5 = comma-separated gene names, '.' if
                      no gene overlaps -> new column `gene`
           --stratum  NAME=FILE, output of `bedtools intersect -c` of the variant
                      BED with a stratification BED: column 5 = number of
                      overlapping intervals -> new column NAME = yes (> 0) / no

Every variant (row number) must appear exactly once in every bedtools result,
otherwise the script stops: a missing row would silently lose an annotation.
"""

import argparse
import csv
import re
import sys

REQUIRED = ["chrom", "pos", "ref"]
NAME_RE = re.compile(r"^[A-Za-z0-9_]+$")


def read_table(path):
    """Header and rows of the variant table; checks the required columns."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        missing = [c for c in REQUIRED if c not in header]
        if missing:
            sys.exit(
                f"ERROR: {path}: missing column(s) {', '.join(missing)} (header: {', '.join(header)})"
            )
        rows = list(reader)
    return header, rows


def to_bed(args):
    _header, rows = read_table(args.table)
    with open(args.output, "w") as fh:
        for i, row in enumerate(rows, start=1):
            # positions are 1-based in the table, BED is 0-based half-open
            try:
                pos = int(row["pos"])
            except ValueError:
                sys.exit(f"ERROR: {args.table}: row {i}: pos '{row['pos']}' is not an integer")
            if not row["ref"]:
                sys.exit(f"ERROR: {args.table}: row {i}: empty ref")
            # the BED interval covers the REF allele; the name is the row number (1-based)
            fh.write(f"{row['chrom']}\t{pos - 1}\t{pos - 1 + len(row['ref'])}\t{i}\n")
    print(f"INFO: {args.table}: {len(rows)} variants written as BED", file=sys.stderr)


def read_bedtools_result(path, n_rows, what):
    """{row number: value of column 5} from a bedtools map/intersect -c result; strict."""
    values = {}
    with open(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                sys.exit(
                    f"ERROR: {path} ({what}): line {line_no}: expected 5 columns (BED4 + "
                    f"result), found {len(fields)}"
                )
            row = int(fields[3])
            if row in values:
                sys.exit(f"ERROR: {path} ({what}): row {row} occurs twice")
            values[row] = fields[4]
    # every variant must be annotated exactly once
    missing = [r for r in range(1, n_rows + 1) if r not in values]
    if missing:
        sys.exit(
            f"ERROR: {path} ({what}): {len(missing)} variant(s) missing (rows "
            f"{', '.join(map(str, missing[:5]))}{', ...' if len(missing) > 5 else ''})"
        )
    return values


def merge(args):
    header, rows = read_table(args.table)
    n = len(rows)
    new_columns = []
    # genes: bedtools map writes '.' if no gene overlaps
    if args.genes:
        genes = read_bedtools_result(args.genes, n, "genes")
        for i, row in enumerate(rows, start=1):
            row["gene"] = genes[i]
        new_columns.append("gene")
        annotated = sum(1 for v in genes.values() if v != ".")
        print(f"INFO: genes: {annotated} of {n} variants overlap a gene", file=sys.stderr)
    # stratifications: bedtools intersect -c writes the number of overlapping intervals
    for item in args.stratum:
        if "=" not in item:
            sys.exit(f"ERROR: --stratum '{item}' must be NAME=FILE")
        name, path = item.split("=", 1)
        if not NAME_RE.match(name) or name in header or name in new_columns:
            sys.exit(f"ERROR: stratum name '{name}' is invalid or already a column")
        counts = read_bedtools_result(path, n, f"stratum {name}")
        for i, row in enumerate(rows, start=1):
            try:
                row[name] = "yes" if int(counts[i]) > 0 else "no"
            except ValueError:
                sys.exit(
                    f"ERROR: {path} (stratum {name}): row {i}: count '{counts[i]}' is not an "
                    "integer; was the file made with bedtools intersect -c?"
                )
        new_columns.append(name)
        n_yes = sum(1 for row in rows if row[name] == "yes")
        print(f"INFO: stratum {name}: {n_yes} of {n} variants inside", file=sys.stderr)
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=header + new_columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    b = sub.add_parser("to-bed", help="variant table -> BED4 (row number as name)")
    b.add_argument("--table", required=True)
    b.add_argument("--output", required=True)
    m = sub.add_parser("merge", help="merge bedtools results into the variant table")
    m.add_argument("--table", required=True)
    m.add_argument("--genes", help="bedtools map result (gene names)")
    m.add_argument("--stratum", nargs="*", default=[], help="NAME=bedtools intersect -c result")
    m.add_argument("--output", required=True)
    args = parser.parse_args()
    to_bed(args) if args.command == "to-bed" else merge(args)


if __name__ == "__main__":
    main()
