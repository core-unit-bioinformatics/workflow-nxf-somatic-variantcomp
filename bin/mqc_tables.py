#!/usr/bin/env python3
"""MultiQC custom-content tables from the summaries of the pipeline steps.

Every step that writes a summary (input checks, harmonization, bcftools norm, regions)
gets one table section in the MultiQC report, so that all intermediate counts can be
followed in one place.

Subcommands (all write <--prefix>_mqc.tsv, a MultiQC custom-content table)
  key-value  summaries with two columns `key value` (check_reference.py,
             known_variants.py prepare, harmonize_vcf.py)
             INPUT   one or more files; the row name is the file name without the given
                     --strip suffix (e.g. 'S1_wes.vcf' of S1_wes.vcf.reference_check.tsv)
             OUTPUT  one row per file, one column per key (first occurrence; the
                     repeated keys `warning` and `error` are counted instead)
  norm-logs  logs of `bcftools norm` (stderr)
             INPUT   one or more logs with the line
                     'Lines total/split/joined/realigned/mismatch_removed/dup_removed/skipped: a/b/...'
             OUTPUT  one row per log: the seven counts
  table      tables with a header (region_summary.py)
             INPUT   one or more TSVs with the same header; the first --id-columns
                     columns form the row name
             OUTPUT  all rows, all columns

Options --section (MultiQC id), --title and --description name the section.

Why a script: the summaries have different formats; turning them into MultiQC tables
is bookkeeping without a tool.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

NORM_COLUMNS = [
    "total",
    "split",
    "joined",
    "realigned",
    "mismatch_removed",
    "dup_removed",
    "skipped",
]
NORM_LINE = re.compile(
    r"total/split/joined/realigned/mismatch_removed/dup_removed/skipped:\s*([\d/]+)"
)


def mqc_header(args):
    """Comment lines that define the MultiQC section."""
    return (
        f"# id: 'variantcomp_{args.section}'\n"
        f"# section_name: '{args.title}'\n"
        f"# description: '{args.description}'\n"
        "# plot_type: 'table'\n"
    )


def row_name(path, strip):
    """Row name of a file: its name without the --strip suffix."""
    name = Path(path).name
    return name[: -len(strip)] if strip and name.endswith(strip) else name


def write_table(path, header_text, columns, rows):
    """MultiQC table: header comments, a header line ('id' + columns), one line per row."""
    with open(path, "w") as fh:
        fh.write(header_text)
        fh.write("\t".join(["id", *columns]) + "\n")
        for name, values in rows:
            fh.write("\t".join([name, *(str(values.get(c, "NA")) for c in columns)]) + "\n")


def key_value(args):
    """One row per key/value summary; warnings and errors counted."""
    rows, columns = [], []
    for path in sorted(args.inputs):
        values = {"warnings": 0, "errors": 0}
        with open(path) as fh:
            lines = fh.read().splitlines()
        if not lines or lines[0].split("\t") != ["key", "value"]:
            sys.exit(f"ERROR: {path}: expected a two-column summary with the header 'key value'")
        for line_no, line in enumerate(lines[1:], start=2):
            fields = line.split("\t", 1)
            if len(fields) != 2:
                sys.exit(f"ERROR: {path}: line {line_no}: expected two tab-separated columns")
            key, value = fields
            # repeated message rows are counted, the other keys keep their first value
            if key in ("warning", "error"):
                values[f"{key}s"] += 1
            else:
                values.setdefault(key, value)
        for key in values:
            if key not in columns:
                columns.append(key)
        rows.append((row_name(path, args.strip), values))
    write_table(f"{args.prefix}_mqc.tsv", mqc_header(args), columns, rows)
    print(f"INFO: {args.section}: {len(rows)} summaries, {len(columns)} columns", file=sys.stderr)


def norm_logs(args):
    """One row per bcftools norm log with its seven counts."""
    rows = []
    for path in sorted(args.inputs):
        match = None
        with open(path) as fh:
            for line in fh:
                match = NORM_LINE.search(line) or match
        if match is None:
            sys.exit(f"ERROR: {path}: no summary line of bcftools norm found")
        counts = match.group(1).split("/")
        if len(counts) != len(NORM_COLUMNS):
            sys.exit(f"ERROR: {path}: unexpected bcftools norm summary '{match.group(1)}'")
        rows.append((row_name(path, args.strip), dict(zip(NORM_COLUMNS, counts, strict=True))))
    write_table(f"{args.prefix}_mqc.tsv", mqc_header(args), NORM_COLUMNS, rows)
    print(f"INFO: {args.section}: {len(rows)} logs", file=sys.stderr)


def table(args):
    """All rows of TSVs with one common header; the first --id-columns form the row name."""
    rows, header = [], None
    for path in sorted(args.inputs):
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if header is None:
                header = reader.fieldnames
            elif reader.fieldnames != header:
                sys.exit(
                    f"ERROR: {path}: header differs from the other tables ({', '.join(header)})"
                )
            for row in reader:
                name = " ".join(row[c] for c in header[: args.id_columns])
                rows.append((name, row))
    columns = header[args.id_columns :] if header else []
    write_table(f"{args.prefix}_mqc.tsv", mqc_header(args), columns, rows)
    print(f"INFO: {args.section}: {len(rows)} rows", file=sys.stderr)


def main():
    """Parse the command line and run the (sub)command."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("key-value", "norm-logs", "table"):
        p = sub.add_parser(name)
        p.add_argument("inputs", nargs="+", help="summary files")
        p.add_argument("--section", required=True, help="MultiQC section id (a-z, 0-9, _)")
        p.add_argument("--title", required=True, help="section title")
        p.add_argument("--description", default="", help="section description")
        p.add_argument("--prefix", required=True, help="output prefix")
        if name != "table":
            p.add_argument("--strip", default="", help="file-name suffix removed for the row name")
        else:
            p.add_argument("--id-columns", type=int, default=1, help="columns of the row name")
    args = parser.parse_args()
    if not re.fullmatch(r"[a-z0-9_]+", args.section):
        parser.error("--section may only contain a-z, 0-9 and _")
    {"key-value": key_value, "norm-logs": norm_logs, "table": table}[args.command](args)


if __name__ == "__main__":
    main()
