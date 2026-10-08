#!/usr/bin/env python3
"""Summarise the 'required' known-variant checks of all samples.

Input: the per-sample results of `known_variants.py check` (<x>.known_variants.tsv
and <x>.known_summary.tsv), each with the list name given as NAME=FILE.

Outputs (--prefix):
  <prefix>.required_variants.tsv  all variants of all lists, samples and technologies
  <prefix>.required_summary.tsv   per list, sample and technology: required, passed,
                                  failed, all_passed
  <prefix>_required_mqc.tsv       the summary as a MultiQC table

Every variant that is not PASS is printed as a WARNING (list, label, position,
sample, technology, status). With --fail-on-missing the script exits with status
1 after writing the outputs, so that the pipeline stops (--fail_on_missing_required).
"""

import argparse
import csv
import sys

VARIANT_COLUMNS = ["sample", "technology", "label", "gene", "chrom", "pos", "ref", "alt", "status"]
SUMMARY_COLUMNS = ["sample", "technology", "required", "passed", "failed", "all_passed"]


def read_tsv(path, required, what):
    """Rows of a TSV; stops if a required column is missing."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in required if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"ERROR: {path} ({what}): missing column(s) {', '.join(missing)}")
        return list(reader)


def parse_pairs(items, option):
    """NAME=FILE arguments -> list of (name, file)."""
    pairs = []
    for item in items:
        if "=" not in item:
            sys.exit(f"ERROR: {option} '{item}' must be LIST_NAME=FILE")
        pairs.append(tuple(item.split("=", 1)))
    return pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--variants", nargs="+", required=True, help="LIST=<x>.known_variants.tsv")
    parser.add_argument("--summary", nargs="+", required=True, help="LIST=<x>.known_summary.tsv")
    parser.add_argument(
        "--fail-on-missing", action="store_true", help="exit 1 if a required variant is missing"
    )
    parser.add_argument("--prefix", default="known")
    args = parser.parse_args()

    # all variants of all lists; the list name becomes the first column
    variants = []
    for name, path in parse_pairs(args.variants, "--variants"):
        for row in read_tsv(path, VARIANT_COLUMNS, "known variants"):
            variants.append({"list": name, **row})
    summary = []
    for name, path in parse_pairs(args.summary, "--summary"):
        for row in read_tsv(path, SUMMARY_COLUMNS, "known summary"):
            summary.append({"list": name, **row})

    # sort for stable outputs
    variants.sort(
        key=lambda r: (r["list"], r["sample"], r["technology"], r["chrom"], int(r["pos"]))
    )
    summary.sort(key=lambda r: (r["list"], r["sample"], r["technology"]))

    var_columns = (
        ["list"] + [c for c in variants[0] if c != "list"]
        if variants
        else ["list", *VARIANT_COLUMNS]
    )
    with open(f"{args.prefix}.required_variants.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=var_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(variants)
    with open(f"{args.prefix}.required_summary.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["list", *SUMMARY_COLUMNS], delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(summary)

    # MultiQC custom content: one table row per list, sample and technology
    with open(f"{args.prefix}_required_mqc.tsv", "w") as fh:
        fh.write("# id: 'variantcomp_required'\n# section_name: 'Required variants'\n"
                 "# description: 'Known variants that must be found in every technology "
                 "(--known_variants, mode required).'\n# plot_type: 'table'\n")  # fmt: skip
        fh.write("id\t" + "\t".join(["list", *SUMMARY_COLUMNS]) + "\n")
        for row in summary:
            row_id = f"{row['list']} {row['sample']} {row['technology']}"
            fh.write(row_id + "\t" + "\t".join(row[c] for c in ["list", *SUMMARY_COLUMNS]) + "\n")

    # one warning per variant that was not found as PASS
    failed = [r for r in variants if r["status"] != "PASS"]
    for r in failed:
        print(f"WARNING: required variant {r['label']} {r['chrom']}:{r['pos']} {r['ref']}>{r['alt']} "
              f"({r['gene']}) of list '{r['list']}' not found in {r['sample']} {r['technology']}: "
              f"{r['status']}", file=sys.stderr)  # fmt: skip
    print(f"INFO: {len(variants)} required variant checks, {len(failed)} not PASS", file=sys.stderr)
    if failed and args.fail_on_missing:
        print(f"ERROR: {len(failed)} required variant check(s) failed, see "
              f"{args.prefix}.required_variants.tsv (--fail_on_missing_required)", file=sys.stderr)  # fmt: skip
        sys.exit(1)


if __name__ == "__main__":
    main()
