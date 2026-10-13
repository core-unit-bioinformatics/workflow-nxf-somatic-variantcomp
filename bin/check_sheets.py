#!/usr/bin/env python3
"""Validate the input sheets of workflow-nxf-somatic-variantcomp.

Subcommands:
  samplesheet     the samplesheet (--input)
  stratification  the stratification sheet (--stratification)
  known           the known-variant sheet (--known_variants)

Every subcommand reads a CSV sheet, checks it completely (all problems are
reported at once), resolves relative paths against --base-dir (the directory of
the sheet) and writes a validated CSV with absolute paths that the pipeline reads.

The script only checks the *content* of the sheets. Whether the files exist is
checked by Nextflow when it stages them (file(..., checkIfExists: true)),
because a containerized process cannot see files that are not staged into it.

Exit status: 0 if the sheet is valid, 1 otherwise. Messages go to stderr with
ERROR/WARNING/INFO prefixes.

Why a script: nf-schema validates the format of each column, but not rules
across rows (duplicates, one truth row per sample, one BED per technology, ...).
"""

import argparse
import csv
import re
import sys
from pathlib import Path

# values that mark an empty optional field
EMPTY = {"", "n/a", "na", "NA", "N/A"}
# URL schemes: such paths are used as given (Nextflow handles remote files)
URL = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.-]*://")

SAMPLESHEET_REQUIRED = ["sample", "technology", "library_type", "tumor_id", "normal_id", "vcf"]
SAMPLESHEET_OPTIONAL = ["caller", "tumor_bam", "normal_bam", "target_bed"]
# columns of earlier versions that were removed (index files are found automatically)
SAMPLESHEET_REMOVED = ["vcf_index", "tumor_bai", "normal_bai"]
LIBRARY_TYPES = {"wgs", "capture", "amplicon"}

SAMPLE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
TECH_RE = re.compile(r"^[a-z0-9][a-z0-9-]*$")
NAME_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
STRATUM_RE = re.compile(r"^[A-Za-z0-9_]+$")


class Errors:
    """Collects errors and warnings, so that all problems are reported at once."""

    def __init__(self, sheet):
        self.sheet = sheet
        self.errors = []
        self.warnings = []

    def error(self, msg):
        self.errors.append(msg)

    def warn(self, msg):
        self.warnings.append(msg)

    def report_and_exit_on_error(self):
        # print warnings first, then errors; exit 1 if there was any error
        for msg in self.warnings:
            print(f"WARNING: {self.sheet}: {msg}", file=sys.stderr)
        for msg in self.errors:
            print(f"ERROR: {self.sheet}: {msg}", file=sys.stderr)
        if self.errors:
            print(
                f"ERROR: {self.sheet}: {len(self.errors)} problem(s) found, see above. "
                "Fix the sheet and run the pipeline again.",
                file=sys.stderr,
            )
            sys.exit(1)


def read_csv(path, errors, required, optional, removed=()):
    """Read a CSV sheet; check the header; return a list of row dicts (with line numbers)."""
    with open(path, newline="") as fh:
        text = fh.read()
    # a tab in the header usually means the sheet is a TSV of an earlier version
    first_line = text.splitlines()[0] if text.strip() else ""
    if "\t" in first_line and "," not in first_line:
        errors.error("the sheet is tab-separated; it must be comma-separated (CSV)")
        errors.report_and_exit_on_error()
    reader = csv.DictReader(text.splitlines())
    header = [h.strip() for h in (reader.fieldnames or [])]
    if not header:
        errors.error("the sheet is empty (no header line)")
        errors.report_and_exit_on_error()
    # removed columns: explain what to do instead
    for col in removed:
        if col in header:
            errors.error(
                f"column '{col}' is no longer supported: index files are taken from next to "
                "the data file or created by the pipeline; remove the column"
            )
    # required columns must be present
    for col in required:
        if col not in header:
            errors.error(f"required column '{col}' is missing (header: {', '.join(header)})")
    # unknown columns are reported, but don't stop the run
    for col in header:
        if col not in required and col not in optional and col not in removed:
            errors.warn(f"unknown column '{col}' is ignored")
    errors.report_and_exit_on_error()
    rows = []
    # line 1 is the header, so data rows start at line 2
    for line_no, row in enumerate(reader, start=2):
        # strip whitespace around every value; skip completely empty lines
        clean = {k.strip(): (v or "").strip() for k, v in row.items() if k is not None}
        if not any(clean.values()):
            continue
        if None in row:
            errors.error(f"line {line_no}: more values than columns in the header")
        clean["_line"] = line_no
        rows.append(clean)
    if not rows:
        errors.error("the sheet has no data rows")
    return rows


def resolve(value, base_dir):
    """Absolute path or URL for a sheet value; '' for an empty field."""
    if value in EMPTY:
        return ""
    if URL.match(value):
        return value
    path = Path(value)
    # relative paths are relative to the sheet's directory
    return str(path if path.is_absolute() else (Path(base_dir) / path).resolve())


def check_extension(errors, line, column, value, extensions):
    """Error if a (non-empty) path does not end with one of the allowed extensions."""
    if value and not value.endswith(tuple(extensions)):
        errors.error(f"line {line}: {column} '{value}' must end with {' or '.join(extensions)}")


def write_csv(path, columns, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in columns})


# ------------------------------------------------------------------ samplesheet


def check_samplesheet(args):
    errors = Errors(args.sheet)
    rows = read_csv(
        args.sheet, errors, SAMPLESHEET_REQUIRED, SAMPLESHEET_OPTIONAL, SAMPLESHEET_REMOVED
    )
    out = []
    for row in rows:
        line = row["_line"]
        sample, tech, lib = (
            row.get("sample", ""),
            row.get("technology", ""),
            row.get("library_type", ""),
        )
        # per-row format checks (nf-schema checks them too; repeated here for clear messages)
        if not SAMPLE_RE.match(sample):
            errors.error(f"line {line}: sample '{sample}' may only contain A-Z a-z 0-9 _ . -")
        if not TECH_RE.match(tech):
            errors.error(f"line {line}: technology '{tech}' may only contain a-z 0-9 -")
        if lib not in LIBRARY_TYPES:
            errors.error(
                f"line {line}: library_type '{lib}' must be one of "
                f"{', '.join(sorted(LIBRARY_TYPES))}"
            )
        tumor_id = row.get("tumor_id", "")
        normal_id = row.get("normal_id", "")
        if tumor_id in EMPTY:
            errors.error(f"line {line}: tumor_id is required (name of the tumor column in the VCF)")
        if normal_id in EMPTY:
            normal_id = ""
            errors.warn(
                f"line {line}: no normal_id (tumor-only calls): germline variants may "
                "appear as false positives"
            )
        elif normal_id == tumor_id:
            errors.error(f"line {line}: tumor_id and normal_id are identical ('{tumor_id}')")
        # paths: resolved against the sheet directory
        vcf = resolve(row.get("vcf", ""), args.base_dir)
        tumor_bam = resolve(row.get("tumor_bam", ""), args.base_dir)
        normal_bam = resolve(row.get("normal_bam", ""), args.base_dir)
        target_bed = resolve(row.get("target_bed", ""), args.base_dir)
        if not vcf:
            errors.error(f"line {line}: vcf is required")
        check_extension(errors, line, "vcf", vcf, [".vcf", ".vcf.gz", ".bcf"])
        check_extension(errors, line, "tumor_bam", tumor_bam, [".bam", ".cram"])
        check_extension(errors, line, "normal_bam", normal_bam, [".bam", ".cram"])
        check_extension(errors, line, "target_bed", target_bed, [".bed", ".bed.gz"])
        if normal_bam and not tumor_bam:
            errors.error(f"line {line}: normal_bam given without tumor_bam")
        if normal_bam and not normal_id:
            errors.error(f"line {line}: normal_bam given, but normal_id is empty")
        if lib in ("capture", "amplicon") and not target_bed:
            errors.error(
                f"line {line}: technology '{tech}' has library_type '{lib}' and "
                "therefore needs a target_bed"
            )
        out.append(
            {
                "sample": sample,
                "technology": tech,
                "library_type": lib,
                "tumor_id": tumor_id,
                "normal_id": normal_id,
                "caller": row.get("caller", "") or "unknown",
                "is_truth": "true" if tech == args.truth_technology else "false",
                "vcf": vcf,
                "tumor_bam": tumor_bam,
                "normal_bam": normal_bam,
                "target_bed": target_bed,
                "_line": line,
            }
        )

    # cross-row checks
    seen = {}
    for row in out:
        key = (row["sample"], row["technology"])
        if key in seen:
            errors.error(
                f"lines {seen[key]} and {row['_line']}: sample '{key[0]}' and "
                f"technology '{key[1]}' occur twice"
            )
        seen.setdefault(key, row["_line"])
    technologies = sorted({r["technology"] for r in out})
    if len(technologies) < 2:
        errors.error(f"at least two technologies are needed, found: {', '.join(technologies)}")
    if args.truth_technology not in technologies:
        errors.error(
            f"truth technology '{args.truth_technology}' (--truth_technology) is not in "
            f"the samplesheet (technologies: {', '.join(technologies)})"
        )
    else:
        # every sample needs the truth technology and at least one other technology
        for sample in sorted({r["sample"] for r in out}):
            techs = {r["technology"] for r in out if r["sample"] == sample}
            if args.truth_technology not in techs:
                errors.error(
                    f"sample '{sample}' has no row for the truth technology "
                    f"'{args.truth_technology}'"
                )
            elif len(techs) < 2:
                errors.error(
                    f"sample '{sample}' has no technology to compare with the truth "
                    f"technology '{args.truth_technology}'"
                )
    # one target BED and one library type per technology
    for tech in technologies:
        tech_rows = [r for r in out if r["technology"] == tech]
        beds = sorted({r["target_bed"] or "n/a" for r in tech_rows})
        if len(beds) > 1:
            errors.error(
                f"technology '{tech}' must use one target_bed for all samples, "
                f"found: {', '.join(beds)}"
            )
        libs = sorted({r["library_type"] for r in tech_rows})
        if len(libs) > 1:
            errors.error(
                f"technology '{tech}' must have one library_type, found: {', '.join(libs)}"
            )
    # alignments: required by --callable_regions / --read_evidence (and normals by
    # --require_normal_coverage); every missing one is listed
    if args.require_tumor_bam:
        for row in out:
            if not row["tumor_bam"]:
                errors.error(
                    f"line {row['_line']}: no tumor_bam for {row['sample']} {row['technology']}, "
                    f"but --{args.require_tumor_bam} needs alignments; give the BAM/CRAM or "
                    "switch the option off"
                )
    if args.require_normal_bam:
        for row in out:
            if row["normal_id"] and not row["normal_bam"]:
                errors.error(
                    f"line {row['_line']}: no normal_bam for {row['sample']} {row['technology']} "
                    f"(normal_id {row['normal_id']}), but --{args.require_normal_bam} needs it; "
                    "give the BAM/CRAM or switch the option off"
                )
    with_bam = sum(1 for r in out if r["tumor_bam"])
    if 0 < with_bam < len(out):
        errors.warn(
            f"alignments are given for {with_bam} of {len(out)} rows; samples without "
            "alignments for every technology are compared in the common target region "
            "only (no callable region, no read evidence)"
        )

    errors.report_and_exit_on_error()
    columns = [
        "sample",
        "technology",
        "library_type",
        "tumor_id",
        "normal_id",
        "caller",
        "is_truth",
        "vcf",
        "tumor_bam",
        "normal_bam",
        "target_bed",
    ]
    write_csv(args.output, columns, out)
    print(
        f"INFO: {args.sheet}: {len(out)} rows, {len({r['sample'] for r in out})} samples, "
        f"technologies: {', '.join(technologies)} (truth: {args.truth_technology})",
        file=sys.stderr,
    )


# ------------------------------------------------------------------ stratification


def check_stratification(args):
    errors = Errors(args.sheet)
    rows = read_csv(args.sheet, errors, ["name", "bed"], [])
    out, names = [], {}
    for row in rows:
        line, name = row["_line"], row.get("name", "")
        # names become column names of the variant table
        if not STRATUM_RE.match(name):
            errors.error(f"line {line}: name '{name}' may only contain A-Z a-z 0-9 _")
        if name in names:
            errors.error(f"lines {names[name]} and {line}: name '{name}' is used twice")
        names.setdefault(name, line)
        bed = resolve(row.get("bed", ""), args.base_dir)
        if not bed:
            errors.error(f"line {line}: bed is required")
        check_extension(errors, line, "bed", bed, [".bed", ".bed.gz"])
        out.append({"name": name, "bed": bed})
    errors.report_and_exit_on_error()
    write_csv(args.output, ["name", "bed"], out)
    print(f"INFO: {args.sheet}: {len(out)} stratifications", file=sys.stderr)


# ------------------------------------------------------------------ known variants


def check_known(args):
    errors = Errors(args.sheet)
    rows = read_csv(args.sheet, errors, ["name", "sample", "mode", "variants"], ["regions"])
    # samples of the (validated) samplesheet
    with open(args.samplesheet, newline="") as fh:
        samples = {r["sample"] for r in csv.DictReader(fh)}
    out, names = [], {}
    for row in rows:
        line, name, sample, mode = (
            row["_line"],
            row.get("name", ""),
            row.get("sample", ""),
            row.get("mode", ""),
        )
        if not NAME_RE.match(name):
            errors.error(f"line {line}: name '{name}' may only contain A-Z a-z 0-9 _ . -")
        if name in names:
            errors.error(f"lines {names[name]} and {line}: name '{name}' is used twice")
        names.setdefault(name, line)
        if sample != "all" and sample not in samples:
            errors.error(
                f"line {line}: sample '{sample}' is not in the samplesheet "
                f"(samples: {', '.join(sorted(samples))}; or use 'all')"
            )
        if mode not in ("required", "benchmark"):
            errors.error(f"line {line}: mode '{mode}' must be 'required' or 'benchmark'")
        # variants: one or more files separated by ';' (or ',' for compatibility)
        files = [
            resolve(v.strip(), args.base_dir)
            for v in re.split(r"[;,]", row.get("variants", ""))
            if v.strip()
        ]
        if not files:
            errors.error(f"line {line}: variants is required (VCF, CSV or TSV file(s))")
        for f in files:
            check_extension(
                errors, line, "variants", f, [".vcf", ".vcf.gz", ".bcf", ".csv", ".tsv"]
            )
        regions = resolve(row.get("regions", ""), args.base_dir)
        check_extension(errors, line, "regions", regions, [".bed", ".bed.gz"])
        if mode == "benchmark" and not regions:
            errors.warn(
                f"line {line}: benchmark list '{name}' has no regions: every variant "
                "outside the truth set counts as a false positive"
            )
        out.append(
            {
                "name": name,
                "sample": sample,
                "mode": mode,
                "variants": ";".join(files),
                "regions": regions,
            }
        )
    errors.report_and_exit_on_error()
    write_csv(args.output, ["name", "sample", "mode", "variants", "regions"], out)
    print(f"INFO: {args.sheet}: {len(out)} known-variant lists", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("samplesheet", "stratification", "known"):
        p = sub.add_parser(name)
        p.add_argument("--sheet", required=True, help="CSV sheet to check")
        p.add_argument(
            "--base-dir", required=True, help="directory relative paths are resolved against"
        )
        p.add_argument("--output", required=True, help="validated CSV")
        if name == "samplesheet":
            p.add_argument("--truth-technology", required=True)
            p.add_argument(
                "--require-tumor-bam",
                metavar="REASON",
                help="every row needs a tumor_bam (REASON: name of the pipeline param that needs it, e.g. callable_regions)",
            )
            p.add_argument(
                "--require-normal-bam",
                metavar="REASON",
                help="every row with a normal_id needs a normal_bam (REASON: the option)",
            )
        if name == "known":
            p.add_argument("--samplesheet", required=True, help="validated samplesheet")
    args = parser.parse_args()
    {
        "samplesheet": check_samplesheet,
        "stratification": check_stratification,
        "known": check_known,
    }[args.command](args)


if __name__ == "__main__":
    main()
