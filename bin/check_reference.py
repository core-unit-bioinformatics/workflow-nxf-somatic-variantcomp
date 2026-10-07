#!/usr/bin/env python3
"""Check that a VCF or BED file is on the configured reference genome.

Checks:
  - contig names: exact match, or a safe rename (chr1 <-> 1, chrM <-> MT/M)
  - contig lengths in the VCF header against the reference .fai
    * nuclear contig with a different length -> different genome build (error)
    * mitochondrial contig with a different length -> header-only issue, fixed
      by re-writing the header from the .fai (REF check still applies)
  - REF alleles of all records (VCF) against the FASTA; mismatch rate above
    --max-ref-mismatch -> error (wrong build or wrong reference)
  - records/intervals on contigs absent from the reference (dropped, reported)
  - records/intervals beyond the end of a contig (error)

Outputs (prefix given by --prefix):
  <prefix>.rename_chrs.tsv      old<TAB>new contig names (for bcftools annotate
                                --rename-chrs); empty if nothing to rename
  <prefix>.keep_contigs.txt     reference contig names present in the input
                                (after renaming), one per line
  <prefix>.reference_check.tsv  key<TAB>value report
  <prefix>.harmonized.bed       (--type bed only) first three columns with
                                renamed contigs, intervals on unknown contigs
                                dropped; input order kept (sort afterwards)

Exit code 0 if the input is usable (possibly after the fixes above), 1 otherwise.
"""

import argparse
import gzip
import sys
from pathlib import Path

MITO = {"chrM", "MT", "M", "chrMT"}


def read_fai(path):
    lengths = {}
    with open(path) as fh:
        for line in fh:
            name, length = line.split("\t")[:2]
            lengths[name] = int(length)
    return lengths


def map_contig(name, reference):
    """Return the reference contig name for `name`, or None."""
    if name in reference:
        return name
    candidates = []
    if name in MITO:
        candidates += ["chrM", "MT", "M", "chrMT"]
    if name.startswith("chr"):
        candidates.append(name[3:])
    else:
        candidates.append("chr" + name)
    for candidate in candidates:
        if candidate in reference:
            return candidate
    return None


class Report:
    def __init__(self):
        self.values = {}
        self.errors = []
        self.warnings = []

    def set(self, key, value):
        self.values[key] = value

    def write(self, path, status):
        with open(path, "w") as fh:
            fh.write("key\tvalue\n")
            fh.write(f"status\t{status}\n")
            for key, value in self.values.items():
                fh.write(f"{key}\t{value}\n")
            for msg in self.warnings:
                fh.write(f"warning\t{msg}\n")
            for msg in self.errors:
                fh.write(f"error\t{msg}\n")


def check_header_lengths(header_lengths, mapping, reference, report):
    fixed = []
    for name, length in header_lengths.items():
        ref_name = mapping.get(name)
        if ref_name is None or length is None or length == reference[ref_name]:
            continue
        if ref_name in MITO:
            fixed.append(f"{name}:{length}->{reference[ref_name]}")
            report.warnings.append(
                f"contig {name} length {length} differs from reference {ref_name} "
                f"({reference[ref_name]}); header will be rewritten from the reference"
            )
        else:
            report.errors.append(
                f"contig {name} has length {length} but reference {ref_name} has "
                f"{reference[ref_name]}: the input is on a different genome build"
            )
    report.set("header_length_fixes", ",".join(fixed) or "none")


def finish_mapping(names, reference, report, prefix):
    mapping = {name: map_contig(name, reference) for name in names}
    renames = {old: new for old, new in mapping.items() if new is not None and old != new}
    unknown = sorted(name for name, new in mapping.items() if new is None)
    with open(f"{prefix}.rename_chrs.tsv", "w") as fh:
        fh.writelines(f"{old}\t{new}\n" for old, new in sorted(renames.items()))
    renamed = ",".join(f"{o}->{n}" for o, n in sorted(renames.items()))
    report.set("renamed_contigs", renamed or "none")
    return mapping, unknown


def check_vcf(args, reference, report):
    import pysam

    vcf = pysam.VariantFile(args.input)
    fasta = pysam.FastaFile(args.fasta)
    header_lengths = {name: rec.length for name, rec in vcf.header.contigs.items()}
    report.set("header_contigs", len(header_lengths))

    mapping = {}
    record_contigs = {}
    n_records = checked = mismatches = beyond_end = 0
    mito_checked = mito_mismatches = 0
    for rec in vcf:
        n_records += 1
        record_contigs[rec.chrom] = record_contigs.get(rec.chrom, 0) + 1
        if rec.chrom not in mapping:
            mapping[rec.chrom] = map_contig(rec.chrom, reference)
        ref_name = mapping[rec.chrom]
        if ref_name is None:
            continue
        ref = rec.ref
        if rec.pos - 1 + len(ref) > reference[ref_name]:
            beyond_end += 1
            continue
        observed = fasta.fetch(ref_name, rec.pos - 1, rec.pos - 1 + len(ref)).upper()
        if "N" in observed:
            continue
        checked += 1
        is_mismatch = observed != ref.upper()
        mismatches += is_mismatch
        if ref_name in MITO:
            mito_checked += 1
            mito_mismatches += is_mismatch
    report.set("records", n_records)

    mapping, unknown = finish_mapping(
        set(header_lengths) | set(record_contigs), reference, report, args.prefix
    )
    check_header_lengths(header_lengths, mapping, reference, report)

    unknown_with_records = {c: record_contigs[c] for c in unknown if c in record_contigs}
    dropped = sum(unknown_with_records.values())
    report.set("records_on_unknown_contigs", dropped)
    if dropped:
        report.warnings.append(
            f"{dropped} records on contigs absent from the reference will be dropped: "
            + ",".join(f"{c}({n})" for c, n in sorted(unknown_with_records.items()))
        )
    if n_records and dropped == n_records:
        report.errors.append("no record is on a reference contig: wrong reference or naming")

    rate = mismatches / checked if checked else 0.0
    report.set("ref_checked", checked)
    report.set("ref_mismatches", mismatches)
    report.set("ref_mismatch_rate", f"{rate:.6f}")
    report.set("records_beyond_contig_end", beyond_end)
    if beyond_end:
        report.errors.append(f"{beyond_end} records lie beyond the end of their contig")
    if rate > args.max_ref_mismatch:
        report.errors.append(
            f"{mismatches} of {checked} REF alleles ({rate:.2%}) do not match the reference "
            f"(threshold {args.max_ref_mismatch:.2%}): wrong genome build or reference"
        )
    if mito_checked and mito_mismatches / mito_checked > args.max_ref_mismatch:
        report.errors.append(
            f"{mito_mismatches} of {mito_checked} mitochondrial REF alleles do not match "
            "the reference: different mitochondrial sequence"
        )
    return mapping


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def check_bed(args, reference, report):
    intervals = []
    with open_text(args.input) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            intervals.append((fields[0], int(fields[1]), int(fields[2])))
    report.set("intervals", len(intervals))
    mapping, unknown = finish_mapping({c for c, _, _ in intervals}, reference, report, args.prefix)
    dropped = sum(1 for c, _, _ in intervals if mapping[c] is None)
    report.set("intervals_on_unknown_contigs", dropped)
    if dropped:
        report.warnings.append(
            f"{dropped} intervals on contigs absent from the reference will be dropped: "
            + ",".join(unknown)
        )
    if intervals and dropped == len(intervals):
        report.errors.append("no interval is on a reference contig: wrong reference or naming")
    with open(f"{args.prefix}.harmonized.bed", "w") as fh:
        fh.writelines(f"{mapping[c]}\t{s}\t{e}\n" for c, s, e in intervals if mapping[c])
    beyond = sum(1 for c, _, e in intervals if mapping[c] and e > reference[mapping[c]])
    report.set("intervals_beyond_contig_end", beyond)
    if beyond:
        report.errors.append(
            f"{beyond} intervals end beyond their contig: the BED is on a different genome build"
        )
    return mapping


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--input", required=True, help="VCF/BCF or BED(.gz) file")
    parser.add_argument("--type", choices=["vcf", "bed"], required=True)
    parser.add_argument("--fai", required=True, help="reference .fai")
    parser.add_argument("--fasta", help="reference FASTA (required for --type vcf)")
    parser.add_argument("--prefix", required=True, help="output prefix")
    parser.add_argument(
        "--max-ref-mismatch",
        type=float,
        default=0.01,
        help="maximum fraction of REF alleles not matching the reference (default: %(default)s)",
    )
    args = parser.parse_args()
    if args.type == "vcf" and not args.fasta:
        parser.error("--fasta is required for --type vcf")

    reference = read_fai(args.fai)
    report = Report()
    report.set("input", Path(args.input).name)
    report.set("reference_contigs", len(reference))
    mapping = (check_vcf if args.type == "vcf" else check_bed)(args, reference, report)

    keep = sorted({new for new in mapping.values() if new is not None}, key=list(reference).index)
    with open(f"{args.prefix}.keep_contigs.txt", "w") as fh:
        fh.writelines(f"{c}\n" for c in keep)

    status = "error" if report.errors else ("fixed" if report.warnings else "ok")
    if status == "ok" and report.values.get("renamed_contigs", "none") != "none":
        status = "fixed"
    report.write(f"{args.prefix}.reference_check.tsv", status)
    for msg in report.warnings:
        print(f"WARNING: {args.input}: {msg}", file=sys.stderr)
    for msg in report.errors:
        print(f"ERROR: {args.input}: {msg}", file=sys.stderr)
    return 1 if report.errors else 0


if __name__ == "__main__":
    sys.exit(main())
