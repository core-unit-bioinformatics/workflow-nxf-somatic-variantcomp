#!/usr/bin/env python3
"""Decide whether an input is on the configured reference genome.

The facts are extracted by tools in the modules before this script runs:

  VCF  header         `bcftools view -h`            -> contig lines and sample names
       record contigs `bcftools query -f '%CHROM\\n'` -> contig of every record
  BAM  contigs        `samtools idxstats`           -> contig names and lengths
  BED  the BED itself (plain text; only the contig names and interval ends are used)
  VCF  REF alleles    log of `bcftools norm --check-ref x` (total and removed
                      REF-mismatch records)

This script makes the decisions (it reads no FASTA and no binary file):

subcommand `contigs` (--type vcf, bed or bam)
  - contig names: exact match with the reference .fai, or a safe rename
    (chr1 <-> 1, chrM <-> MT/M/chrMT)
  - contig lengths (VCF header, BAM): a nuclear contig with another length means
    another genome build -> error; a mitochondrial contig with another length is a
    known header-only problem (e.g. chrM 16570 vs 16571) -> fixed by rewriting the
    VCF header from the .fai (the REF check still applies to its records)
  - records/intervals on contigs absent from the reference -> dropped and reported
  - BED intervals that end beyond their contig -> error (another genome build)
  - VCF: the sample columns given with --samples must exist -> error otherwise
  Outputs (--prefix):
    <prefix>.rename_chrs.tsv      old<TAB>new names for `bcftools annotate --rename-chrs`
    <prefix>.keep_contigs.txt     reference contigs present in the input
    <prefix>.reference_check.tsv  key<TAB>value report (status ok / fixed / error)
    <prefix>.harmonized.bed       (bed) all columns, renamed, unknown contigs dropped
    <prefix>.regions.bed          (bam, with --regions) a reference-named BED translated
                                  to the alignment's contig names (for mosdepth --by)

subcommand `ref-alleles`
  - reads the summary line of `bcftools norm --check-ref x`
    ("Lines total/split/joined/realigned/mismatch_removed/dup_removed/skipped: ...")
  - mismatch rate = mismatch_removed / total; above --max-ref-mismatch -> error
    (the input was called on another reference or build)

Exit status 1 on any error; messages are prefixed ERROR/WARNING/INFO on stderr.

Tools considered: `bcftools norm --check-ref` and `bcftools +fixref` check REF
alleles (the first is used); `bcftools reheader --fai` and `bcftools annotate
--rename-chrs` apply the fixes. No tool decides which contig names correspond and
whether a length difference is a build mismatch, so this logic is in this script.
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

# names of the mitochondrial contig in common references
MITO = {"chrM", "MT", "M", "chrMT"}
# ##contig=<ID=chr1,length=248956422,...>
CONTIG_LINE = re.compile(r"^##contig=<(.*)>$")
# summary line of bcftools norm
NORM_LINE = re.compile(
    r"total/split/joined/realigned/mismatch_removed/dup_removed/skipped:\s*([\d/]+)"
)


class Report:
    """Collects key/value results, warnings and errors of one check."""

    def __init__(self, name):
        self.name = name
        self.values = {}
        self.errors = []
        self.warnings = []

    def set(self, key, value):
        self.values[key] = value

    def write(self, path):
        # status: error > fixed (anything was renamed/rewritten/dropped) > ok
        if self.errors:
            status = "error"
        elif self.warnings or self.values.get("renamed_contigs", "none") != "none":
            status = "fixed"
        else:
            status = "ok"
        with open(path, "w") as fh:
            fh.write("key\tvalue\n")
            fh.write(f"status\t{status}\n")
            for key, value in self.values.items():
                fh.write(f"{key}\t{value}\n")
            for msg in self.warnings:
                fh.write(f"warning\t{msg}\n")
            for msg in self.errors:
                fh.write(f"error\t{msg}\n")
        return status

    def print_and_exit(self):
        # print all messages; exit 1 if there was an error
        for msg in self.warnings:
            print(f"WARNING: {self.name}: {msg}", file=sys.stderr)
        for msg in self.errors:
            print(f"ERROR: {self.name}: {msg}", file=sys.stderr)
        if self.errors:
            sys.exit(1)


def read_fai(path):
    """{contig: length} of the reference, in reference order."""
    lengths = {}
    with open(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[1].isdigit():
                sys.exit(f"ERROR: {path}: line {line_no}: not a FASTA index (.fai) line")
            lengths[fields[0]] = int(fields[1])
    if not lengths:
        sys.exit(f"ERROR: {path}: the reference index is empty")
    return lengths


def map_contig(name, reference):
    """Reference contig name for `name` (exact or safe rename), or None."""
    if name in reference:
        return name
    candidates = []
    # mitochondrial contig: any of the common names
    if name in MITO:
        candidates += ["chrM", "MT", "M", "chrMT"]
    # add or remove the 'chr' prefix
    candidates.append(name[3:] if name.startswith("chr") else "chr" + name)
    for candidate in candidates:
        if candidate in reference:
            return candidate
    return None


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def check_lengths(lengths, mapping, reference, report):
    """Compare contig lengths of the input with the reference."""
    fixes = []
    for name, length in lengths.items():
        ref_name = mapping.get(name)
        # contigs not in the reference, or without length, can't be compared
        if ref_name is None or length is None or length == reference[ref_name]:
            continue
        if ref_name in MITO:
            # header-only problem: the header is rewritten from the .fai
            fixes.append(f"{name}:{length}->{reference[ref_name]}")
            report.warnings.append(
                f"contig {name} has length {length}, the reference {ref_name} has "
                f"{reference[ref_name]}; the header is rewritten from the reference"
            )
        else:
            report.errors.append(
                f"contig {name} has length {length}, but the reference {ref_name} has "
                f"{reference[ref_name]}: the input is on another genome build than --fasta"
            )
    report.set("header_length_fixes", ",".join(fixes) or "none")


def write_mapping(names, reference, report, prefix):
    """Map all input contig names; write the rename map; return (mapping, unknown names)."""
    mapping = {name: map_contig(name, reference) for name in names}
    renames = {old: new for old, new in mapping.items() if new is not None and old != new}
    unknown = sorted(name for name, new in mapping.items() if new is None)
    # two input names for one reference contig (e.g. '1' and 'chr1') can't both be renamed
    targets = {}
    for old, new in mapping.items():
        if new is not None:
            targets.setdefault(new, []).append(old)
    for new, olds in sorted(targets.items()):
        if len(olds) > 1:
            report.errors.append(
                f"contigs {', '.join(sorted(olds))} all correspond to the reference "
                f"contig {new}; the input mixes naming styles"
            )
    # nothing on the reference at all: wrong reference or wrong file
    if names and len(unknown) == len(mapping):
        report.errors.append(
            f"none of the {len(names)} contigs ({', '.join(sorted(names)[:5])}"
            f"{', ...' if len(names) > 5 else ''}) is in the reference; is the "
            "input on another genome?"
        )
    with open(f"{prefix}.rename_chrs.tsv", "w") as fh:
        fh.writelines(f"{old}\t{new}\n" for old, new in sorted(renames.items()))
    report.set(
        "renamed_contigs", ",".join(f"{o}->{n}" for o, n in sorted(renames.items())) or "none"
    )
    # reference contigs present in the input, in reference order
    keep = [c for c in reference if c in set(mapping.values())]
    with open(f"{prefix}.keep_contigs.txt", "w") as fh:
        fh.writelines(f"{c}\n" for c in keep)
    return mapping, unknown


def parse_vcf_header(path):
    """({contig: length or None}, [samples]) from the output of `bcftools view -h`."""
    contigs, samples, has_chrom_line = {}, [], False
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            match = CONTIG_LINE.match(line)
            if match:
                # key=value pairs of the contig line (values may contain no commas here)
                fields = dict(kv.split("=", 1) for kv in match.group(1).split(",") if "=" in kv)
                if "ID" not in fields:
                    sys.exit(f"ERROR: {path}: contig line without ID: {line}")
                length = fields.get("length")
                contigs[fields["ID"]] = int(length) if length and length.isdigit() else None
            elif line.startswith("#CHROM"):
                has_chrom_line = True
                # sample names follow the 9 fixed columns
                samples = line.split("\t")[9:]
    if not has_chrom_line:
        sys.exit(f"ERROR: {path}: no #CHROM line; is this the output of `bcftools view -h`?")
    return contigs, samples


def count_record_contigs(path):
    """{contig: number of records} from the output of `bcftools query -f '%CHROM\\n'`."""
    counts = {}
    with open_text(path) as fh:
        for line in fh:
            chrom = line.strip()
            if chrom:
                counts[chrom] = counts.get(chrom, 0) + 1
    return counts


def contigs_vcf(args, reference, report):
    if not args.vcf_header or not args.record_contigs:
        sys.exit("ERROR: --type vcf needs --vcf-header and --record-contigs")
    header_lengths, samples = parse_vcf_header(args.vcf_header)
    record_contigs = count_record_contigs(args.record_contigs)
    n_records = sum(record_contigs.values())
    report.set("header_contigs", len(header_lengths))
    report.set("records", n_records)
    report.set("samples", ",".join(samples) or "none")
    # the samplesheet's tumor/normal sample columns must exist in the VCF
    if args.samples:
        missing = [s for s in args.samples.split(",") if s not in samples]
        if missing:
            report.errors.append(
                f"sample(s) {', '.join(missing)} not in the VCF (samples in the VCF: "
                f"{', '.join(samples) or 'none'}); check tumor_id/normal_id in the samplesheet"
            )
    mapping, unknown = write_mapping(
        set(header_lengths) | set(record_contigs), reference, report, args.prefix
    )
    check_lengths(header_lengths, mapping, reference, report)
    # records on contigs the reference doesn't have (decoys, alt contigs, ...) are dropped later
    dropped = {c: record_contigs[c] for c in unknown if c in record_contigs}
    report.set("records_on_unknown_contigs", sum(dropped.values()))
    if dropped:
        report.warnings.append(
            f"{sum(dropped.values())} records on contigs absent from the reference will be "
            f"dropped: {', '.join(f'{c} ({n})' for c, n in sorted(dropped.items()))}"
        )
    if n_records == 0:
        report.warnings.append("the VCF has no records")
    elif sum(dropped.values()) == n_records:
        report.errors.append(
            "no record is on a contig of the reference: wrong reference or contig naming "
            f"(record contigs: {', '.join(sorted(record_contigs)[:10])})"
        )


def contigs_bed(args, reference, report):
    if not args.bed:
        sys.exit("ERROR: --type bed needs --bed")
    intervals = []
    with open_text(args.bed) as fh:
        for line_no, line in enumerate(fh, start=1):
            # skip empty, comment and UCSC header lines
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                sys.exit(
                    f"ERROR: {args.bed}: line {line_no}: expected at least 3 tab-separated columns"
                )
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                sys.exit(
                    f"ERROR: {args.bed}: line {line_no}: start/end '{fields[1]}'/'{fields[2]}' are not integers"
                )
            if start < 0 or end < start:
                sys.exit(f"ERROR: {args.bed}: line {line_no}: invalid interval {start}-{end}")
            intervals.append((fields[0], start, end, fields[3:]))
    report.set("intervals", len(intervals))
    if not intervals:
        report.errors.append("the BED has no intervals")
        return
    mapping, unknown = write_mapping({c for c, *_ in intervals}, reference, report, args.prefix)
    dropped = sum(1 for c, *_ in intervals if mapping[c] is None)
    report.set("intervals_on_unknown_contigs", dropped)
    if dropped:
        report.warnings.append(
            f"{dropped} intervals on contigs absent from the reference are dropped: {', '.join(unknown)}"
        )
    if dropped == len(intervals):
        report.errors.append(
            "no interval is on a contig of the reference: wrong reference or contig naming"
        )
    # an interval beyond the end of its contig means another genome build
    beyond = [(c, e) for c, _, e, _ in intervals if mapping[c] and e > reference[mapping[c]]]
    report.set("intervals_beyond_contig_end", len(beyond))
    if beyond:
        c, e = beyond[0]
        report.errors.append(
            f"{len(beyond)} intervals end beyond their contig (e.g. {c} at {e}, reference length "
            f"{reference[mapping[c]]}): the BED is on another genome build than --fasta"
        )
    # harmonized BED: all columns, renamed contigs, unknown contigs dropped
    with open(f"{args.prefix}.harmonized.bed", "w") as fh:
        for c, s, e, rest in intervals:
            if mapping[c]:
                fh.write("\t".join([mapping[c], str(s), str(e), *rest]) + "\n")


def contigs_bam(args, reference, report):
    if not args.idxstats:
        sys.exit("ERROR: --type bam needs --idxstats")
    lengths = {}
    with open(args.idxstats) as fh:
        for line_no, line in enumerate(fh, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                sys.exit(
                    f"ERROR: {args.idxstats}: line {line_no}: expected 4 columns of samtools idxstats"
                )
            # the last line ('*') counts unplaced reads; it is no contig
            if fields[0] != "*":
                lengths[fields[0]] = int(fields[1])
    report.set("header_contigs", len(lengths))
    mapping, _unknown = write_mapping(set(lengths), reference, report, args.prefix)
    check_lengths(lengths, mapping, reference, report)
    if not any(mapping.values()):
        report.errors.append(
            "no contig of the alignment is in the reference: wrong reference or naming"
        )
    # translate a reference-named BED to the alignment's names (for mosdepth --by)
    if args.regions:
        to_alignment = {new: old for old, new in mapping.items() if new is not None}
        n_written = 0
        with open_text(args.regions) as fin, open(f"{args.prefix}.regions.bed", "w") as fout:
            for line in fin:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3 and fields[0] in to_alignment:
                    fout.write("\t".join([to_alignment[fields[0]], *fields[1:3]]) + "\n")
                    n_written += 1
        report.set("regions_translated", n_written)
        if n_written == 0:
            report.errors.append(
                f"none of the regions in {Path(args.regions).name} is on a contig of the alignment"
            )


def contigs(args):
    reference = read_fai(args.fai)
    report = Report(args.name)
    report.set("input", args.name)
    report.set("reference_contigs", len(reference))
    {"vcf": contigs_vcf, "bed": contigs_bed, "bam": contigs_bam}[args.type](args, reference, report)
    status = report.write(f"{args.prefix}.reference_check.tsv")
    print(f"INFO: {args.name}: reference check {status}", file=sys.stderr)
    report.print_and_exit()


def ref_alleles(args):
    report = Report(args.name)
    match = None
    with open(args.norm_log) as fh:
        for line in fh:
            match = NORM_LINE.search(line) or match
    if match is None:
        sys.exit(
            f"ERROR: {args.name}: {args.norm_log}: no summary line of bcftools norm found; "
            "the log must contain the stderr of `bcftools norm --check-ref x`"
        )
    counts = [int(x) for x in match.group(1).split("/")]
    if len(counts) != 7:
        sys.exit(f"ERROR: {args.name}: unexpected bcftools norm summary '{match.group(1)}'")
    total, _split, _joined, realigned, mismatch, dups, _skipped = counts
    rate = mismatch / total if total else 0.0
    report.set("records", total)
    report.set("ref_mismatches_removed", mismatch)
    report.set("ref_mismatch_rate", f"{rate:.6f}")
    report.set("realigned", realigned)
    report.set("duplicates_removed", dups)
    if mismatch:
        report.warnings.append(
            f"{mismatch} of {total} records had a REF allele different from the reference and were removed"
        )
    if rate > args.max_ref_mismatch:
        report.errors.append(
            f"{mismatch} of {total} REF alleles ({rate:.2%}) differ from the reference (threshold "
            f"--max_ref_mismatch {args.max_ref_mismatch:.2%}): the VCF was called on another "
            "reference or genome build than --fasta"
        )
    status = report.write(f"{args.prefix}.ref_alleles.tsv")
    print(
        f"INFO: {args.name}: REF allele check {status} ({mismatch} of {total} mismatching)",
        file=sys.stderr,
    )
    report.print_and_exit()


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    c = sub.add_parser("contigs", help="contig names, lengths and sample names")
    c.add_argument("--type", choices=["vcf", "bed", "bam"], required=True)
    c.add_argument("--name", required=True, help="name of the input for messages")
    c.add_argument("--fai", required=True, help="reference .fai")
    c.add_argument("--prefix", required=True, help="output prefix")
    c.add_argument("--vcf-header", help="(vcf) output of bcftools view -h")
    c.add_argument("--record-contigs", help="(vcf) output of bcftools query -f '%%CHROM\\n'")
    c.add_argument("--samples", help="(vcf) comma-separated sample names that must exist")
    c.add_argument("--bed", help="(bed) the BED file")
    c.add_argument("--idxstats", help="(bam) output of samtools idxstats")
    c.add_argument(
        "--regions", help="(bam) reference-named BED to translate to the alignment's names"
    )
    r = sub.add_parser("ref-alleles", help="REF allele mismatch rate from the bcftools norm log")
    r.add_argument("--name", required=True)
    r.add_argument("--norm-log", required=True, help="stderr of bcftools norm --check-ref x")
    r.add_argument("--max-ref-mismatch", type=float, default=0.01)
    r.add_argument("--prefix", required=True)
    args = parser.parse_args()
    contigs(args) if args.command == "contigs" else ref_alleles(args)


if __name__ == "__main__":
    main()
