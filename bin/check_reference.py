#!/usr/bin/env python3
"""Check that an input (VCF, alignment, BED) is on the configured reference genome.

Runs at the start of the pipeline (subworkflow CHECK_INPUTS), so that a wrong input
stops the run before any long step. Files are read directly with pysam.

Subcommands, inputs and checks
  vcf   --vcf <VCF/BCF> --fai <.fai> [--fasta <FASTA>] [--samples T,N]
        - contig names: exact match with the reference or a safe rename
          (chr1 <-> 1, chrM <-> MT/M/chrMT); records on unknown contigs are dropped
          later and reported; no record on a reference contig -> error
        - contig lengths in the VCF header: any difference -> error (another build);
          lengths are never "fixed"
        - mitochondrial contig: dropped from the comparison and reported, its length is
          not checked (--keep-mito keeps it with the normal checks)
        - sample columns given with --samples must exist -> error otherwise
        - REF alleles (with --fasta, unless --skip-ref-alleles): REF of every record on a reference contig is
          compared with the FASTA; mismatch rate > --max-ref-mismatch -> error
  bam   --alignment <BAM/CRAM> --fai <.fai> [--fasta <FASTA> for CRAM]
        - contig names and lengths of the alignment header (same rules as VCF)
  bed   --bed <BED> --fai <.fai>
        - BED format (read_bed: tab-separated, integer 0 <= start <= end)
        - contig names (renamed / unknown dropped / mito dropped)
        - intervals ending beyond their contig -> error (another build)
  rename-bed --bed <reference-named BED> --rename <rename map> --output <BED>
        - translates a BED to the input's contig names (reverse of the rename map),
          e.g. the common region for mosdepth --by on an alignment with other names

--warn-only turns all errors into warnings (the user switched the check off);
the outputs are written in any case.

Outputs (--prefix)
  <prefix>.rename_chrs.tsv       two columns, no header: input contig name, reference
                                 contig name (only names that differ; for
                                 `bcftools annotate --rename-chrs` and mosdepth_callable.py)
  <prefix>.reference_check.tsv   two columns with header `key value`: status (ok,
                                 fixed = something was renamed/dropped, error), then
                                 counts and lists (see the docs, docs/output.md), one
                                 `warning` / `error` row per message
  <prefix>.harmonized.bed        (bed) all columns of the input BED, contigs renamed,
                                 unknown and mitochondrial contigs removed

Exit status 1 on any error (unless --warn-only); messages are prefixed
ERROR/WARNING/INFO on stderr.

Why a script: `bcftools norm --check-ref` and `bcftools +fixref` check REF alleles,
`bcftools reheader --fai` and `bcftools annotate --rename-chrs` apply renames, but no
tool decides which contig names correspond, whether a length difference means another
build, or checks all inputs of a run before it starts.
"""

import argparse
import gzip
import sys
from typing import NamedTuple

# names of the mitochondrial contig in common references (GRCh38/UCSC, Ensembl, ...)
MITO = ("chrM", "MT", "M", "chrMT")


class BedInterval(NamedTuple):
    """One BED line: contig, 0-based start, exclusive end, further columns, line number."""

    chrom: str
    start: int
    end: int
    extra: list
    line_no: int


class Report:
    """Collects key/value results, warnings and errors of one check."""

    def __init__(self, name, warn_only=False):
        self.name = name
        self.warn_only = warn_only
        self.values = {}
        self.errors = []
        self.warnings = []

    def set(self, key, value):
        """Store one result value (written as a `key value` row)."""
        self.values[key] = value

    def error(self, message):
        """Record an error; with --warn-only it is recorded as a warning instead."""
        if self.warn_only:
            self.warnings.append(f"(check switched off) {message}")
        else:
            self.errors.append(message)

    def write(self, path):
        """Write the report; return the status (error > fixed > ok)."""
        renamed = self.values.get("renamed_contigs", "none") != "none"
        dropped = any(
            str(self.values.get(key, 0)) not in ("0", "none")
            for key in ("dropped_mito", "records_on_unknown_contigs", "intervals_dropped")
        )
        if self.errors:
            status = "error"
        elif self.warnings or renamed or dropped:
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
        """Print all messages to stderr; exit 1 if there was an error."""
        for msg in self.warnings:
            print(f"WARNING: {self.name}: {msg}", file=sys.stderr)
        for msg in self.errors:
            print(f"ERROR: {self.name}: {msg}", file=sys.stderr)
        if self.errors:
            sys.exit(1)


# ------------------------------------------------------------------ readers


def open_text(path):
    """Text handle of a plain or gzipped file."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fai(path):
    """{contig: length} of a FASTA index (.fai), in reference order; strict."""
    lengths = {}
    with open(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            fields = line.rstrip("\n").split("\t")
            # a .fai line starts with the contig name and its length
            if len(fields) < 2 or not fields[1].isdigit():
                sys.exit(f"ERROR: {path}: line {line_no}: not a FASTA index (.fai) line")
            if fields[0] in lengths:
                sys.exit(f"ERROR: {path}: line {line_no}: contig {fields[0]} occurs twice")
            lengths[fields[0]] = int(fields[1])
    if not lengths:
        sys.exit(f"ERROR: {path}: the reference index is empty")
    return lengths


def read_bed(path):
    """All intervals of a BED file; stops with a line-specific error if it is malformed.

    Header lines ('#', 'track', 'browser') and empty lines are skipped. Columns: contig,
    start (0-based, integer >= 0), end (exclusive, integer >= start), further columns.
    """
    intervals = []
    with open_text(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            # skip empty, comment and UCSC header lines
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            # BED is tab-separated: a space-separated line has a single column
            if len(fields) < 3:
                hint = " (space-separated? BED needs tabs)" if len(line.split()) >= 3 else ""
                sys.exit(
                    f"ERROR: {path}: line {line_no}: expected at least 3 tab-separated "
                    f"columns{hint}"
                )
            # start and end are integers, 0-based half-open: 0 <= start <= end
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                sys.exit(
                    f"ERROR: {path}: line {line_no}: start/end '{fields[1]}'/'{fields[2]}' "
                    "are not integers"
                )
            if start < 0 or end < start:
                sys.exit(
                    f"ERROR: {path}: line {line_no}: invalid interval {start}-{end} "
                    "(need 0 <= start <= end)"
                )
            intervals.append(BedInterval(fields[0], start, end, fields[3:], line_no))
    return intervals


def open_vcf(path):
    """pysam.VariantFile of a VCF/BCF; stops with a clear error if it can't be read."""
    import pysam

    try:
        return pysam.VariantFile(str(path))
    except (OSError, ValueError) as err:
        sys.exit(f"ERROR: {path}: can't be read as VCF/BCF ({err})")


def open_alignment(path, fasta=None):
    """pysam.AlignmentFile of a BAM/CRAM (a CRAM needs the reference FASTA)."""
    import pysam

    try:
        return pysam.AlignmentFile(str(path), reference_filename=fasta)
    except (OSError, ValueError) as err:
        sys.exit(f"ERROR: {path}: can't be read as BAM/CRAM ({err})")


def read_rename_map(path):
    """{input name: reference name} of a rename map (two columns, may be empty)."""
    mapping = {}
    with open(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                sys.exit(f"ERROR: {path}: line {line_no}: expected two tab-separated columns")
            mapping[fields[0]] = fields[1]
    return mapping


# ------------------------------------------------------------------ contig decisions


def is_mito(name):
    """True for the usual mitochondrial contig names."""
    return name in MITO


def map_contig(name, reference):
    """Reference contig name for `name` (exact or safe rename), or None if there is none."""
    # the same name exists in the reference
    if name in reference:
        return name
    candidates = []
    # mitochondrial contig: any of the common names
    if is_mito(name):
        candidates += list(MITO)
    # add or remove the 'chr' prefix (1 <-> chr1)
    candidates.append(name[3:] if name.startswith("chr") else "chr" + name)
    for candidate in candidates:
        if candidate in reference:
            return candidate
    return None


def build_mapping(names, reference, report, keep_mito):
    """Map input contig names to the reference.

    Returns {input name: reference name or None}; None means the contig is not compared
    (absent from the reference, or mitochondrial without --keep-mito). Errors: two
    input names for one reference contig; no input contig on the reference.
    """
    mapping = {name: map_contig(name, reference) for name in names}
    # mitochondrial contigs are not compared (unless --keep-mito)
    mito = sorted(name for name, ref in mapping.items() if ref and is_mito(ref))
    if not keep_mito:
        for name in mito:
            mapping[name] = None
    report.set("dropped_mito", ",".join(mito) if mito and not keep_mito else "none")
    # two input names for one reference contig (e.g. '1' and 'chr1') can't both be renamed
    targets = {}
    for old, new in mapping.items():
        if new is not None:
            targets.setdefault(new, []).append(old)
    for new, olds in sorted(targets.items()):
        if len(olds) > 1:
            report.error(
                f"contigs {', '.join(sorted(olds))} all correspond to the reference contig "
                f"{new}; the input mixes naming styles"
            )
    # nothing on the reference at all: wrong reference or wrong file
    unknown = sorted(n for n, ref in mapping.items() if ref is None and n not in mito)
    if names and len(unknown) == len(names):
        shown = ", ".join(sorted(names)[:5]) + (", ..." if len(names) > 5 else "")
        report.error(
            f"none of the {len(names)} contigs ({shown}) is in the reference; is the input "
            "on another genome?"
        )
    renames = {old: new for old, new in mapping.items() if new is not None and old != new}
    report.set(
        "renamed_contigs", ",".join(f"{o}->{n}" for o, n in sorted(renames.items())) or "none"
    )
    report.set("unknown_contigs", ",".join(unknown) or "none")
    return mapping


def write_rename_map(mapping, prefix):
    """Write the contigs whose name differs from the reference (input<TAB>reference)."""
    with open(f"{prefix}.rename_chrs.tsv", "w") as fh:
        for old, new in sorted(mapping.items()):
            if new is not None and old != new:
                fh.write(f"{old}\t{new}\n")


def check_lengths(lengths, mapping, reference, report):
    """Compare contig lengths of the input with the reference; any difference is an error.

    Contigs without a length (VCF contig lines without length=) and contigs that are not
    compared (mapping None) are skipped. Lengths are never rewritten.
    """
    differ = []
    for name, length in sorted(lengths.items()):
        ref_name = mapping.get(name)
        if ref_name is None or length is None or length == reference[ref_name]:
            continue
        differ.append(f"{name}:{length} (reference {ref_name}:{reference[ref_name]})")
    report.set("length_mismatches", ",".join(differ) or "none")
    if differ:
        report.error(
            f"{len(differ)} contig(s) have another length than in the reference: "
            f"{', '.join(differ[:5])}; the input is on another genome build than --fasta"
        )


# ------------------------------------------------------------------ subcommands


def ref_alleles(vcf, mapping, fasta_path, report, max_mismatch):
    """Compare the REF allele of every compared VCF record with the reference FASTA.

    Records on contigs that are not compared are skipped. REF alleles with characters
    other than A/C/G/T/N are counted as mismatches. Reports records checked, mismatches
    and the rate; a rate above `max_mismatch` is an error (another reference/build).
    """
    import pysam

    fasta = pysam.FastaFile(fasta_path)
    checked = mismatches = 0
    examples = []
    for rec in vcf:
        ref_chrom = mapping.get(rec.chrom)
        if ref_chrom is None:
            continue
        checked += 1
        # the reference sequence at the REF position (0-based start, REF length)
        expected = fasta.fetch(ref_chrom, rec.start, rec.start + len(rec.ref)).upper()
        if rec.ref.upper() != expected:
            mismatches += 1
            if len(examples) < 3:
                examples.append(f"{rec.chrom}:{rec.pos} REF {rec.ref}, reference {expected}")
    rate = mismatches / checked if checked else 0.0
    report.set("ref_alleles_checked", checked)
    report.set("ref_mismatches", mismatches)
    report.set("ref_mismatch_rate", f"{rate:.6f}")
    if mismatches:
        report.warnings.append(
            f"{mismatches} of {checked} REF alleles differ from the reference (e.g. "
            f"{'; '.join(examples)}); these records are removed by bcftools norm"
        )
    if rate > max_mismatch:
        report.error(
            f"{mismatches} of {checked} REF alleles ({rate:.2%}) differ from the reference "
            f"(threshold --max_ref_mismatch {max_mismatch:.2%}): the VCF was called on another "
            "reference or genome build than --fasta"
        )


def check_vcf(args, reference, report):
    """Contigs, lengths, samples and (with --fasta) REF alleles of a VCF."""
    vcf = open_vcf(args.vcf)
    header_lengths = {name: contig.length for name, contig in vcf.header.contigs.items()}
    samples = list(vcf.header.samples)
    report.set("header_contigs", len(header_lengths))
    report.set("samples", ",".join(samples) or "none")
    # the samplesheet's tumor/normal sample columns must exist in the VCF
    if args.samples:
        missing = [s for s in args.samples.split(",") if s not in samples]
        if missing:
            report.error(
                f"sample(s) {', '.join(missing)} not in the VCF (samples in the VCF: "
                f"{', '.join(samples) or 'none'}); check tumor_id/normal_id in the samplesheet"
            )
    # contigs of all records (records may use contigs without header line)
    record_contigs = {}
    for rec in vcf:
        record_contigs[rec.chrom] = record_contigs.get(rec.chrom, 0) + 1
    n_records = sum(record_contigs.values())
    report.set("records", n_records)
    mapping = build_mapping(
        set(header_lengths) | set(record_contigs), reference, report, args.keep_mito
    )
    check_lengths(header_lengths, mapping, reference, report)
    write_rename_map(mapping, args.prefix)
    # records that are not compared: unknown contigs (warning) and mitochondrial ones
    unknown = {c: n for c, n in record_contigs.items() if mapping[c] is None and not is_mito(c)}
    mito = {c: n for c, n in record_contigs.items() if mapping[c] is None and is_mito(c)}
    report.set("records_on_unknown_contigs", sum(unknown.values()))
    report.set("records_on_mito", sum(mito.values()))
    if unknown:
        listed = ", ".join(f"{c} ({n})" for c, n in sorted(unknown.items()))
        report.warnings.append(
            f"{sum(unknown.values())} records on contigs absent from the reference are "
            f"dropped: {listed}"
        )
    if mito:
        report.warnings.append(
            f"{sum(mito.values())} records on the mitochondrial contig are not compared "
            "(--keep_mito false)"
        )
    if n_records == 0:
        report.warnings.append("the VCF has no records")
    elif sum(unknown.values()) + sum(mito.values()) == n_records:
        report.error(
            "no record is on a compared contig of the reference: wrong reference or contig "
            f"naming (record contigs: {', '.join(sorted(record_contigs)[:10])})"
        )
    # REF alleles against the FASTA (second pass over the records)
    if args.fasta and not args.skip_ref_alleles:
        vcf.close()
        ref_alleles(open_vcf(args.vcf), mapping, args.fasta, report, args.max_ref_mismatch)


def check_bam(args, reference, report):
    """Contig names and lengths of an alignment header."""
    alignment = open_alignment(args.alignment, args.fasta)
    lengths = dict(zip(alignment.references, alignment.lengths, strict=True))
    report.set("header_contigs", len(lengths))
    mapping = build_mapping(set(lengths), reference, report, args.keep_mito)
    check_lengths(lengths, mapping, reference, report)
    write_rename_map(mapping, args.prefix)


def check_bed(args, reference, report):
    """Format, contig names and interval ends of a BED; writes the harmonized BED."""
    intervals = read_bed(args.bed)
    report.set("intervals", len(intervals))
    if not intervals:
        report.error("the BED has no intervals")
    mapping = build_mapping({iv.chrom for iv in intervals}, reference, report, args.keep_mito)
    write_rename_map(mapping, args.prefix)
    # intervals that are not compared (unknown or mitochondrial contigs)
    dropped = [iv for iv in intervals if mapping[iv.chrom] is None]
    report.set("intervals_dropped", len(dropped))
    if dropped:
        contigs = sorted({iv.chrom for iv in dropped})
        report.warnings.append(
            f"{len(dropped)} intervals on contigs that are not compared are dropped: "
            f"{', '.join(contigs[:10])}"
        )
    if intervals and len(dropped) == len(intervals):
        report.error("no interval is on a compared contig of the reference")
    # an interval beyond the end of its contig means another genome build
    beyond = [iv for iv in intervals if mapping[iv.chrom] and iv.end > reference[mapping[iv.chrom]]]
    report.set("intervals_beyond_contig_end", len(beyond))
    if beyond:
        first = beyond[0]
        report.error(
            f"{len(beyond)} intervals end beyond their contig (e.g. line {first.line_no}: "
            f"{first.chrom} at {first.end}, reference length "
            f"{reference[mapping[first.chrom]]}): the BED is on another genome build than --fasta"
        )
    # harmonized BED: all columns, renamed contigs, dropped contigs removed
    with open(f"{args.prefix}.harmonized.bed", "w") as fh:
        for iv in intervals:
            if mapping[iv.chrom]:
                fh.write("\t".join([mapping[iv.chrom], str(iv.start), str(iv.end), *iv.extra]))
                fh.write("\n")


def rename_bed(args):
    """Translate a reference-named BED to an input's contig names (reverse rename map)."""
    to_input = {new: old for old, new in read_rename_map(args.rename).items()}
    written = 0
    with open(args.output, "w") as out:
        for iv in read_bed(args.bed):
            # contigs without rename keep their (reference) name
            out.write(f"{to_input.get(iv.chrom, iv.chrom)}\t{iv.start}\t{iv.end}\n")
            written += 1
    print(f"INFO: {args.bed}: {written} intervals written to {args.output}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    for name, help_text in (
        ("vcf", "contigs, lengths, samples and REF alleles of a VCF"),
        ("bam", "contigs and lengths of a BAM/CRAM header"),
        ("bed", "format, contigs and interval ends of a BED"),
    ):
        p = sub.add_parser(name, help=help_text)
        p.add_argument("--name", required=True, help="name of the input for messages")
        p.add_argument("--fai", required=True, help="reference .fai")
        p.add_argument("--prefix", required=True, help="output prefix")
        p.add_argument(
            "--keep-mito", action="store_true", help="compare the mitochondrial contig, too"
        )
        p.add_argument(
            "--warn-only", action="store_true", help="report errors as warnings (check off)"
        )
        if name == "vcf":
            p.add_argument("--vcf", required=True, help="VCF/BCF (plain, gzip or bgzip)")
            p.add_argument("--fasta", help="reference FASTA: check REF alleles")
            p.add_argument("--samples", help="comma-separated sample names that must exist")
            p.add_argument("--max-ref-mismatch", type=float, default=0.01)
            p.add_argument(
                "--skip-ref-alleles", action="store_true", help="don't compare REF alleles"
            )
        elif name == "bam":
            p.add_argument("--alignment", required=True, help="BAM or CRAM")
            p.add_argument("--fasta", help="reference FASTA (needed for CRAM)")
        else:
            p.add_argument("--bed", required=True, help="BED file")
    r = sub.add_parser("rename-bed", help="translate a BED to an input's contig names")
    r.add_argument("--bed", required=True, help="BED with reference contig names")
    r.add_argument("--rename", required=True, help="rename map of the input (input, reference)")
    r.add_argument("--output", required=True)
    args = parser.parse_args()

    if args.command == "rename-bed":
        rename_bed(args)
        return
    reference = read_fai(args.fai)
    report = Report(args.name, warn_only=args.warn_only)
    report.set("input", args.name)
    report.set("reference_contigs", len(reference))
    {"vcf": check_vcf, "bam": check_bam, "bed": check_bed}[args.command](args, reference, report)
    status = report.write(f"{args.prefix}.reference_check.tsv")
    print(f"INFO: {args.name}: reference check {status}", file=sys.stderr)
    report.print_and_exit()


if __name__ == "__main__":
    main()
