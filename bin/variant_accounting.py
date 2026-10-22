#!/usr/bin/env python3
"""Accounting of one input VCF: how many records each preprocessing step kept or removed, and why.

The steps follow the normalization chain of PREPROCESS_VCF (docs/output.md):

  step  name           removed (reason)                          from
  1     input          -                                         --input-vcf
  2     pass_split     filter_not_pass (FILTER not PASS or '.')  --input-vcf (--pass-only)
  3     contigs        mitochondrial, contig_not_compared         --input-vcf + --rename +
                                                                  --genome-bed
  4     gvcf_blocks    gvcf_block (only ALT <NON_REF> or <*>)     --input-vcf
  5     normalization  ref_mismatch, duplicate; added: split      --norm-log (+ harmonize)
  6     harmonization  symbolic_allele                            --harmonize-summary
  7     filters        outside_region, vaf_below_min,             --harmonized + --regions +
                       depth_below_min, alt_reads_below_min        thresholds; --compare

Steps 2-4 are recomputed from the input records (the same rules as the bcftools steps);
step 7 is recomputed from the harmonized records with the uniform filters and the
comparison region (bcftools view --targets-overlap record): the records predicted to be
kept must be exactly the records of the compared VCF.

INPUT
  --input-vcf          the VCF of the samplesheet (VCF/BCF, plain/gzip/bgzip)
  --rename             rename map of the reference check (input contig, reference contig)
  --genome-bed         BED of the compared reference contigs (chrM excluded unless --keep_mito)
  --norm-log           bcftools norm log of the compared branch ('Lines total/split/...')
  --harmonize-summary  harmonize_vcf.py summary of the compared branch (key, value)
  --harmonized         harmonized VCF of the compared branch (before the uniform filters)
  --compare            compared VCF (after the uniform filters and the region)
  --regions            comparison region of the entry (BED, 0-based)
  --min-vaf, --min-depth, --min-alt-reads  the uniform filters (0 = off)
OUTPUT (--prefix)
  <prefix>.accounting.tsv   one row per step and outcome: sample, technology, step (1-7),
                            step_name, outcome (in, kept, removed, added), reason, records,
                            snv, indel, other (records by variant type: before normalization
                            from the alleles, after it from INFO/VC_TYPE; SNV = SNV + MNV,
                            indel = INS, DEL, COMPLEX; other = symbolic or mixed)
  <prefix>.removed.tsv      every record removed in steps 2-4 and 7: sample, technology,
                            step_name, reason, chrom, pos (1-based), ref, alt, filter
  <prefix>.checks.tsv       consistency checks: sample, technology, check, expected,
                            observed, status (ok / FAIL)

The script never stops on a failed check (the checks are evaluated at the end of the
pipeline, after the report); it stops on unreadable or inconsistent input files.

Why a script: bcftools reports how many records a step writes, but not why records were
removed (several reasons per step), and not whether the steps add up.
"""

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_reference import is_mito, read_bed  # noqa: E402

NORM_LINE = re.compile(
    r"total/split/joined/realigned/mismatch_removed/dup_removed/skipped:\s*([\d/]+)"
)
GVCF_ALTS = {"<NON_REF>", "<*>"}
TYPE_GROUP = {"SNV": "snv", "MNV": "snv", "INS": "indel", "DEL": "indel", "COMPLEX": "indel"}


def allele_type(ref, alts):
    """snv, indel or other of a record from its alleles (before normalization)."""
    alts = [a for a in alts or () if a is not None]
    if not alts or any(a.startswith("<") or a in ("*", ".") for a in alts):
        return "other"
    if all(len(a) == len(ref) for a in alts):
        return "snv"
    return "indel"


def harmonized_type(rec):
    """snv, indel or other of a harmonized record (INFO/VC_TYPE of harmonize_vcf.py)."""
    vc_type = rec.info.get("VC_TYPE")
    if isinstance(vc_type, tuple):
        vc_type = vc_type[0] if vc_type else None
    return TYPE_GROUP.get(vc_type, allele_type(rec.ref, rec.alts))


def read_rename(path):
    """{input contig: reference contig} of a rename map (empty file: no renaming)."""
    rename = {}
    for line_no, line in enumerate(Path(path).read_text().splitlines(), start=1):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            sys.exit(f"ERROR: {path}: line {line_no}: expected two tab-separated columns")
        rename[fields[0]] = fields[1]
    return rename


def read_norm_log(path):
    """{total, split, mismatch_removed, dup_removed} of a bcftools norm log."""
    match = None
    for line in Path(path).read_text().splitlines():
        match = NORM_LINE.search(line) or match
    if match is None:
        sys.exit(f"ERROR: {path}: no 'Lines total/split/...' line of bcftools norm found")
    keys = ["total", "split", "joined", "realigned", "mismatch_removed", "dup_removed", "skipped"]
    return dict(zip(keys, (int(v) for v in match.group(1).split("/")), strict=True))


def read_summary(path):
    """{key: value} of a harmonize_vcf.py summary (header 'key value')."""
    lines = Path(path).read_text().splitlines()
    if not lines or lines[0].split("\t") != ["key", "value"]:
        sys.exit(f"ERROR: {path}: expected a summary with the header 'key value'")
    return dict(line.split("\t", 1) for line in lines[1:] if "\t" in line)


class Counts:
    """Records per (step, outcome, reason), by variant type."""

    def __init__(self):
        self.rows = {}

    def add(self, step, name, outcome, reason, vtype, n=1):
        """Count n records of a type for a step, outcome and reason."""
        key = (step, name, outcome, reason)
        self.rows.setdefault(key, Counter())[vtype] += n

    def total(self, step, outcome, reason=None):
        """Records of a step and outcome (all reasons, or one)."""
        return sum(
            sum(c.values())
            for (s, _n, o, r), c in self.rows.items()
            if s == step and o == outcome and (reason is None or r == reason)
        )


def classify_input(args, counts, removed):
    """Steps 1-4 from the input records; returns the number of records entering norm."""
    rename = read_rename(args.rename)
    compared = {iv.chrom for iv in read_bed(args.genome_bed)}
    entering = 0
    with pysam.VariantFile(args.input_vcf) as vcf:
        for rec in vcf:
            vtype = allele_type(rec.ref, rec.alts)
            counts.add(1, "input", "in", "", vtype)
            filters = list(rec.filter.keys())
            row = (rec.chrom, rec.pos, rec.ref, ",".join(rec.alts or ()), ";".join(filters) or ".")
            # step 2: only PASS / '.' records are compared (--filter_pass_only)
            if args.pass_only and filters not in ([], ["PASS"]):
                counts.add(2, "pass_split", "removed", "filter_not_pass", vtype)
                removed.append(("pass_split", "filter_not_pass", *row))
                continue
            # step 3: contig renamed to the reference; only compared contigs are kept
            contig = rename.get(rec.chrom, rec.chrom)
            if contig not in compared:
                reason = "mitochondrial" if is_mito(contig) else "contig_not_compared"
                counts.add(3, "contigs", "removed", reason, vtype)
                removed.append(("contigs", reason, *row))
                continue
            # step 4: gVCF reference blocks (the only ALT is <NON_REF> or <*>)
            alts = rec.alts or ()
            if len(alts) == 1 and alts[0] in GVCF_ALTS:
                counts.add(4, "gvcf_blocks", "removed", "gvcf_block", vtype)
                removed.append(("gvcf_blocks", "gvcf_block", *row))
                continue
            entering += 1
    return entering


def classify_normalization(args, counts, entering, checks):
    """Steps 5-6 from the norm log and the harmonization summary; returns their output."""
    norm = read_norm_log(args.norm_log)
    summary = read_summary(args.harmonize_summary)
    norm_out = int(summary["records_read"])
    # records added by splitting multiallelic records / atomizing MNPs
    added = norm_out - norm["total"] + norm["mismatch_removed"] + norm["dup_removed"]
    counts.add(5, "normalization", "removed", "ref_mismatch", "all", norm["mismatch_removed"])
    counts.add(5, "normalization", "removed", "duplicate", "all", norm["dup_removed"])
    counts.add(5, "normalization", "added", "split", "all", added)
    counts.add(
        6,
        "harmonization",
        "removed",
        "symbolic_allele",
        "all",
        int(summary["records_dropped_symbolic"]),
    )
    checks.append(
        ("records entering bcftools norm = records left after steps 1-4", entering, norm["total"])
    )
    checks.append(
        (
            "records added by bcftools norm are not negative",
            "yes",
            "yes" if added >= 0 else f"no ({added})",
        )
    )
    return int(summary["records_written"])


def in_regions(index, chrom, start, end):
    """Does the 0-based half-open span [start, end) overlap an interval of the BED index?"""
    return any(s < end and e > start for s, e in index.get(chrom, ()))


def filter_reason(rec, args, regions):
    """First reason why a harmonized record is not compared, or None (kept)."""
    # --targets-overlap record: the REF span must overlap the comparison region
    if not in_regions(regions, rec.chrom, rec.pos - 1, rec.pos - 1 + len(rec.ref)):
        return "outside_region"
    # uniform filters (bcftools --include): a missing value fails an active filter
    for threshold, tag, reason in (
        (args.min_vaf, "VC_T_VAF", "vaf_below_min"),
        (args.min_depth, "VC_T_DP", "depth_below_min"),
        (args.min_alt_reads, "VC_T_ALT", "alt_reads_below_min"),
    ):
        if threshold > 0:
            value = rec.info.get(tag)
            if value is None or value < threshold:
                return reason
    return None


def classify_filters(args, counts, removed, checks):
    """Step 7: uniform filters and region, recomputed and compared with the compared VCF."""
    regions = {}
    for iv in read_bed(args.regions):
        regions.setdefault(iv.chrom, []).append((iv.start, iv.end))
    predicted, harmonized = set(), 0
    with pysam.VariantFile(args.harmonized) as vcf:
        for rec in vcf:
            harmonized += 1
            vtype = harmonized_type(rec)
            key = (rec.chrom, rec.pos, rec.ref, rec.alts)
            reason = filter_reason(rec, args, regions)
            if reason is None:
                predicted.add(key)
                counts.add(7, "filters", "kept", "", vtype)
                continue
            counts.add(7, "filters", "removed", reason, vtype)
            removed.append(
                (
                    "filters",
                    reason,
                    rec.chrom,
                    rec.pos,
                    rec.ref,
                    ",".join(rec.alts or ()),
                    ";".join(rec.filter.keys()) or ".",
                )
            )
    with pysam.VariantFile(args.compare) as vcf:
        compared = {(rec.chrom, rec.pos, rec.ref, rec.alts) for rec in vcf}
    checks.append(("records of the compared VCF", len(predicted), len(compared)))
    checks.append(
        (
            "compared VCF = harmonized records passing the filters and region",
            0,
            len(predicted ^ compared),
        )
    )
    return harmonized


def write_outputs(args, counts, removed, checks):
    """accounting.tsv, removed.tsv and checks.tsv."""
    with open(f"{args.prefix}.accounting.tsv", "w") as fh:
        fh.write(
            "sample\ttechnology\tstep\tstep_name\toutcome\treason\trecords\tsnv\tindel\tother\n"
        )
        for (step, name, outcome, reason), c in sorted(counts.rows.items()):
            total = sum(c.values())
            # counts without a type split (from logs) are written as NA per type
            per_type = ["NA"] * 3 if "all" in c else [str(c[t]) for t in ("snv", "indel", "other")]
            fh.write(
                "\t".join(
                    [
                        args.sample,
                        args.technology,
                        str(step),
                        name,
                        outcome,
                        reason or ".",
                        str(total),
                        *per_type,
                    ]
                )
                + "\n"
            )
    with open(f"{args.prefix}.removed.tsv", "w") as fh:
        fh.write("sample\ttechnology\tstep_name\treason\tchrom\tpos\tref\talt\tfilter\n")
        for row in removed:
            fh.write("\t".join([args.sample, args.technology, *map(str, row)]) + "\n")
    with open(f"{args.prefix}.checks.tsv", "w") as fh:
        fh.write("sample\ttechnology\tcheck\texpected\tobserved\tstatus\n")
        for name, expected, observed in checks:
            status = "ok" if str(expected) == str(observed) else "FAIL"
            fh.write(
                f"{args.sample}\t{args.technology}\t{name}\t{expected}\t{observed}\t{status}\n"
            )
            if status == "FAIL":
                print(
                    f"WARNING: {args.sample} {args.technology}: check failed: {name}: expected {expected}, observed {observed}",
                    file=sys.stderr,
                )


def main():
    """Parse the command line and account for every preprocessing step."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--sample", required=True)
    parser.add_argument("--technology", required=True)
    for name in (
        "input-vcf",
        "rename",
        "genome-bed",
        "norm-log",
        "harmonize-summary",
        "harmonized",
        "compare",
        "regions",
    ):
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument(
        "--pass-only", action="store_true", help="only PASS/'.' records are compared"
    )
    parser.add_argument("--min-vaf", type=float, default=0)
    parser.add_argument("--min-depth", type=float, default=0)
    parser.add_argument("--min-alt-reads", type=float, default=0)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    counts, removed, checks = Counts(), [], []
    entering = classify_input(args, counts, removed)
    harmonized_expected = classify_normalization(args, counts, entering, checks)
    harmonized = classify_filters(args, counts, removed, checks)
    checks.append(
        (
            "records of the harmonized VCF = harmonize summary records_written",
            harmonized_expected,
            harmonized,
        )
    )
    # step 1 = sum of the records removed in steps 2-4 + records entering norm
    removed_before_norm = sum(counts.total(s, "removed") for s in (2, 3, 4))
    checks.append(
        (
            "input records = removed in steps 2-4 + entering norm",
            counts.total(1, "in"),
            removed_before_norm + entering,
        )
    )
    write_outputs(args, counts, removed, checks)
    print(
        f"INFO: {args.sample} {args.technology}: {counts.total(1, 'in')} input records, "
        f"{entering} normalized, {harmonized} harmonized, {counts.total(7, 'kept')} compared; "
        f"{len(removed)} removed records listed",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
