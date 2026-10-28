#!/usr/bin/env python3
"""Explain the discordant variants of one truth/query comparison.

Input is the per-variant table of benchmark_metrics.py (status TP, FN or FP).
A variant missing in one technology ("other" = the query for an FN, the truth for an FP)
gets reasons from the view of that other technology.

The main reason (column `class`) describes the call and the reads of the other
technology. Exactly one applies; they exclude each other:
  called, but not compared (the call is in its normalized VCF, same CHROM, POS, REF, ALT):
    filtered_non_pass      FILTER is not PASS (value in other_filter); with
                           --filter_pass_only only PASS calls are compared
    present_but_filtered   FILTER PASS, but the call failed the uniform thresholds
                           (--min_vaf, --min_depth, --min_alt_reads) or lies outside the
                           comparison region (the region is the same for both
                           technologies of a sample, so there is no separate
                           `outside_callable` reason)
  not called (no record with these alleles in the normalized VCF):
    absent                 no tumor alignment to look further (none given, or
                           --read_evidence false)
    low_coverage           tumor depth at the variant < --min-coverage
    subthreshold_evidence  depth >= --min-coverage and >= 1 read with the alt allele
    no_evidence            depth >= --min-coverage and no read with the alt allele
TP variants get `concordant`.

Further reasons that can come in addition to the main reason (column `reasons` = the
main reason first, then these, comma-separated):
    germline_in_normal     the normal alignments (summed) show the alt allele (alt reads
                           >= --min-alt-reads and VAF >= --germline-vaf): likely a
                           germline variant (or a contaminated normal)
Context reasons (indel, low VAF, repeat, GIAB and other stratifications) need the
annotation and are added by `annotate_variants.py merge --context-reasons`.

Read evidence comes from read_evidence.py (one VCF per alignment, one sample each, one
record per variant: FORMAT/DP, AD = REF and ALT reads, OT = reads with another allele):
reads and depth are matched by (CHROM, POS, REF, ALT); a variant without a record has
depth 0.
Indels are counted by their haplotype sequence, so the count doesn't depend on how the
aligner placed the indel.

For truth variants (TP and FN), the script also computes the detection power of
the query: P(X >= --power-min-alt) with X ~ Binomial(query tumor depth, truth VAF).
Variants with power >= --detection-power are "detectable"; the power-adjusted
recall only counts these. Without a query tumor alignment, power is NA.

With --no-read-evidence (pipeline --read_evidence false) no alignments are used:
the classes are concordant, present_but_filtered or absent, and the read-evidence
columns (other_tumor_dp, other_tumor_alt, normal_*, query_power, detectable) are not
written.

Outputs (--prefix):
  <prefix>.discordance.tsv  input table + class, reasons and evidence columns
                            (other_tumor_dp, other_tumor_alt, other_tumor_other: depth, alt
                            reads and reads with another allele of the other technology's
                            tumor alignment; normal_*: summed over the normal alignments)
  <prefix>.power.tsv        power-adjusted recall (one metrics row)

Why a script: the classification combines the normalized calls, the read evidence of
up to four alignments and the binomial power model; no tool does this.
"""

import argparse
import csv
import math
import sys
from collections import Counter

import pysam

REQUIRED_COLUMNS = ["sample", "truth", "query", "status", "chrom", "pos", "ref", "alt", "truth_vaf"]
# columns that need read evidence (left out with --no-read-evidence)
EVIDENCE_COLUMNS = [
    "other_tumor_dp",
    "other_tumor_alt",
    "other_tumor_other",
    "normal_dp",
    "normal_alt",
    "normal_vaf",
    "query_power",
    "detectable",
]
EXTRA_COLUMNS = [
    "class",
    "reasons",
    "other_in_normalized",
    "other_filter",
    "other_vaf_normalized",
    "other_tumor_dp",
    "other_tumor_alt",
    "other_tumor_other",
    "normal_dp",
    "normal_alt",
    "normal_vaf",
    "query_power",
    "detectable",
]
# the main reasons (exclusive) and the further reasons of this script
CALLED = ("filtered_non_pass", "present_but_filtered")
NOT_CALLED = ("absent", "low_coverage", "subthreshold_evidence", "no_evidence")


def detection_power(depth, vaf, min_alt):
    """P(X >= min_alt) for X ~ Binomial(depth, vaf), exact (min_alt is small)."""
    # no reads or no variant allele: the variant can't be seen
    if depth <= 0 or vaf <= 0:
        return 0.0
    vaf = min(vaf, 1.0)
    # P(X < min_alt) = sum of the binomial probabilities of 0 .. min_alt - 1 alt reads
    below = sum(
        math.comb(depth, i) * vaf**i * (1 - vaf) ** (depth - i)
        for i in range(min(min_alt, depth + 1))
    )
    return max(0.0, 1.0 - below)


def load_calls(path):
    """{(chrom, pos, ref, alt): (filter, vaf)} of a normalized (harmonized) VCF."""
    calls = {}
    with pysam.VariantFile(path) as vcf:
        # the harmonized VAF tag must be defined (output of harmonize_vcf.py)
        if "VC_T_VAF" not in vcf.header.info:
            sys.exit(
                f"ERROR: {path}: INFO/VC_T_VAF is not defined; expected a VCF harmonized "
                "by harmonize_vcf.py"
            )
        collisions = 0
        for rec in vcf:
            for alt in rec.alts or ():
                key = (rec.chrom, rec.pos, rec.ref, alt)
                # FILTER as written in the VCF ('PASS', or the failed filters joined by ';')
                call = (";".join(rec.filter.keys()) or ".", rec.info.get("VC_T_VAF"))
                if key in calls:
                    # the same allele twice (e.g. a PASS and a non-PASS record after
                    # atomizing): a PASS call wins, so the variant counts as called PASS
                    collisions += 1
                    if calls[key][0] in ("PASS", "."):
                        continue
                calls[key] = call
    print(
        f"INFO: {path}: {len(calls)} normalized calls ({collisions} alleles in more than one "
        "record; PASS kept)",
        file=sys.stderr,
    )
    return calls


def load_evidence(path):
    """Read evidence of one alignment (read_evidence.py): ({(chrom, pos, ref, alt): alt
    reads}, {(chrom, pos, ref, alt): depth}, {(chrom, pos, ref, alt): reads with another
    allele}). The depth is the variant's own (an indel counts the reads spanning its
    window), not the largest depth at the position."""
    alts, depth, others = {}, {}, {}
    with pysam.VariantFile(path) as vcf:
        # one alignment = one sample with FORMAT/AD (read_evidence.py)
        samples = list(vcf.header.samples)
        if len(samples) != 1:
            sys.exit(f"ERROR: {path}: expected one sample (one alignment), found {len(samples)}")
        if "AD" not in vcf.header.formats:
            sys.exit(f"ERROR: {path}: FORMAT/AD is not defined; expected read_evidence.py output")
        sample = samples[0]
        for rec in vcf:
            fmt = rec.samples[sample]
            ad = fmt.get("AD") if "AD" in rec.format else None
            dp = fmt.get("DP") if "DP" in rec.format else None
            # depth: FORMAT/DP, else the sum of the allelic depths
            if dp is None and ad is not None:
                dp = sum(x for x in ad if x is not None)
            # depth of this variant (one ALT per record in read_evidence.py)
            if dp is not None:
                for alt in rec.alts or ():
                    depth[(rec.chrom, rec.pos, rec.ref, alt)] = dp
            # alt reads of each non-symbolic ALT allele (AD index i = ALT i)
            for i, alt in enumerate(rec.alts or (), start=1):
                if alt.startswith("<") or ad is None or i >= len(ad) or ad[i] is None:
                    continue
                alts[(rec.chrom, rec.pos, rec.ref, alt)] = ad[i]
                # reads with another allele (FORMAT/OT of read_evidence.py, one ALT per record)
                if i == 1 and "OT" in rec.format and fmt.get("OT") is not None:
                    others[(rec.chrom, rec.pos, rec.ref, alt)] = fmt.get("OT")
    print(f"INFO: {path}: evidence at {len(depth)} positions", file=sys.stderr)
    return alts, depth, others


def evidence_at(evidence, key):
    """(depth, alt reads, other-allele reads) at a variant; 0 if the evidence has no record
    there; (None, None, None) without alignment."""
    if evidence is None:
        return None, None, None
    alts, depth, others = evidence
    return depth.get(key, 0), alts.get(key, 0), others.get(key, 0)


def merge_evidence(paths):
    """Sum the normal evidence over several alignments (normals of truth and query)."""
    if not paths:
        return None
    merged_alts, merged_depth, merged_others = Counter(), Counter(), Counter()
    for path in paths:
        alts, depth, others = load_evidence(path)
        merged_alts.update(alts)
        merged_depth.update(depth)
        merged_others.update(others)
    return dict(merged_alts), dict(merged_depth), dict(merged_others)


def fmt(value, digits=4):
    """Table value: NA for missing, fixed digits for floats."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def to_float(text, where):
    """Float of a table value; NA/empty -> None; anything else is an error."""
    if text in (None, "", "NA", "."):
        return None
    try:
        return float(text)
    except ValueError:
        sys.exit(f"ERROR: {where}: '{text}' is not a number")


def read_variants(path, sample, truth, query):
    """Rows of the variant table; checks columns, status and the comparison identity."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(
                f"ERROR: {path}: missing column(s) {', '.join(missing)}; expected the "
                "variants table of benchmark_metrics.py"
            )
        rows = list(reader)
    for i, row in enumerate(rows, start=2):
        # every row must belong to this comparison and have a valid status
        if (row["sample"], row["truth"], row["query"]) != (sample, truth, query):
            sys.exit(
                f"ERROR: {path}: line {i}: row of {row['sample']} {row['truth']} vs "
                f"{row['query']}, expected {sample} {truth} vs {query}"
            )
        if row["status"] not in ("TP", "FN", "FP"):
            sys.exit(f"ERROR: {path}: line {i}: status '{row['status']}' is not TP, FN or FP")
        if not row["pos"].isdigit():
            sys.exit(f"ERROR: {path}: line {i}: pos '{row['pos']}' is not a positive integer")
    return rows


def classify(status, in_other, other_filter, tumor_dp, tumor_alt, args):
    """Main reason of one variant (exactly one; see the module docstring)."""
    if status == "TP":
        return "concordant"
    # the other technology called it: the FILTER decides between the two called reasons
    if in_other:
        return "filtered_non_pass" if other_filter not in ("PASS", ".") else "present_but_filtered"
    # not called: what the other technology's tumor reads show
    if tumor_dp is None:
        return "absent"
    if tumor_dp < args.min_coverage:
        return "low_coverage"
    if tumor_alt > 0:
        return "subthreshold_evidence"
    return "no_evidence"


def reasons(main_reason, normal_alt, normal_vaf, args):
    """The reasons of one variant: the main reason first, then germline_in_normal if the
    normals show the alt allele; ['concordant'] for TP."""
    if main_reason == "concordant":
        return ["concordant"]
    found = [main_reason]
    # alt reads in the normal(s): likely germline (in addition to the main reason)
    if (normal_alt or 0) >= args.min_alt_reads and (normal_vaf or 0) >= args.germline_vaf:
        found.append("germline_in_normal")
    return found


def main():
    """Parse the command line and run the (sub)command."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--variants", required=True, help="benchmark_metrics.py variants table")
    parser.add_argument("--sample", required=True, help="sample of the comparison")
    parser.add_argument("--truth", required=True, help="truth technology")
    parser.add_argument("--query", required=True, help="query technology")
    parser.add_argument("--truth-normalized", required=True, help="normalized VCF of the truth")
    parser.add_argument("--query-normalized", required=True, help="normalized VCF of the query")
    parser.add_argument("--truth-tumor-evidence", help="pileup VCF of the truth tumor alignment")
    parser.add_argument("--query-tumor-evidence", help="pileup VCF of the query tumor alignment")
    parser.add_argument("--normal-evidence", nargs="*", default=[], help="pileup VCFs of normals")
    parser.add_argument(
        "--min-coverage", type=int, default=10, help="low_coverage below this depth"
    )
    parser.add_argument(
        "--min-alt-reads",
        type=int,
        default=2,
        help="alt reads in the normals for germline_in_normal",
    )
    parser.add_argument(
        "--germline-vaf", type=float, default=0.1, help="normal VAF for germline_in_normal"
    )
    parser.add_argument(
        "--power-min-alt", type=int, default=3, help="alt reads a caller needs to call a variant"
    )
    parser.add_argument(
        "--detection-power",
        type=float,
        default=0.95,
        help="power needed to count a truth variant as detectable",
    )
    parser.add_argument("--prefix", required=True)
    parser.add_argument(
        "--no-read-evidence",
        action="store_true",
        help="no alignments: classes without read evidence, no evidence columns",
    )
    args = parser.parse_args()
    # parameter ranges
    if args.min_coverage < 0 or args.min_alt_reads < 1 or args.power_min_alt < 1:
        parser.error("--min-coverage must be >= 0, --min-alt-reads and --power-min-alt >= 1")
    if not 0 <= args.germline_vaf <= 1 or not 0 <= args.detection_power <= 1:
        parser.error("--germline-vaf and --detection-power must be within [0, 1]")

    rows = read_variants(args.variants, args.sample, args.truth, args.query)
    truth_calls = load_calls(args.truth_normalized)
    query_calls = load_calls(args.query_normalized)
    # evidence files and --no-read-evidence exclude each other
    given = [args.truth_tumor_evidence, args.query_tumor_evidence, *args.normal_evidence]
    if args.no_read_evidence and any(given):
        sys.exit("ERROR: --no-read-evidence was given together with read-evidence VCFs")
    truth_tumor = load_evidence(args.truth_tumor_evidence) if args.truth_tumor_evidence else None
    query_tumor = load_evidence(args.query_tumor_evidence) if args.query_tumor_evidence else None
    normals = merge_evidence(args.normal_evidence)

    detectable_tp = detectable_fn = 0
    classes = Counter()
    for i, row in enumerate(rows, start=2):
        key = (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
        status = row["status"]
        # the technology that misses the variant: the query for TP/FN, the truth for FP
        other_calls, other_tumor = (
            (query_calls, query_tumor) if status in ("TP", "FN") else (truth_calls, truth_tumor)
        )
        in_other = key in other_calls
        other_filter, other_vaf = other_calls.get(key, ("NA", None))
        tumor_dp, tumor_alt, tumor_other = evidence_at(other_tumor, key)
        normal_dp, normal_alt, _ = evidence_at(normals, key)
        normal_vaf = normal_alt / normal_dp if normal_dp else None
        # main reason (exactly one), then the further reasons
        row_class = classify(status, in_other, other_filter, tumor_dp, tumor_alt, args)
        classes[(status, row_class)] += 1
        row_reasons = reasons(row_class, normal_alt, normal_vaf, args)

        # detection power of the query for truth variants (needs the query tumor depth)
        power, detectable = None, None
        if status in ("TP", "FN"):
            query_dp, _, _ = evidence_at(query_tumor, key)
            truth_vaf = to_float(row["truth_vaf"], f"{args.variants}: line {i}: truth_vaf")
            if query_dp is not None and truth_vaf is not None:
                power = detection_power(query_dp, truth_vaf, args.power_min_alt)
                detectable = power >= args.detection_power
                if detectable:
                    detectable_tp += status == "TP"
                    detectable_fn += status == "FN"

        row.update(
            {
                "class": row_class,
                "reasons": ",".join(row_reasons),
                "other_in_normalized": "yes" if in_other else "no",
                "other_filter": other_filter,
                "other_vaf_normalized": fmt(other_vaf),
                "other_tumor_dp": fmt(tumor_dp),
                "other_tumor_alt": fmt(tumor_alt),
                "other_tumor_other": fmt(tumor_other),
                "normal_dp": fmt(normal_dp),
                "normal_alt": fmt(normal_alt),
                "normal_vaf": fmt(normal_vaf),
                "query_power": fmt(power),
                "detectable": "NA" if detectable is None else ("yes" if detectable else "no"),
            }
        )

    # input columns first, then the new ones (an input that already has them keeps its order)
    with open(args.variants, newline="") as fh:
        columns = list(csv.DictReader(fh, delimiter="\t").fieldnames)
    columns += [c for c in EXTRA_COLUMNS if c not in columns]
    # without read evidence, its columns are not written at all
    if args.no_read_evidence:
        columns = [c for c in columns if c not in EVIDENCE_COLUMNS]
    with open(f"{args.prefix}.discordance.tsv", "w", newline="") as fh:
        # extrasaction='ignore': the evidence values are dropped with --no-read-evidence
        writer = csv.DictWriter(
            fh, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)

    # power-adjusted recall: NA without a query tumor alignment
    has_power = query_tumor is not None
    total = detectable_tp + detectable_fn
    recall = detectable_tp / total if total else math.nan
    with open(f"{args.prefix}.power.tsv", "w") as fh:
        fh.write(
            "sample\ttruth\tquery\tstratum_type\tstratum\ttp_baseline\tfn\ttp_call\tfp\t"
            "recall\tprecision\tf1\n"
        )
        fh.write(
            "\t".join(
                [
                    args.sample,
                    args.truth,
                    args.query,
                    "power_adjusted",
                    "detectable",
                    str(detectable_tp) if has_power else "NA",
                    str(detectable_fn) if has_power else "NA",
                    "NA",
                    "NA",
                    fmt(recall) if has_power else "NA",
                    "NA",
                    "NA",
                ]
            )
            + "\n"
        )

    # log the class counts per status
    summary = ", ".join(f"{s} {c}: {n}" for (s, c), n in sorted(classes.items()))
    print(
        f"INFO: {args.sample} {args.truth} vs {args.query}: {len(rows)} variants ({summary})",
        file=sys.stderr,
    )
    if has_power:
        print(
            f"INFO: detectable truth variants: {total} (TP {detectable_tp}, FN {detectable_fn})",
            file=sys.stderr,
        )
    else:
        print("INFO: no query tumor alignment: power-adjusted recall is NA", file=sys.stderr)


if __name__ == "__main__":
    main()
