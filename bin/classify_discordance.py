#!/usr/bin/env python3
"""Explain the discordant variants of one truth/query comparison.

Input is the per-variant table of benchmark_metrics.py (status TP, FN or FP).
For a variant missing in one technology ("other" = query for FN, truth for FP),
the script looks for the reason in this order:

  present_but_filtered   the other technology has the call in its normalized
                         VCF (all records), but it failed FILTER or the uniform
                         thresholds and therefore was not compared
  germline_in_normal     a normal alignment has alt reads (VAF >= --germline-vaf
                         and >= --min-alt-reads): likely a germline variant,
                         so the somatic call is expected to be missing
  low_coverage           the other technology's tumor depth < --min-coverage
  subthreshold_evidence  the other tumor alignment has >= 1 alt read, but the
                         variant wasn't called
  no_evidence            the other tumor alignment has no alt read
  absent                 the other technology has no call and no alignment was
                         given to look further

Read evidence comes from `bcftools mpileup -a AD,DP | bcftools norm -m-any`
VCFs of the alignments (one per alignment, matched by CHROM/POS/REF/ALT).

For truth variants (TP and FN), the script also computes the detection power of
the query: P(X >= --power-min-alt) with X ~ Binomial(query depth, truth VAF).
Variants with power >= --detection-power are "detectable"; the power-adjusted
recall only counts these. Without a query alignment, power is NA.

Outputs:
  <prefix>.discordance.tsv  input table + class and evidence columns
  <prefix>.power.tsv        power-adjusted recall (one metrics row)
"""

import argparse
import csv
import math

import pysam

EXTRA_COLUMNS = [
    "class", "other_in_normalized", "other_filter", "other_vaf_normalized",
    "other_tumor_dp", "other_tumor_alt", "normal_dp", "normal_alt", "normal_vaf",
    "query_power", "detectable",
]  # fmt: skip


def detection_power(depth, vaf, min_alt):
    """P(X >= min_alt) for X ~ Binomial(depth, vaf), exact (min_alt is small)."""
    if depth <= 0 or vaf <= 0:
        return 0.0
    vaf = min(vaf, 1.0)
    below = sum(
        math.comb(depth, i) * vaf**i * (1 - vaf) ** (depth - i)
        for i in range(min(min_alt, depth + 1))
    )
    return max(0.0, 1.0 - below)


def load_calls(path):
    """{(chrom, pos, ref, alt): (filter, vaf)} of a normalized VCF."""
    calls = {}
    if not path:
        return calls
    with pysam.VariantFile(path) as vcf:
        for rec in vcf:
            for alt in rec.alts or ():
                calls[(rec.chrom, rec.pos, rec.ref, alt)] = (
                    ";".join(rec.filter.keys()) or ".",
                    rec.info.get("VC_T_VAF"),
                )
    return calls


def load_evidence(path):
    """Pileup evidence: ({(chrom, pos, ref, alt): alt reads}, {(chrom, pos): depth})."""
    alts, depth = {}, {}
    if not path:
        return None
    with pysam.VariantFile(path) as vcf:
        sample = next(iter(vcf.header.samples))
        for rec in vcf:
            fmt = rec.samples[sample]
            ad = fmt.get("AD") if "AD" in rec.format else None
            dp = fmt.get("DP") if "DP" in rec.format else None
            if dp is None and ad is not None:
                dp = sum(x for x in ad if x is not None)
            if dp is not None:
                depth[(rec.chrom, rec.pos)] = max(depth.get((rec.chrom, rec.pos), 0), dp)
            for i, alt in enumerate(rec.alts or (), start=1):
                if alt.startswith("<") or ad is None or i >= len(ad) or ad[i] is None:
                    continue
                alts[(rec.chrom, rec.pos, rec.ref, alt)] = ad[i]
    return alts, depth


def evidence_at(evidence, key):
    """(depth, alt reads) at a variant; (None, None) without alignment."""
    if evidence is None:
        return None, None
    alts, depth = evidence
    return depth.get(key[:2], 0), alts.get(key, 0)


def merge_evidence(paths):
    """Sum normal evidence over several alignments (e.g. normals of truth and query)."""
    merged_alts, merged_depth, any_file = {}, {}, False
    for path in paths:
        loaded = load_evidence(path)
        if loaded is None:
            continue
        any_file = True
        for key, value in loaded[0].items():
            merged_alts[key] = merged_alts.get(key, 0) + value
        for key, value in loaded[1].items():
            merged_depth[key] = merged_depth.get(key, 0) + value
    return (merged_alts, merged_depth) if any_file else None


def fmt(value, digits=4):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def to_float(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--variants", required=True, help="benchmark_metrics.py variants table")
    parser.add_argument("--truth-normalized", required=True)
    parser.add_argument("--query-normalized", required=True)
    parser.add_argument("--truth-tumor-evidence", help="pileup VCF of the truth tumor alignment")
    parser.add_argument("--query-tumor-evidence", help="pileup VCF of the query tumor alignment")
    parser.add_argument("--normal-evidence", nargs="*", default=[], help="pileup VCFs of normals")
    parser.add_argument("--min-coverage", type=int, default=10)
    parser.add_argument("--min-alt-reads", type=int, default=2,
                        help="alt reads in a normal for germline_in_normal")  # fmt: skip
    parser.add_argument("--germline-vaf", type=float, default=0.1)
    parser.add_argument("--power-min-alt", type=int, default=3,
                        help="alt reads a caller needs to call a variant")  # fmt: skip
    parser.add_argument("--detection-power", type=float, default=0.95)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    truth_calls = load_calls(args.truth_normalized)
    query_calls = load_calls(args.query_normalized)
    truth_tumor = load_evidence(args.truth_tumor_evidence)
    query_tumor = load_evidence(args.query_tumor_evidence)
    normals = merge_evidence(args.normal_evidence)

    with open(args.variants) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    detectable_tp = detectable_fn = 0
    for row in rows:
        key = (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
        status = row["status"]
        other_calls, other_tumor = (
            (query_calls, query_tumor) if status in ("TP", "FN") else (truth_calls, truth_tumor)
        )
        in_other = key in other_calls
        other_filter, other_vaf = other_calls.get(key, ("NA", None))
        tumor_dp, tumor_alt = evidence_at(other_tumor, key)
        normal_dp, normal_alt = evidence_at(normals, key)
        normal_vaf = normal_alt / normal_dp if normal_dp else None

        row_class = "concordant"
        if status != "TP":
            if in_other:
                row_class = "present_but_filtered"
            elif (normal_alt or 0) >= args.min_alt_reads and (normal_vaf or 0) >= args.germline_vaf:
                row_class = "germline_in_normal"
            elif tumor_dp is None:
                row_class = "absent"
            elif tumor_dp < args.min_coverage:
                row_class = "low_coverage"
            elif tumor_alt > 0:
                row_class = "subthreshold_evidence"
            else:
                row_class = "no_evidence"

        power, detectable = None, None
        if status in ("TP", "FN"):
            query_dp, _ = evidence_at(query_tumor, key)
            truth_vaf = to_float(row.get("truth_vaf"))
            if query_dp is not None and truth_vaf is not None:
                power = detection_power(query_dp, truth_vaf, args.power_min_alt)
                detectable = power >= args.detection_power
                if detectable:
                    detectable_tp += status == "TP"
                    detectable_fn += status == "FN"

        row.update({
            "class": row_class,
            "other_in_normalized": "yes" if in_other else "no",
            "other_filter": other_filter,
            "other_vaf_normalized": fmt(other_vaf),
            "other_tumor_dp": fmt(tumor_dp),
            "other_tumor_alt": fmt(tumor_alt),
            "normal_dp": fmt(normal_dp),
            "normal_alt": fmt(normal_alt),
            "normal_vaf": fmt(normal_vaf),
            "query_power": fmt(power),
            "detectable": "NA" if detectable is None else ("yes" if detectable else "no"),
        })  # fmt: skip

    columns = (list(rows[0].keys()) if rows else []) or []
    for col in EXTRA_COLUMNS:
        if col not in columns:
            columns.append(col)
    with open(f"{args.prefix}.discordance.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    sample = rows[0]["sample"] if rows else "NA"
    truth = rows[0]["truth"] if rows else "NA"
    query = rows[0]["query"] if rows else "NA"
    has_power = query_tumor is not None
    total = detectable_tp + detectable_fn
    recall = detectable_tp / total if total else math.nan
    with open(f"{args.prefix}.power.tsv", "w") as fh:
        fh.write("sample\ttruth\tquery\tstratum_type\tstratum\ttp_baseline\tfn\ttp_call\tfp\t"
                 "recall\tprecision\tf1\n")  # fmt: skip
        fh.write("\t".join([sample, truth, query, "power_adjusted", "detectable",
                            str(detectable_tp) if has_power else "NA",
                            str(detectable_fn) if has_power else "NA", "NA", "NA",
                            fmt(recall) if has_power else "NA", "NA", "NA"]) + "\n")  # fmt: skip


if __name__ == "__main__":
    main()
