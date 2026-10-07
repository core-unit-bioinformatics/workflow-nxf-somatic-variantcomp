#!/usr/bin/env python3
"""Recall, precision and F1 of one truth/query comparison.

Reads the output of `rtg vcfeval` (tp-baseline, fn, tp, fp) for a truth and a
query technology of the same sample, plus the two compared VCFs, and writes:

  <prefix>.metrics.tsv   one row per stratum
      stratum_type  stratum            counts used
      overall       all                TP_baseline, FN, TP_call, FP
      variant_type  SNV | INDEL        (INDEL = INS, DEL, COMPLEX; MNV counted as SNV)
      truth_vaf     [lo,hi)            TP_baseline, FN   -> recall only
      query_vaf     [lo,hi)            TP_call, FP       -> precision only
      match         exact | position   set comparison of the two VCFs (see below)
  <prefix>.variants.tsv  one row per variant: status (TP, FN, FP), type, length,
                         VAF/depth in truth and query (where called)

vcfeval matching is haplotype-aware. The `match` strata compare the two VCFs as
sets instead: `exact` uses (CHROM, POS, REF, ALT), `position` uses (CHROM, POS)
only. Their difference to `overall` shows how much discordance comes from variant
representation and how much from different alleles at the same position.

With --regions (one or more BEDs), the `match` strata only use records inside
all given regions (the same regions vcfeval was restricted to).

recall = TP_baseline / (TP_baseline + FN); precision = TP_call / (TP_call + FP);
F1 = harmonic mean. Undefined values are written as NA.
"""

import argparse
import csv
import math

import pysam

TYPE_GROUP = {"SNV": "SNV", "MNV": "SNV", "INS": "INDEL", "DEL": "INDEL", "COMPLEX": "INDEL"}
VARIANT_COLUMNS = [
    "sample", "truth", "query", "status", "chrom", "pos", "ref", "alt", "type", "length",
    "truth_vaf", "truth_dp", "truth_alt", "query_vaf", "query_dp", "query_alt", "query_filter",
]  # fmt: skip
METRIC_COLUMNS = [
    "sample", "truth", "query", "stratum_type", "stratum", "tp_baseline", "fn", "tp_call", "fp",
    "recall", "precision", "f1",
]  # fmt: skip


def ratio(num, den):
    return num / den if den else math.nan


def f1(recall, precision):
    if math.isnan(recall) or math.isnan(precision) or recall + precision == 0:
        return math.nan
    return 2 * recall * precision / (recall + precision)


def fmt(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def variant_type(ref, alt):
    if len(ref) == len(alt):
        return "SNV" if len(ref) == 1 else "MNV"
    if len(ref) == 1 and alt.startswith(ref):
        return "INS"
    if len(alt) == 1 and ref.startswith(alt):
        return "DEL"
    return "COMPLEX"


def read_records(path):
    """Yield dicts for every biallelic record of a VCF (empty if path is None)."""
    if not path:
        return
    with pysam.VariantFile(path) as vcf:
        tags = set(vcf.header.info)
        for rec in vcf:
            info = {t: rec.info.get(t) for t in ("VC_TYPE", "VC_LEN", "VC_T_VAF", "VC_T_DP",
                                                 "VC_T_ALT") if t in tags}  # fmt: skip
            for alt in rec.alts or ():
                vtype = info.get("VC_TYPE") or variant_type(rec.ref, alt)
                yield {
                    "chrom": rec.chrom,
                    "pos": rec.pos,
                    "ref": rec.ref,
                    "alt": alt,
                    "type": vtype,
                    "length": info.get("VC_LEN")
                    if info.get("VC_LEN") is not None
                    else len(alt) - len(rec.ref),
                    "vaf": info.get("VC_T_VAF"),
                    "dp": info.get("VC_T_DP"),
                    "alt_count": info.get("VC_T_ALT"),
                    "filter": ";".join(rec.filter.keys()) or ".",
                }


def load_regions(path):
    regions = {}
    with open(path) as fh:
        for line in fh:
            fields = line.split("\t")
            if len(fields) >= 3 and not line.startswith(("#", "track", "browser")):
                regions.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    return regions


def inside(record, region_sets):
    start, end = record["pos"] - 1, record["pos"] - 1 + len(record["ref"])
    return all(
        any(s < end and start < e for s, e in regions.get(record["chrom"], ()))
        for regions in region_sets
    )


def parse_bins(text):
    edges = [float(x) for x in text.split(",")]
    if edges != sorted(edges) or len(edges) < 2:
        raise SystemExit("--vaf-bins must be at least two increasing numbers")
    return list(zip(edges[:-1], edges[1:], strict=True))


def bin_label(lo, hi):
    return f"[{lo:g},{hi:g})"


def find_bin(value, bins):
    if value is None:
        return "unknown"
    for i, (lo, hi) in enumerate(bins):
        last = i == len(bins) - 1
        if lo <= value < hi or (last and value == hi):
            return bin_label(lo, hi)
    return "unknown"


def metric_row(keys, stratum_type, stratum, tpb, fn, tpc, fp):
    recall = ratio(tpb, tpb + fn) if tpb is not None else math.nan
    precision = ratio(tpc, tpc + fp) if tpc is not None else math.nan
    row = dict(keys, stratum_type=stratum_type, stratum=stratum, tp_baseline=tpb, fn=fn,
               tp_call=tpc, fp=fp, recall=recall, precision=precision,
               f1=f1(recall, precision))  # fmt: skip
    return {k: fmt(v) if k not in keys else v for k, v in row.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--tp-baseline", required=True)
    parser.add_argument("--fn", required=True)
    parser.add_argument("--tp", required=True)
    parser.add_argument("--fp", required=True)
    parser.add_argument("--truth-vcf", required=True, help="compared truth VCF")
    parser.add_argument("--query-vcf", required=True, help="compared query VCF")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--truth", required=True, help="truth technology")
    parser.add_argument("--query", required=True, help="query technology")
    parser.add_argument("--vaf-bins", default="0,0.05,0.1,0.2,0.5,1")
    parser.add_argument("--regions", nargs="*", default=[], help="restrict match strata (BEDs)")
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    keys = {"sample": args.sample, "truth": args.truth, "query": args.query}
    bins = parse_bins(args.vaf_bins)
    tpb = list(read_records(args.tp_baseline))
    fn = list(read_records(args.fn))
    tpc = list(read_records(args.tp))
    fp = list(read_records(args.fp))

    rows = [metric_row(keys, "overall", "all", len(tpb), len(fn), len(tpc), len(fp))]
    for group in ("SNV", "INDEL"):

        def count(records, group=group):
            return sum(TYPE_GROUP.get(r["type"], "INDEL") == group for r in records)

        rows.append(metric_row(keys, "variant_type", group, count(tpb), count(fn),
                               count(tpc), count(fp)))  # fmt: skip
    labels = [bin_label(lo, hi) for lo, hi in bins] + ["unknown"]
    for label in labels:

        def in_bin(records, label=label):
            return sum(find_bin(r["vaf"], bins) == label for r in records)

        if label != "unknown" or in_bin(tpb) + in_bin(fn):
            rows.append(metric_row(keys, "truth_vaf", label, in_bin(tpb), in_bin(fn), None, None))
        if label != "unknown" or in_bin(tpc) + in_bin(fp):
            rows.append(metric_row(keys, "query_vaf", label, None, None, in_bin(tpc), in_bin(fp)))

    truth_all = list(read_records(args.truth_vcf))
    query_all = list(read_records(args.query_vcf))
    if args.regions:
        region_sets = [load_regions(p) for p in args.regions]
        truth_all = [r for r in truth_all if inside(r, region_sets)]
        query_all = [r for r in query_all if inside(r, region_sets)]
    for name, key in (("exact", ("chrom", "pos", "ref", "alt")), ("position", ("chrom", "pos"))):
        truth_set = {tuple(r[k] for k in key) for r in truth_all}
        query_set = {tuple(r[k] for k in key) for r in query_all}
        shared = len(truth_set & query_set)
        rows.append(metric_row(keys, "match", name, shared, len(truth_set - query_set), shared,
                               len(query_set - truth_set)))  # fmt: skip

    with open(f"{args.prefix}.metrics.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=METRIC_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    # per-variant table: TP rows pair the baseline and call representation by position
    calls_by_pos = {}
    for r in tpc:
        calls_by_pos.setdefault((r["chrom"], r["pos"]), []).append(r)
    out = []
    for status, records in (("TP", tpb), ("FN", fn)):
        for r in records:
            call = None
            if status == "TP":
                candidates = calls_by_pos.get((r["chrom"], r["pos"]), [])
                call = next((c for c in candidates if c["alt"] == r["alt"]), None) or (
                    candidates[0] if candidates else None
                )
            out.append(dict(keys, status=status, chrom=r["chrom"], pos=r["pos"], ref=r["ref"],
                            alt=r["alt"], type=r["type"], length=r["length"],
                            truth_vaf=fmt(r["vaf"]), truth_dp=fmt(r["dp"]),
                            truth_alt=fmt(r["alt_count"]),
                            query_vaf=fmt(call["vaf"]) if call else "NA",
                            query_dp=fmt(call["dp"]) if call else "NA",
                            query_alt=fmt(call["alt_count"]) if call else "NA",
                            query_filter=call["filter"] if call else "NA"))  # fmt: skip
    for r in fp:
        out.append(dict(keys, status="FP", chrom=r["chrom"], pos=r["pos"], ref=r["ref"],
                        alt=r["alt"], type=r["type"], length=r["length"], truth_vaf="NA",
                        truth_dp="NA", truth_alt="NA", query_vaf=fmt(r["vaf"]),
                        query_dp=fmt(r["dp"]), query_alt=fmt(r["alt_count"]),
                        query_filter=r["filter"]))  # fmt: skip
    out.sort(key=lambda r: (r["chrom"], r["pos"], r["ref"], r["alt"], r["status"]))
    with open(f"{args.prefix}.variants.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=VARIANT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out)


if __name__ == "__main__":
    main()
