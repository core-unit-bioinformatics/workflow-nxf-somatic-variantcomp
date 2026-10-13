#!/usr/bin/env python3
"""Recall, precision and F1 of one truth/query comparison.

Inputs (all produced by tools before this script):
  --tp-baseline, --fn, --tp, --fp
        the four VCFs of `rtg vcfeval` (true positives in the baseline/truth
        representation, false negatives, true positives in the call/query
        representation, false positives)
  --isec-exact, --isec-position
        the output folders of `bcftools isec -c none` (same CHROM, POS, REF, ALT) and
        `bcftools isec -c all` (same CHROM, POS, any allele) of the truth and the query
        VCF; only their sites.txt is read (CHROM, POS, REF, ALT, presence flags
        "10" = truth only, "01" = query only, "11" = both)

Outputs (--prefix):
  <prefix>.metrics.tsv   one row per stratum
      stratum_type  stratum            counts used
      overall       all                TP_baseline, FN, TP_call, FP
      variant_type  SNV | INDEL        (INDEL = INS, DEL, COMPLEX; MNV counted as SNV)
      truth_vaf     [lo,hi)            TP_baseline, FN   -> recall only
      query_vaf     [lo,hi)            TP_call, FP       -> precision only
      match         exact | position   bcftools isec counts (see below)
  <prefix>.variants.tsv  one row per variant: status (TP, FN, FP), type, length,
                         VAF/depth/alt reads in truth and query (where called)

vcfeval matching is haplotype-aware, so different representations of the same
variant match. The `match` strata compare the two VCFs record by record instead:
`exact` by (CHROM, POS, REF, ALT), `position` by (CHROM, POS) only (counted per
position: several alleles at one position count once). Their
difference to `overall` shows how much discordance comes from variant
representation and how much from different alleles at the same position.
With --regions (one or more BEDs), the `match` strata only count sites inside
all given regions (the regions vcfeval was restricted to).

VAF, depth, alt reads and the variant type are taken from the harmonized INFO tags
of harmonize_vcf.py (VC_T_VAF, VC_T_DP, VC_T_ALT, VC_TYPE, VC_LEN). A truth set
without these tags (known-variant benchmark) is accepted: type and length are
derived from REF/ALT, VAF-based strata are 'unknown' (counted and reported).
A VCF with the header line `##variantcomp_vaf=unavailable` (too many records without
VAF, see harmonize_vcf.py) is treated the same way: all its VAFs count as unknown.

recall = TP_baseline / (TP_baseline + FN); precision = TP_call / (TP_call + FP);
F1 = harmonic mean. Undefined values are written as NA.
"""

import argparse
import csv
import math
import sys
from pathlib import Path

import pysam

# variant types of harmonize_vcf.py grouped for the variant_type stratum
TYPE_GROUP = {"SNV": "SNV", "MNV": "SNV", "INS": "INDEL", "DEL": "INDEL", "COMPLEX": "INDEL"}
VARIANT_COLUMNS = [
    "sample",
    "truth",
    "query",
    "status",
    "chrom",
    "pos",
    "ref",
    "alt",
    "type",
    "length",
    "truth_vaf",
    "truth_dp",
    "truth_alt",
    "query_vaf",
    "query_dp",
    "query_alt",
    "query_filter",
]
METRIC_COLUMNS = [
    "sample",
    "truth",
    "query",
    "stratum_type",
    "stratum",
    "tp_baseline",
    "fn",
    "tp_call",
    "fp",
    "recall",
    "precision",
    "f1",
]
HARMONIZED_TAGS = ("VC_TYPE", "VC_LEN", "VC_T_VAF", "VC_T_DP", "VC_T_ALT")


def ratio(num, den):
    """num / den, NaN if den is 0."""
    return num / den if den else math.nan


def f1(recall, precision):
    """Harmonic mean of recall and precision; NaN if undefined."""
    if math.isnan(recall) or math.isnan(precision) or recall + precision == 0:
        return math.nan
    return 2 * recall * precision / (recall + precision)


def fmt(value):
    """Table representation: NA for missing/NaN, 4 decimals for floats."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def variant_type(ref, alt):
    """SNV, MNV, INS, DEL or COMPLEX from the alleles (same rules as harmonize_vcf.py)."""
    if len(ref) == len(alt):
        return "SNV" if len(ref) == 1 else "MNV"
    if len(ref) == 1 and alt.startswith(ref):
        return "INS"
    if len(alt) == 1 and ref.startswith(alt):
        return "DEL"
    return "COMPLEX"


def read_records(path, what):
    """List of dicts (one per ALT allele) of a vcfeval output VCF; reports missing tags."""
    if not Path(path).exists():
        sys.exit(f"ERROR: {what}: {path} does not exist (vcfeval output missing?)")
    records = []
    with pysam.VariantFile(path) as vcf:
        # only tags declared in the header can be read (pysam raises otherwise)
        tags = set(vcf.header.info)
        missing_tags = [t for t in HARMONIZED_TAGS if t not in tags]
        # harmonize_vcf.py marks VCFs with too many records without VAF: no VAF strata
        vaf_unavailable = any(
            r.key == "variantcomp_vaf" and r.value == "unavailable" for r in vcf.header.records
        )
        if vaf_unavailable:
            print(
                f"WARNING: {what}: {path}: VAF unavailable (harmonize_vcf.py); its variants are "
                "counted in the VAF stratum 'unknown'",
                file=sys.stderr,
            )
        for rec in vcf:
            info = {t: rec.info.get(t) for t in HARMONIZED_TAGS if t in tags}
            for alt in rec.alts or ():
                records.append(
                    {
                        "chrom": rec.chrom,
                        "pos": rec.pos,
                        "ref": rec.ref,
                        "alt": alt,
                        # type and length from the harmonized tags, else from the alleles
                        "type": info.get("VC_TYPE") or variant_type(rec.ref, alt),
                        "length": (
                            info["VC_LEN"]
                            if info.get("VC_LEN") is not None
                            else len(alt) - len(rec.ref)
                        ),
                        "vaf": None if vaf_unavailable else info.get("VC_T_VAF"),
                        "dp": info.get("VC_T_DP"),
                        "alt_count": info.get("VC_T_ALT"),
                        "filter": ";".join(rec.filter.keys()) or ".",
                    }
                )
    if records and missing_tags:
        print(
            f"INFO: {what}: {path}: no {', '.join(missing_tags)} in the header (e.g. an external "
            "truth set); type/length from REF/ALT, VAF strata 'unknown'",
            file=sys.stderr,
        )
    return records


def load_regions(path):
    """{chrom: [(start, end)]} of a BED."""
    regions = {}
    with open(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                sys.exit(f"ERROR: {path}: line {line_no}: expected at least 3 columns")
            regions.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    return regions


def inside(chrom, pos, ref, region_sets):
    """True if the REF span overlaps an interval of every region set."""
    start, end = pos - 1, pos - 1 + len(ref)
    return all(
        any(s < end and start < e for s, e in regions.get(chrom, ())) for regions in region_sets
    )


def read_isec_sites(folder, region_sets, by_position=False):
    """(truth only, query only, both) counts from a bcftools isec sites.txt.

    by_position=False: one count per sites.txt line (exact allele match, isec -c none).
    by_position=True:  one count per (CHROM, POS). `bcftools isec -c all` still writes
    a separate line for an extra allele at a shared position (e.g. a split
    multiallelic call A>C,G against the truth A>C gives 'A C 11' and 'A G 01'); the
    flags of all lines at a position are combined, so the position counts as found
    in both VCFs if any truth and any query record is there.
    """
    path = Path(folder) / "sites.txt"
    if not path.exists():
        sys.exit(
            f"ERROR: {path} does not exist; --isec-* must be a `bcftools isec -p` output folder"
        )
    # key (line number or position) -> [in truth, in query]
    presence = {}
    with open(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            fields = line.rstrip("\n").split("\t")
            # CHROM, POS, REF, ALT, flags (one character per input VCF: truth, query)
            if len(fields) < 5 or fields[4] not in ("10", "01", "11"):
                sys.exit(
                    f"ERROR: {path}: line {line_no}: expected CHROM POS REF ALT and flags "
                    "10/01/11 (two input VCFs: truth, query)"
                )
            # only sites inside all comparison regions (vcfeval was restricted to them)
            if region_sets and not inside(fields[0], int(fields[1]), fields[2], region_sets):
                continue
            key = (fields[0], fields[1]) if by_position else line_no
            flags = presence.setdefault(key, [False, False])
            # combine the presence flags of all lines with this key
            flags[0] |= fields[4][0] == "1"
            flags[1] |= fields[4][1] == "1"
    truth_only = sum(1 for t, q in presence.values() if t and not q)
    query_only = sum(1 for t, q in presence.values() if q and not t)
    both = sum(1 for t, q in presence.values() if t and q)
    return truth_only, query_only, both


def parse_bins(text):
    """VAF bin edges '0,0.05,...' -> [(lo, hi), ...]; strictly increasing."""
    try:
        edges = [float(x) for x in text.split(",")]
    except ValueError:
        sys.exit(f"ERROR: --vaf-bins '{text}' must be comma-separated numbers")
    if len(edges) < 2 or any(b <= a for a, b in zip(edges, edges[1:], strict=False)):
        sys.exit(f"ERROR: --vaf-bins '{text}' must be at least two strictly increasing numbers")
    return list(zip(edges[:-1], edges[1:], strict=False))


def bin_label(lo, hi):
    return f"[{lo:g},{hi:g})"


def find_bin(value, bins):
    """Label of the VAF bin of value; the upper edge of the last bin is included."""
    if value is None:
        return "unknown"
    for i, (lo, hi) in enumerate(bins):
        last = i == len(bins) - 1
        if lo <= value < hi or (last and value == hi):
            return bin_label(lo, hi)
    return "unknown"


def metric_row(keys, stratum_type, stratum, tpb, fn, tpc, fp):
    """One metrics row; recall needs TP_baseline/FN, precision TP_call/FP (None = not used)."""
    recall = ratio(tpb, tpb + fn) if tpb is not None else math.nan
    precision = ratio(tpc, tpc + fp) if tpc is not None else math.nan
    row = dict(
        keys,
        stratum_type=stratum_type,
        stratum=stratum,
        tp_baseline=tpb,
        fn=fn,
        tp_call=tpc,
        fp=fp,
        recall=recall,
        precision=precision,
        f1=f1(recall, precision),
    )
    return {k: fmt(v) if k not in keys else v for k, v in row.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--tp-baseline", required=True, help="vcfeval tp-baseline.vcf.gz")
    parser.add_argument("--fn", required=True, help="vcfeval fn.vcf.gz")
    parser.add_argument("--tp", required=True, help="vcfeval tp.vcf.gz")
    parser.add_argument("--fp", required=True, help="vcfeval fp.vcf.gz")
    parser.add_argument("--isec-exact", required=True, help="bcftools isec -c none output folder")
    parser.add_argument("--isec-position", required=True, help="bcftools isec -c all output folder")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--truth", required=True, help="truth technology (or known-variant list)")
    parser.add_argument("--query", required=True, help="query technology")
    parser.add_argument("--vaf-bins", default="0,0.05,0.1,0.2,0.5,1")
    parser.add_argument("--regions", nargs="*", default=[], help="restrict match strata (BEDs)")
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    keys = {"sample": args.sample, "truth": args.truth, "query": args.query}
    name = f"{args.sample} {args.truth} vs {args.query}"
    bins = parse_bins(args.vaf_bins)

    # the four vcfeval outputs
    tpb = read_records(args.tp_baseline, f"{name} (tp-baseline)")
    fn = read_records(args.fn, f"{name} (fn)")
    tpc = read_records(args.tp, f"{name} (tp)")
    fp = read_records(args.fp, f"{name} (fp)")
    print(
        f"INFO: {name}: vcfeval TP_baseline {len(tpb)}, FN {len(fn)}, TP_call {len(tpc)}, "
        f"FP {len(fp)}",
        file=sys.stderr,
    )

    # overall
    rows = [metric_row(keys, "overall", "all", len(tpb), len(fn), len(tpc), len(fp))]
    # by variant type: types not in TYPE_GROUP would be a harmonization bug
    for records in (tpb, fn, tpc, fp):
        unknown = {r["type"] for r in records} - set(TYPE_GROUP)
        if unknown:
            sys.exit(f"ERROR: {name}: unknown variant type(s) {', '.join(sorted(unknown))}")
    for group in ("SNV", "INDEL"):

        def count(records, group=group):
            return sum(TYPE_GROUP[r["type"]] == group for r in records)

        rows.append(
            metric_row(keys, "variant_type", group, count(tpb), count(fn), count(tpc), count(fp))
        )
    # by VAF bin: truth VAF for recall, query VAF for precision; 'unknown' only if present
    labels = [bin_label(lo, hi) for lo, hi in bins] + ["unknown"]
    for label in labels:

        def in_bin(records, label=label):
            return sum(find_bin(r["vaf"], bins) == label for r in records)

        if label != "unknown" or in_bin(tpb) + in_bin(fn):
            rows.append(metric_row(keys, "truth_vaf", label, in_bin(tpb), in_bin(fn), None, None))
        if label != "unknown" or in_bin(tpc) + in_bin(fp):
            rows.append(metric_row(keys, "query_vaf", label, None, None, in_bin(tpc), in_bin(fp)))
    n_unknown_vaf = sum(r["vaf"] is None for r in tpc + fp)
    if n_unknown_vaf:
        print(
            f"WARNING: {name}: {n_unknown_vaf} query calls without VC_T_VAF (VAF stratum "
            "'unknown')",
            file=sys.stderr,
        )

    # exact and position-only matching (bcftools isec), optionally inside the regions
    region_sets = [load_regions(p) for p in args.regions]
    for stratum, folder in (("exact", args.isec_exact), ("position", args.isec_position)):
        truth_only, query_only, both = read_isec_sites(folder, region_sets, stratum == "position")
        rows.append(metric_row(keys, "match", stratum, both, truth_only, both, query_only))

    with open(f"{args.prefix}.metrics.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=METRIC_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    # per-variant table; a TP is reported in the truth representation, with the query
    # values of the query call at the same position (same ALT if possible)
    calls_by_pos = {}
    for r in tpc:
        calls_by_pos.setdefault((r["chrom"], r["pos"]), []).append(r)
    out = []
    n_tp_unpaired = 0
    for status, records in (("TP", tpb), ("FN", fn)):
        for r in records:
            call = None
            if status == "TP":
                candidates = calls_by_pos.get((r["chrom"], r["pos"]), [])
                call = next((c for c in candidates if c["alt"] == r["alt"]), None) or (
                    candidates[0] if candidates else None
                )
                # different representations (e.g. MNV vs SNVs) may have no call at this position
                n_tp_unpaired += call is None
            out.append(
                dict(
                    keys,
                    status=status,
                    chrom=r["chrom"],
                    pos=r["pos"],
                    ref=r["ref"],
                    alt=r["alt"],
                    type=r["type"],
                    length=r["length"],
                    truth_vaf=fmt(r["vaf"]),
                    truth_dp=fmt(r["dp"]),
                    truth_alt=fmt(r["alt_count"]),
                    query_vaf=fmt(call["vaf"]) if call else "NA",
                    query_dp=fmt(call["dp"]) if call else "NA",
                    query_alt=fmt(call["alt_count"]) if call else "NA",
                    query_filter=call["filter"] if call else "NA",
                )
            )
    for r in fp:
        out.append(
            dict(
                keys,
                status="FP",
                chrom=r["chrom"],
                pos=r["pos"],
                ref=r["ref"],
                alt=r["alt"],
                type=r["type"],
                length=r["length"],
                truth_vaf="NA",
                truth_dp="NA",
                truth_alt="NA",
                query_vaf=fmt(r["vaf"]),
                query_dp=fmt(r["dp"]),
                query_alt=fmt(r["alt_count"]),
                query_filter=r["filter"],
            )
        )
    if n_tp_unpaired:
        print(
            f"INFO: {name}: {n_tp_unpaired} TP in the truth representation without a query call "
            "at the same position (other representation); query values NA",
            file=sys.stderr,
        )
    out.sort(key=lambda r: (r["chrom"], r["pos"], r["ref"], r["alt"], r["status"]))
    with open(f"{args.prefix}.variants.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=VARIANT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out)


if __name__ == "__main__":
    main()
