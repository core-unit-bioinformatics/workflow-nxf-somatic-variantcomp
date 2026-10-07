#!/usr/bin/env python3
"""Aggregate per-pair metrics over samples.

Input: any number of <prefix>.metrics.tsv files from benchmark_metrics.py.
Output:
  <prefix>.per_pair.tsv     all input rows, sorted
  <prefix>.aggregated.tsv   per truth/query technology and stratum:
      micro_*  recall/precision/F1 from the summed counts of all samples
      macro_*  mean and standard deviation of the per-sample values
               (samples with an undefined value are skipped)
      n_samples number of samples
"""

import argparse
import math

import pandas as pd

COUNTS = ["tp_baseline", "fn", "tp_call", "fp"]
METRICS = ["recall", "precision", "f1"]


def safe_ratio(num, den):
    return num / den if den else math.nan


def micro(group):
    sums = {c: pd.to_numeric(group[c], errors="coerce").sum(min_count=1) for c in COUNTS}
    recall = safe_ratio(sums["tp_baseline"], sums["tp_baseline"] + sums["fn"])
    precision = safe_ratio(sums["tp_call"], sums["tp_call"] + sums["fp"])
    if math.isnan(recall) or math.isnan(precision) or recall + precision == 0:
        f1 = math.nan
    else:
        f1 = 2 * recall * precision / (recall + precision)
    return sums, {"recall": recall, "precision": precision, "f1": f1}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("metrics", nargs="+", help="per-pair metrics TSV files")
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    frames = [pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False) for p in args.metrics]
    data = pd.concat(frames, ignore_index=True)
    order = {"overall": 0, "variant_type": 1, "truth_vaf": 2, "query_vaf": 3, "match": 4}
    data["_order"] = data["stratum_type"].map(order).fillna(9)
    data = data.sort_values(["truth", "query", "_order", "stratum", "sample"], kind="stable")
    data.drop(columns="_order").to_csv(f"{args.prefix}.per_pair.tsv", sep="\t", index=False)

    rows = []
    keys = ["truth", "query", "stratum_type", "stratum"]
    for key, group in data.groupby(keys, sort=False):
        sums, micro_values = micro(group)
        row = dict(zip(keys, key, strict=True))
        row["n_samples"] = group["sample"].nunique()
        for c in COUNTS:
            row[c] = "NA" if pd.isna(sums[c]) else int(sums[c])
        for m in METRICS:
            row[f"micro_{m}"] = micro_values[m]
            values = pd.to_numeric(group[m].replace("NA", None), errors="coerce").dropna()
            row[f"macro_{m}_mean"] = values.mean() if len(values) else math.nan
            row[f"macro_{m}_sd"] = values.std(ddof=1) if len(values) > 1 else math.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(f"{args.prefix}.aggregated.tsv", sep="\t", index=False, na_rep="NA",
               float_format="%.4f")  # fmt: skip


if __name__ == "__main__":
    main()
