#!/usr/bin/env python3
"""Aggregate per-pair metrics over samples.

Input: any number of metrics tables with the columns of benchmark_metrics.py
(<prefix>.metrics.tsv) or classify_discordance.py (<prefix>.power.tsv):
sample, truth, query, stratum_type, stratum, tp_baseline, fn, tp_call, fp, recall,
precision, f1. Counts and metrics are integers / numbers or NA.

Outputs (--prefix):
  <prefix>.per_pair.tsv     all input rows, sorted (truth, query, stratum type, stratum, sample)
  <prefix>.aggregated.tsv   one row per truth/query technology and stratum:
      n_samples   number of samples with a row
      tp_baseline, fn, tp_call, fp
                  counts summed over the samples (NA if no sample has a count)
      micro_*     recall/precision/F1 from the summed counts
      macro_*     mean and standard deviation (ddof 1) of the per-sample values;
                  samples with an undefined (NA) value are skipped, sd needs >= 2

Checks: required columns, one row per sample and stratum (duplicates would be
counted twice), counts are non-negative integers or NA, metrics are within [0, 1]
or NA.

Why a script: summing counts and averaging per-sample metrics over groups is a
pandas groupby; no bioinformatics tool aggregates these tables.
"""

import argparse
import math
import sys

import pandas as pd

KEYS = ["truth", "query", "stratum_type", "stratum"]
COUNTS = ["tp_baseline", "fn", "tp_call", "fp"]
METRICS = ["recall", "precision", "f1"]
REQUIRED = ["sample", *KEYS, *COUNTS, *METRICS]
# order of the stratum types in the outputs; unknown types go last
ORDER = {"overall": 0, "variant_type": 1, "truth_vaf": 2, "query_vaf": 3, "match": 4,
         "power_adjusted": 5}  # fmt: skip


def safe_ratio(num, den):
    """num / den, NaN if den is 0 or NaN."""
    return num / den if den and not math.isnan(den) else math.nan


def read_metrics(path):
    """One metrics table as strings, with checked columns and values."""
    data = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in REQUIRED if c not in data.columns]
    if missing:
        sys.exit(f"ERROR: {path}: missing column(s) {', '.join(missing)}; expected a metrics table "
                 "of benchmark_metrics.py or classify_discordance.py")  # fmt: skip
    for column in COUNTS:
        # counts: non-negative integers or NA
        bad = data[~data[column].str.fullmatch(r"\d+|NA")]
        if len(bad):
            line = bad.index[0] + 2
            sys.exit(f"ERROR: {path}: line {line}: {column} '{bad[column].iloc[0]}' is not a "
                     "non-negative integer or NA")  # fmt: skip
    for column in METRICS:
        # metrics: numbers within [0, 1] or NA
        values = pd.to_numeric(data[column].replace("NA", None), errors="coerce")
        bad = data[(values.isna() & (data[column] != "NA")) | (values < 0) | (values > 1)]
        if len(bad):
            line = bad.index[0] + 2
            sys.exit(f"ERROR: {path}: line {line}: {column} '{bad[column].iloc[0]}' is not a "
                     "number within [0, 1] or NA")  # fmt: skip
    return data


def micro(group):
    """Summed counts and the recall/precision/F1 computed from them."""
    # NA counts are skipped; min_count=1 keeps NA if no sample has a count
    sums = {c: pd.to_numeric(group[c].replace("NA", None)).sum(min_count=1) for c in COUNTS}
    recall = safe_ratio(sums["tp_baseline"], sums["tp_baseline"] + sums["fn"])
    precision = safe_ratio(sums["tp_call"], sums["tp_call"] + sums["fp"])
    # F1 is undefined if one of recall/precision is undefined or both are 0
    if math.isnan(recall) or math.isnan(precision) or recall + precision == 0:
        f1 = math.nan
    else:
        f1 = 2 * recall * precision / (recall + precision)
    return sums, {"recall": recall, "precision": precision, "f1": f1}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("metrics", nargs="+", help="per-pair metrics TSV files")
    parser.add_argument("--prefix", required=True, help="output prefix")
    args = parser.parse_args()

    # all input rows in one table
    data = pd.concat([read_metrics(p) for p in args.metrics], ignore_index=True)
    if data.empty:
        sys.exit(f"ERROR: no metrics rows in {len(args.metrics)} input file(s)")
    # one row per sample and stratum: a duplicate would be counted twice
    duplicated = data[data.duplicated(["sample", *KEYS], keep=False)]
    if len(duplicated):
        first = duplicated.iloc[0]
        sys.exit(f"ERROR: {len(duplicated)} rows for the same sample and stratum, e.g. "
                 f"{first['sample']} {first['truth']} vs {first['query']} {first['stratum_type']} "
                 f"{first['stratum']}; was a metrics file given twice?")  # fmt: skip

    # per-pair table, sorted by the stratum type order above
    data["_order"] = data["stratum_type"].map(ORDER).fillna(len(ORDER))
    data = data.sort_values(["truth", "query", "_order", "stratum", "sample"], kind="stable")
    data.drop(columns="_order").to_csv(f"{args.prefix}.per_pair.tsv", sep="\t", index=False)

    rows = []
    for key, group in data.groupby(KEYS, sort=False):
        sums, micro_values = micro(group)
        row = dict(zip(KEYS, key, strict=True))
        row["n_samples"] = group["sample"].nunique()
        for c in COUNTS:
            row[c] = "NA" if pd.isna(sums[c]) else int(sums[c])
        for m in METRICS:
            row[f"micro_{m}"] = micro_values[m]
            # per-sample values without the undefined ones
            values = pd.to_numeric(group[m].replace("NA", None)).dropna()
            row[f"macro_{m}_mean"] = values.mean() if len(values) else math.nan
            row[f"macro_{m}_sd"] = values.std(ddof=1) if len(values) > 1 else math.nan
        rows.append(row)
    pd.DataFrame(rows).to_csv(f"{args.prefix}.aggregated.tsv", sep="\t", index=False, na_rep="NA",
                              float_format="%.4f")  # fmt: skip
    print(f"INFO: {len(args.metrics)} files, {len(data)} rows, {data['sample'].nunique()} samples "
          f"-> {len(rows)} aggregated rows", file=sys.stderr)  # fmt: skip


if __name__ == "__main__":
    main()
