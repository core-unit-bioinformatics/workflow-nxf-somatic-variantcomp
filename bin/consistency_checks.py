#!/usr/bin/env python3
"""Consistency checks of a variantcomp run: do the counts of the steps add up?

Subcommands
  check   evaluate all checks and write one table (never stops: the report shows it)
  gate    stop with exit status 1 if the table has a failed check (run after the report)

Checks of `check` (scope = the sample/technology or comparison the check is about)
  accounting   the checks of variant_accounting.py, one table per VCF: input = removed +
               kept per step, bcftools norm and harmonization counts, compared VCF = records
               passing the filters and region
  benchmark    per pair: TP_baseline + FN >= compared truth records and TP_call + FP >=
               compared query records (equal unless rtg vcfeval --decompose splits complex
               records); FN and FP rows of the variant table = FN and FP of the metrics
  isec         per pair, bcftools isec exact: the isec status of the truth records in the
               variant table (TP, FN) adds up to the TP_baseline and FN of the isec metrics
               (skipped if vcfeval decomposed records isec can't know; position mode counts
               positions, not records, so it has no such check)
  aggregation  per comparison: the summed (micro) counts of metrics.aggregated.tsv = the sum
               of the per-pair counts
  discordance  per comparison and status: the main reasons (classes) add up to FN resp.
               FP; every FN/FP row has a known main reason
  reasons      per comparison, recomputed from the columns of the variant table (number
               of rows breaking a rule must be 0): exactly one main reason per FN/FP, the
               first of `reasons` and equal to `class`; called reasons (filtered_non_pass,
               present_but_filtered) if and only if the other technology has the call;
               absent if and only if there is no tumor read evidence; low_coverage,
               subthreshold_evidence, no_evidence agree with depth (--min-coverage) and alt
               reads; germline_in_normal agrees with the normal reads (--min-alt-reads,
               --germline-vaf); indel, low_vaf (--low-vaf), in_repeat and every
               stratification agree with their columns; no unknown reason; TP rows have
               only `concordant`
  evidence     per comparison and variant group (SNV, indel): the share of TP variants with
               alt reads in the query tumor alignment >= --evidence-min-fraction (a TP was
               called by both technologies, so its reads must show it; a low share points
               to a read-counting problem, e.g. indel alleles not found)
  regions      the last region step = the common region of the target summary; every
               sample's callable region <= the common region
  known        per required list, sample and technology: passed + failed = required

INPUT
  --accounting      <x>.checks.tsv and <x>.accounting.tsv of variant_accounting.py (several)
  --per-pair        benchmark/metrics.per_pair.tsv     --aggregated  metrics.aggregated.tsv
  --variants        discordance/variants.annotated.tsv --classes     report classes table
  --region-steps    regions/common_targets.steps.tsv   --targets     common_targets.summary.tsv
  --callable        coverage/callable_summary.tsv (optional)
  --required        known_variants/required_summary.tsv (optional)
OUTPUT
  --output          TSV: category, scope, check, expected, observed, status (ok / FAIL)

Why a script: the checks combine tables of different steps; no tool does this.
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict

# main reasons of classify_discordance.py (exactly one per FN/FP)
CALLED = ["filtered_non_pass", "present_but_filtered"]
NOT_CALLED = ["absent", "low_coverage", "subthreshold_evidence", "no_evidence"]
CLASSES = CALLED + NOT_CALLED
# further reasons with a fixed meaning; every other reason must be a stratification column
FURTHER = ["germline_in_normal", "indel", "low_vaf", "in_repeat"]
INDEL_TYPES = ("INS", "DEL", "COMPLEX")
# yes/no columns of the variant table that are not stratifications
NOT_STRATA = {"other_in_normalized", "detectable", "in_repeat"}


def read_tsv(path):
    """Rows of a TSV with header as dicts (empty list if path is None)."""
    if not path:
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class Checks:
    """The collected checks."""

    def __init__(self):
        self.rows = []

    def equal(self, category, scope, name, expected, observed):
        """A check that two values are equal."""
        status = "ok" if str(expected) == str(observed) else "FAIL"
        self.rows.append((category, scope, name, str(expected), str(observed), status))

    def at_least(self, category, scope, name, minimum, observed):
        """A check that a value is at least a minimum."""
        status = "ok" if int(observed) >= int(minimum) else "FAIL"
        self.rows.append((category, scope, name, f">= {minimum}", str(observed), status))


def accounting_checks(paths, checks):
    """Copy the checks of variant_accounting.py; return compared records per sample/technology."""
    compared = {}
    for path in paths:
        rows = read_tsv(path)
        if path.endswith(".checks.tsv"):
            for r in rows:
                scope = f"{r['sample']} {r['technology']}"
                checks.rows.append(
                    ("accounting", scope, r["check"], r["expected"], r["observed"], r["status"])
                )
        elif path.endswith(".accounting.tsv"):
            for r in rows:
                # step 7, outcome kept: the records of the compared VCF
                if r["step"] == "7" and r["outcome"] == "kept":
                    compared[(r["sample"], r["technology"])] = int(r["records"])
        else:
            sys.exit(
                f"ERROR: {path}: expected a .checks.tsv or .accounting.tsv of variant_accounting.py"
            )
    return compared


def benchmark_checks(per_pair, variants, compared, checks):
    """Per pair: metrics vs. compared records and vs. the variant table."""
    status_rows = Counter((v["sample"], v["truth"], v["query"], v["status"]) for v in variants)
    for r in per_pair:
        if (r["stratum_type"], r["stratum"]) != ("overall", "all"):
            continue
        key = (r["sample"], r["truth"], r["query"])
        scope = f"{r['sample']} {r['truth']} vs {r['query']}"
        truth = compared.get((r["sample"], r["truth"]))
        query = compared.get((r["sample"], r["query"]))
        if truth is not None:
            checks.at_least(
                "benchmark",
                scope,
                "TP_baseline + FN vs. compared truth records",
                truth,
                int(r["tp_baseline"]) + int(r["fn"]),
            )
        if query is not None:
            checks.at_least(
                "benchmark",
                scope,
                "TP_call + FP vs. compared query records",
                query,
                int(r["tp_call"]) + int(r["fp"]),
            )
        checks.equal(
            "benchmark",
            scope,
            "FN rows of the variant table = FN of the metrics",
            r["fn"],
            status_rows[(*key, "FN")],
        )
        checks.equal(
            "benchmark",
            scope,
            "FP rows of the variant table = FP of the metrics",
            r["fp"],
            status_rows[(*key, "FP")],
        )


def isec_checks(per_pair, variants, checks):
    """Per pair and isec method: the isec status of the truth records adds up to the isec
    metrics (verifies that the isec counts reach the metrics unchanged)."""
    if not variants or "isec_exact" not in variants[0]:
        return
    metrics = {
        (r["sample"], r["truth"], r["query"], r["stratum"]): r
        for r in per_pair
        if r["stratum_type"] == "match"
    }
    for method, column in (("exact", "isec_exact"), ("position", "isec_position")):
        counts = defaultdict(Counter)
        for v in variants:
            # truth records: the TP and FN rows of vcfeval
            if v["status"] in ("TP", "FN"):
                counts[(v["sample"], v["truth"], v["query"])][v[column]] += 1
        for pair, c in sorted(counts.items()):
            m = metrics.get((*pair, method))
            if m is None or c["NA"]:
                continue
            scope = f"{pair[0]} {pair[1]} vs {pair[2]}"
            # position mode counts positions, the variant table records: several truth
            # records at one position count once in the metrics
            if method == "position":
                continue
            checks.equal(
                "isec",
                scope,
                f"isec {method}: truth records found by the query (variant table) = TP_baseline",
                m["tp_baseline"],
                c["TP"],
            )
            checks.equal(
                "isec",
                scope,
                f"isec {method}: truth records missed by the query (variant table) = FN",
                m["fn"],
                c["FN"],
            )


def aggregation_checks(per_pair, aggregated, checks):
    """Per comparison and stratum: micro counts = sum of the per-pair counts."""
    sums = defaultdict(Counter)
    for r in per_pair:
        # only the overall rows are compared (strata like truth_vaf have NA for calls)
        if r["stratum_type"] != "overall":
            continue
        key = (r["truth"], r["query"], r["stratum_type"], r["stratum"])
        for col in ("tp_baseline", "fn", "tp_call", "fp"):
            sums[key][col] += int(r[col])
    for r in aggregated:
        key = (r["truth"], r["query"], r["stratum_type"], r["stratum"])
        if r["stratum_type"] != "overall":
            continue
        for col in ("tp_baseline", "fn", "tp_call", "fp"):
            checks.equal(
                "aggregation",
                f"{r['truth']} vs {r['query']}",
                f"{col}: aggregated = sum over samples",
                r[col],
                sums[key][col],
            )


def discordance_checks(per_pair, variants, classes, checks):
    """Per comparison and status: classes add up to FN / FP; every FN/FP has a known class."""
    totals = Counter()
    for r in per_pair:
        if (r["stratum_type"], r["stratum"]) == ("overall", "all"):
            totals[(r["truth"], r["query"], "FN")] += int(r["fn"])
            totals[(r["truth"], r["query"], "FP")] += int(r["fp"])
    for r in classes:
        observed = sum(int(r.get(c) or 0) for c in CLASSES)
        checks.equal(
            "discordance",
            f"{r['truth']} vs {r['query']} {r['status']}",
            "discordance classes add up to the variants",
            totals[(r["truth"], r["query"], r["status"])],
            observed,
        )
    unknown = [v for v in variants if v["status"] in ("FN", "FP") and v.get("class") not in CLASSES]
    checks.equal(
        "discordance", "all comparisons", "FN/FP rows without a known class", 0, len(unknown)
    )


def number(text):
    """A table value as float, None for NA / empty / '.'."""
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def reason_errors(row, strata, args):
    """The rules of the reason model that one FN/FP row breaks (list of short texts)."""
    errors = []
    reasons = row["reasons"].split(",")
    main = [r for r in reasons if r in CLASSES]
    # exactly one main reason, first in the list and equal to class
    if len(main) != 1 or reasons[0] != main[0] or row["class"] != main[0]:
        return ["main reason"]
    main = main[0]
    # called reasons if and only if the other technology has the call
    if (main in CALLED) != (row.get("other_in_normalized") == "yes"):
        errors.append("called")
    has_evidence = "other_tumor_dp" in row
    dp, alt = number(row.get("other_tumor_dp")), number(row.get("other_tumor_alt"))
    if main in NOT_CALLED and has_evidence:
        # absent if and only if there is no tumor read evidence
        if (main == "absent") != (dp is None):
            errors.append("absent")
        elif dp is not None:
            expected = (
                "low_coverage"
                if dp < args.min_coverage
                else ("subthreshold_evidence" if alt else "no_evidence")
            )
            if main != expected:
                errors.append("evidence")
    elif main in NOT_CALLED and main != "absent":
        # without read evidence columns only 'absent' is possible
        errors.append("evidence")
    # germline_in_normal agrees with the normal reads
    n_alt, n_vaf = number(row.get("normal_alt")), number(row.get("normal_vaf"))
    germline = (n_alt or 0) >= args.min_alt_reads and (n_vaf or 0) >= args.germline_vaf
    if germline != ("germline_in_normal" in reasons):
        errors.append("germline")
    # context reasons agree with their columns (only if the column exists)
    if "type" in row and (row["type"] in INDEL_TYPES) != ("indel" in reasons):
        errors.append("indel")
    vaf = number(row.get("truth_vaf" if row["status"] == "FN" else "query_vaf"))
    if (vaf is not None and vaf < args.low_vaf) != ("low_vaf" in reasons):
        errors.append("low_vaf")
    if "in_repeat" in row and (row["in_repeat"] == "yes") != ("in_repeat" in reasons):
        errors.append("in_repeat")
    for name in strata:
        if (row.get(name) == "yes") != (name in reasons):
            errors.append(name)
    # no reason outside the vocabulary
    if [r for r in reasons if r not in CLASSES + FURTHER + strata]:
        errors.append("unknown reason")
    return errors


def reason_checks(variants, checks, args):
    """Per comparison: number of rows that break a rule of the reason model (must be 0)."""
    if not variants or "reasons" not in variants[0]:
        checks.equal(
            "reasons", "all comparisons", "column reasons in the variant table", "yes", "no"
        )
        return
    columns = list(variants[0])
    # stratification columns: yes/no columns other than the fixed ones
    strata = [
        c
        for c in columns
        if c not in NOT_STRATA
        and all(v[c] in ("yes", "no") for v in variants)
        and not c.startswith("found_")
    ]
    bad = defaultdict(Counter)
    totals = Counter()
    for v in variants:
        scope = f"{v['sample']} {v['truth']} vs {v['query']}"
        totals[scope] += 1
        if v["status"] == "TP":
            # concordant variants have no reason
            if v["reasons"] != "concordant" or v["class"] != "concordant":
                bad[scope]["TP not concordant"] += 1
            continue
        for error in reason_errors(v, strata, args):
            bad[scope][error] += 1
    for scope in sorted(totals):
        errors = bad[scope]
        detail = ", ".join(f"{k}: {n}" for k, n in sorted(errors.items())) or "none"
        checks.equal(
            "reasons",
            scope,
            f"rows breaking a rule of the reason model ({detail})",
            0,
            sum(errors.values()),
        )


def evidence_checks(variants, checks, args):
    """Per comparison and variant group: share of TP variants with alt reads in the query
    tumor alignment (only with read evidence)."""
    groups = defaultdict(lambda: [0, 0])
    for v in variants:
        alt = number(v.get("other_tumor_alt"))
        if v["status"] != "TP" or alt is None:
            continue
        group = "indel" if v.get("type") in INDEL_TYPES else "SNV"
        key = (f"{v['sample']} {v['truth']} vs {v['query']}", group)
        groups[key][0] += 1
        groups[key][1] += alt > 0
    for (scope, group), (n, with_alt) in sorted(groups.items()):
        share = with_alt / n
        checks.rows.append(
            (
                "evidence",
                scope,
                f"{group} TP with alt reads in the query tumor alignment ({with_alt} of {n})",
                f">= {args.evidence_min_fraction:g}",
                f"{share:.3f}",
                "ok" if share >= args.evidence_min_fraction else "FAIL",
            )
        )


def region_checks(steps, targets, callable_rows, checks):
    """Region steps end with the common region; callable regions within it."""
    common = next((int(r["bp"]) for r in targets if r["region"] == "common"), None)
    if steps and common is not None:
        checks.equal(
            "regions",
            "common region",
            "last region step = common region (bp)",
            common,
            steps[-1]["bp"],
        )
    for r in callable_rows:
        checks.at_least(
            "regions",
            r["sample"],
            "common region bp >= callable bp",
            r["callable_bp"],
            r["common_bp"],
        )


def known_checks(required, checks):
    """Required lists: passed + failed = required."""
    for r in required:
        scope = f"{r['list']} {r['sample']} {r['technology']}"
        checks.equal(
            "known",
            scope,
            "passed + failed = required",
            r["required"],
            int(r["passed"]) + int(r["failed"]),
        )


def check(args):
    """Evaluate all checks and write the table."""
    checks = Checks()
    compared = accounting_checks(args.accounting, checks)
    per_pair, variants = read_tsv(args.per_pair), read_tsv(args.variants)
    benchmark_checks(per_pair, variants, compared, checks)
    isec_checks(per_pair, variants, checks)
    aggregation_checks(per_pair, read_tsv(args.aggregated), checks)
    discordance_checks(per_pair, variants, read_tsv(args.classes), checks)
    reason_checks(variants, checks, args)
    evidence_checks(variants, checks, args)
    region_checks(
        read_tsv(args.region_steps), read_tsv(args.targets), read_tsv(args.callable), checks
    )
    known_checks(read_tsv(args.required), checks)
    with open(args.output, "w") as fh:
        fh.write("category\tscope\tcheck\texpected\tobserved\tstatus\n")
        fh.writelines("\t".join(row) + "\n" for row in checks.rows)
    failed = [r for r in checks.rows if r[5] == "FAIL"]
    print(f"INFO: {len(checks.rows)} consistency checks, {len(failed)} failed", file=sys.stderr)
    for r in failed:
        print(
            f"WARNING: check failed: [{r[0]}] {r[1]}: {r[2]}: expected {r[3]}, observed {r[4]}",
            file=sys.stderr,
        )


def gate(args):
    """Stop with an error if a check failed (after the report was written)."""
    failed = [r for r in read_tsv(args.checks) if r["status"] != "ok"]
    if failed:
        for r in failed:
            print(
                f"ERROR: consistency check failed: [{r['category']}] {r['scope']}: {r['check']}: "
                f"expected {r['expected']}, observed {r['observed']}",
                file=sys.stderr,
            )
        sys.exit(
            f"ERROR: {len(failed)} consistency checks failed (see report/consistency_checks.tsv "
            "and the report); the results may be wrong"
        )
    print("INFO: all consistency checks passed", file=sys.stderr)


def main():
    """Parse the command line and run the subcommand."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    c = sub.add_parser("check")
    c.add_argument("--accounting", nargs="*", default=[])
    for name in ("per-pair", "aggregated", "variants", "classes", "region-steps", "targets"):
        c.add_argument(f"--{name}", required=True)
    c.add_argument("--callable")
    c.add_argument("--required")
    c.add_argument("--output", required=True)
    # thresholds of the reason model (pipeline params of the same name)
    c.add_argument("--min-coverage", type=float, default=10, help="--min_coverage_tumor")
    c.add_argument("--min-alt-reads", type=float, default=2, help="--evidence_min_alt_reads")
    c.add_argument("--germline-vaf", type=float, default=0.1, help="--germline_vaf")
    c.add_argument("--low-vaf", type=float, default=0.1, help="--low_vaf")
    c.add_argument(
        "--evidence-min-fraction",
        type=float,
        default=0.8,
        help="--evidence_check_min_fraction: share of TP with alt reads in the query tumor",
    )
    g = sub.add_parser("gate")
    g.add_argument("--checks", required=True)
    args = parser.parse_args()
    check(args) if args.command == "check" else gate(args)


if __name__ == "__main__":
    main()
