"""Unit tests for bin/benchmark_metrics.py and bin/aggregate_metrics.py."""

import csv
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[2] / "bin"
HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##INFO=<ID=VC_T_VAF,Number=1,Type=Float,Description="vaf">
##INFO=<ID=VC_T_DP,Number=1,Type=Integer,Description="dp">
##INFO=<ID=VC_T_ALT,Number=1,Type=Integer,Description="alt">
##INFO=<ID=VC_TYPE,Number=1,Type=String,Description="type">
##INFO=<ID=VC_LEN,Number=1,Type=Integer,Description="len">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def vcf(path, records):
    lines = []
    for pos, ref, alt, vaf in records:
        vtype = "SNV" if len(ref) == len(alt) else ("INS" if len(alt) > len(ref) else "DEL")
        info = (
            f"VC_T_VAF={vaf};VC_T_DP=100;VC_T_ALT={int(vaf * 100)};"
            f"VC_TYPE={vtype};VC_LEN={len(alt) - len(ref)}"
        )
        lines.append(f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\n")
    path.write_text(HEADER + "".join(lines))
    return path


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture()
def metrics(tmp_path):
    truth = [(100, "A", "C", 0.3), (200, "A", "G", 0.04), (300, "AT", "A", 0.3)]
    query = [(100, "A", "C", 0.28), (200, "A", "T", 0.05), (400, "G", "T", 0.6)]
    files = {
        "tp-baseline": vcf(tmp_path / "tpb.vcf", [truth[0]]),
        "fn": vcf(tmp_path / "fn.vcf", truth[1:]),
        "tp": vcf(tmp_path / "tp.vcf", [query[0]]),
        "fp": vcf(tmp_path / "fp.vcf", query[1:]),
        "truth-vcf": vcf(tmp_path / "truth.vcf", truth),
        "query-vcf": vcf(tmp_path / "query.vcf", query),
    }
    cmd = [sys.executable, str(BIN / "benchmark_metrics.py"), "--sample", "S1", "--truth", "wgs",
           "--query", "wes", "--prefix", str(tmp_path / "S1_wes")]  # fmt: skip
    for key, path in files.items():
        cmd += [f"--{key}", str(path)]
    subprocess.run(cmd, check=True)
    return tmp_path / "S1_wes"


def row(rows, stratum_type, stratum):
    return next(r for r in rows if r["stratum_type"] == stratum_type and r["stratum"] == stratum)


def test_overall_and_types(metrics):
    rows = read_tsv(f"{metrics}.metrics.tsv")
    overall = row(rows, "overall", "all")
    assert (overall["tp_baseline"], overall["fn"], overall["tp_call"], overall["fp"]) == (
        "1", "2", "1", "2"
    )  # fmt: skip
    assert overall["recall"] == "0.3333"
    assert overall["precision"] == "0.3333"
    indel = row(rows, "variant_type", "INDEL")
    assert (indel["tp_baseline"], indel["fn"], indel["recall"]) == ("0", "1", "0.0000")


def test_vaf_bins(metrics):
    rows = read_tsv(f"{metrics}.metrics.tsv")
    low = row(rows, "truth_vaf", "[0,0.05)")
    assert (low["tp_baseline"], low["fn"], low["recall"], low["precision"]) == (
        "0",
        "1",
        "0.0000",
        "NA",
    )
    high = row(rows, "query_vaf", "[0.5,1)")
    assert (high["tp_call"], high["fp"], high["precision"]) == ("0", "1", "0.0000")


def test_match_strata(metrics):
    rows = read_tsv(f"{metrics}.metrics.tsv")
    exact = row(rows, "match", "exact")
    position = row(rows, "match", "position")
    assert (exact["tp_baseline"], exact["fn"], exact["fp"]) == ("1", "2", "2")
    # chr1:200 has different alleles: matched by position only
    assert (position["tp_baseline"], position["fn"], position["fp"]) == ("2", "1", "1")


def test_variant_table(metrics):
    rows = read_tsv(f"{metrics}.variants.tsv")
    assert [(r["pos"], r["status"]) for r in rows] == [
        ("100", "TP"), ("200", "FN"), ("200", "FP"), ("300", "FN"), ("400", "FP")
    ]  # fmt: skip
    tp = rows[0]
    assert (tp["truth_vaf"], tp["query_vaf"]) == ("0.3000", "0.2800")


def test_aggregate(metrics, tmp_path):
    second = tmp_path / "S2_wes.metrics.tsv"
    text = Path(f"{metrics}.metrics.tsv").read_text().replace("S1\t", "S2\t")
    second.write_text(text)
    prefix = tmp_path / "agg"
    subprocess.run([sys.executable, str(BIN / "aggregate_metrics.py"), f"{metrics}.metrics.tsv",
                    str(second), "--prefix", str(prefix)], check=True)  # fmt: skip
    agg = read_tsv(f"{prefix}.aggregated.tsv")
    overall = row(agg, "overall", "all")
    assert overall["n_samples"] == "2"
    assert (overall["tp_baseline"], overall["fn"]) == ("2", "4")
    assert overall["micro_recall"] == "0.3333"
    assert overall["macro_recall_mean"] == "0.3333"
    assert overall["macro_recall_sd"] == "0.0000"
    assert len(read_tsv(f"{prefix}.per_pair.tsv")) == 2 * len(read_tsv(f"{metrics}.metrics.tsv"))
