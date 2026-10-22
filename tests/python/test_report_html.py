"""Unit tests for bin/report_html.py."""

import json
import os
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "report_html.py"
sys.path.insert(0, str(SCRIPT.parent))

import report_html  # noqa: E402


def tsv(path, header, rows):
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))
    return str(path)


METRIC_COLUMNS = ["tp_baseline", "fn", "tp_call", "fp"]


def inputs(tmp_path, check_status="ok"):
    """A small run: one sample (S1), truth wgs, query wes; one FN and one FP."""
    per_pair = tsv(
        tmp_path / "pp.tsv",
        [
            "sample",
            "truth",
            "query",
            "stratum_type",
            "stratum",
            *METRIC_COLUMNS,
            "recall",
            "precision",
            "f1",
        ],
        [["S1", "wgs", "wes", "overall", "all", 2, 1, 2, 1, 0.667, 0.667, 0.667]],
    )
    metrics = tsv(
        tmp_path / "agg.tsv",
        [
            "truth",
            "query",
            "stratum_type",
            "stratum",
            "n_samples",
            *METRIC_COLUMNS,
            "micro_recall",
            "micro_precision",
            "micro_f1",
        ],
        [
            ["wgs", "wes", "overall", "all", 1, 2, 1, 2, 1, 0.667, 0.667, 0.667],
            ["wgs", "wes", "variant_type", "SNV", 1, 2, 0, 2, 1, 1.0, 0.667, 0.8],
            ["wgs", "wes", "truth_vaf", "(0.2,0.5]", 1, 2, 1, "NA", "NA", 0.667, "NA", "NA"],
        ],
    )
    variant_columns = [
        "sample",
        "truth",
        "query",
        "status",
        "class",
        "chrom",
        "pos",
        "ref",
        "alt",
        "type",
    ]
    variants = tsv(
        tmp_path / "v.tsv",
        variant_columns,
        [
            ["S1", "wgs", "wes", "TP", "concordant", "chr1", 100, "A", "C", "SNV"],
            ["S1", "wgs", "wes", "TP", "concordant", "chr1", 200, "G", "T", "SNV"],
            ["S1", "wgs", "wes", "FN", "low_coverage", "chr1", 300, "AT", "A", "DEL"],
            ["S1", "wgs", "wes", "FP", "no_evidence", "chr1", 400, "C", "G", "SNV"],
        ],
    )
    truth = tsv(
        tmp_path / "truth.tsv",
        [
            "sample",
            "truth",
            "chrom",
            "pos",
            "ref",
            "alt",
            "type",
            "found_by",
            "queries",
            "found_wes",
            "class_wes",
        ],
        [
            ["S1", "wgs", "chr1", 100, "A", "C", "SNV", 1, 1, "yes", "concordant"],
            ["S1", "wgs", "chr1", 200, "G", "T", "SNV", 1, 1, "yes", "concordant"],
            ["S1", "wgs", "chr1", 300, "AT", "A", "DEL", 0, 1, "no", "low_coverage"],
        ],
    )
    classes = tsv(
        tmp_path / "classes.tsv",
        ["truth", "query", "status", *report_html.CLASS_ORDER],
        [
            ["wgs", "wes", "FN", 0, 0, 0, 1, 0, 0, 0],
            ["wgs", "wes", "FP", 0, 0, 0, 0, 0, 1, 0],
        ],
    )
    accounting = tsv(
        tmp_path / "S1_wes.accounting.tsv",
        [
            "sample",
            "technology",
            "step",
            "step_name",
            "outcome",
            "reason",
            "records",
            "snv",
            "indel",
            "other",
        ],
        [
            ["S1", "wes", 1, "input", "input", ".", 5, 4, 1, 0],
            ["S1", "wes", 2, "pass_split", "removed", "filter_not_pass", 1, 1, 0, 0],
            ["S1", "wes", 7, "filters", "removed", "outside_region", 1, 1, 0, 0],
            ["S1", "wes", 7, "filters", "kept", ".", 3, 2, 1, 0],
        ],
    )
    removed = tsv(
        tmp_path / "S1_wes.removed.tsv",
        ["sample", "technology", "step_name", "reason", "chrom", "pos", "ref", "alt", "filter"],
        [["S1", "wes", "pass_split", "filter_not_pass", "chr1", 500, "A", "G", "weak_evidence"]],
    )
    steps = tsv(
        tmp_path / "steps.tsv",
        [
            "step",
            "name",
            "intervals",
            "bp",
            "fraction_of_genome",
            "fraction_of_previous",
            "removed_bp",
        ],
        [
            [1, "1_genome", 1, 10000, 1, 1, "NA"],
            [2, "2_targets_intersection", 2, 2000, 0.2, 0.2, 8000],
        ],
    )
    checks = tsv(
        tmp_path / "checks.tsv",
        ["category", "scope", "check", "expected", "observed", "status"],
        [["benchmark", "S1 wgs vs wes", "FN rows = FN", 1, 1, check_status]],
    )
    sheet = tmp_path / "sheet.csv"
    sheet.write_text("sample,technology,library_type\nS1,wgs,wgs\nS1,wes,capture\n")
    versions = tmp_path / "versions.yml"
    versions.write_text("RTGTOOLS_VCFEVAL:\n  rtg-tools: 3.12.1\nWorkflow:\n  Nextflow: 25.10.4\n")
    return [
        "--per-pair", per_pair, "--metrics", metrics, "--variants", variants,
        "--truth-variants", truth, "--classes", classes, "--accounting", accounting,
        "--removed", removed, "--region-steps", steps, "--checks", checks,
        "--samplesheet", str(sheet), "--versions", str(versions),
    ]  # fmt: skip


def run(tmp_path, *args):
    env = dict(os.environ, MPLCONFIGDIR=str(tmp_path / ".mpl"))
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args], capture_output=True, text=True, env=env, cwd=tmp_path
    )


def test_html_has_sections_plots_and_banner(tmp_path):
    result = run(tmp_path, *inputs(tmp_path), "--html", "report.html")
    assert result.returncode == 0, result.stderr
    text = (tmp_path / "report.html").read_text()
    for anchor in (
        "summary",
        "introduction",
        "results",
        "regions",
        "inputs",
        "metrics",
        "truth",
        "discordance",
        "methods",
        "checks",
    ):
        assert f"id='{anchor}'" in text, anchor  # fmt: skip
    # plots are inline SVG; no external scripts
    assert text.count("<svg") >= 4 and "<script" not in text
    assert "All 1 consistency checks passed" in text
    # every discordance class is explained, with its counts
    assert "low_coverage: 1 FN, 0 FP" in text and "no_evidence: 0 FN, 1 FP" in text
    # region flow with the readable step name; tool versions in the methods
    assert "intersection of all target BEDs" in text and "rtg-tools 3.12.1" in text


def test_failed_check_banner(tmp_path):
    result = run(tmp_path, *inputs(tmp_path, check_status="FAIL"), "--html", "report.html")
    assert result.returncode == 0, result.stderr
    assert "1 of 1 consistency checks FAILED" in (tmp_path / "report.html").read_text()


def test_datavzrd_tables_and_links(tmp_path):
    result = run(tmp_path, *inputs(tmp_path), "--datavzrd", "dv")
    assert result.returncode == 0, result.stderr
    config = json.loads((tmp_path / "dv" / "datavzrd.yaml").read_text())
    assert config["default-view"] == "truth-variants"
    # one table per discordance class that occurs, linked to the truth variants
    assert {"discordant-low-coverage", "discordant-no-evidence"} <= set(config["views"])
    link = config["datasets"]["compared-variants"]["links"]["truth variant"]
    assert link == {
        "column": "variant_id",
        "table-row": "truth-variants/variant_id",
        "optional": True,
    }
    for dataset in config["datasets"].values():
        assert (tmp_path / "dv" / dataset["path"]).exists()
    header = (tmp_path / "dv" / "truth-variants.tsv").read_text().splitlines()[1]
    assert header.startswith("S1:chr1:100:A>C\t")


def test_needs_an_output(tmp_path):
    result = run(tmp_path, *inputs(tmp_path))
    assert result.returncode == 1 and "nothing to write" in result.stderr


def test_read_versions(tmp_path):
    path = tmp_path / "v.yml"
    path.write_text("A:\n  bcftools: 1.22\n  'samtools': '1.21'\n\nB:\n  python: 3.12.7\n")
    assert report_html.read_versions(path) == {
        "A": {"bcftools": "1.22", "samtools": "1.21"},
        "B": {"python": "3.12.7"},
    }
