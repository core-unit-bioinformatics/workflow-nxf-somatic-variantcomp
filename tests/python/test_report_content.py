"""Unit tests for bin/report_content.py (MultiQC sections and datavzrd tables)."""

import json
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "report_content.py"
sys.path.insert(0, str(SCRIPT.parent))

import report_content  # noqa: E402


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
        [
            ["S1", "wgs", "wes", "overall", "all", 2, 1, 2, 1, 0.667, 0.667, 0.667],
            ["S1", "wgs", "wes", "match", "exact", 2, 1, 2, 1, 0.667, 0.667, 0.667],
            ["S1", "wgs", "wes", "match", "position", 3, 0, 3, 0, 1.0, 1.0, 1.0],
            ["S1", "wgs", "wes", "variant_type", "SNV", 2, 0, 2, 1, 1.0, 0.667, 0.8],
            ["S1", "wgs", "wes", "truth_vaf", "[0.2,0.5)", 2, 1, "NA", "NA", 0.667, "NA", "NA"],
        ],
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
            ["wgs", "wes", "match", "exact", 1, 2, 1, 2, 1, 0.667, 0.667, 0.667],
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
        "truth_vaf",
        "query_vaf",
        "reasons",
        "gene",
        "hgvsp",
        "consequence",
        "consequence_group",
        "impact",
        "dbsnp",
    ]
    variants = tsv(
        tmp_path / "v.tsv",
        variant_columns,
        [
            ["S1", "wgs", "wes", "TP", "concordant", "chr1", 100, "A", "C", "SNV", 0.3, 0.3,
             "concordant", "KRAS", "p.Gly12Val", "missense_variant", "non_synonymous", "MODERATE",
             "rs121913529"],
            ["S1", "wgs", "wes", "TP", "concordant", "chr1", 200, "G", "T", "SNV", 0.4, 0.4,
             "concordant", ".", ".", "intergenic_variant", "other", "MODIFIER", "."],
            ["S1", "wgs", "wes", "FN", "low_coverage", "chr1", 300, "AT", "A", "DEL", 0.1, "NA",
             "low_coverage,indel,in_repeat", "TP53", "p.Arg175fs", "frameshift_variant",
             "non_synonymous", "HIGH", "."],
            ["S1", "wgs", "wes", "FP", "no_evidence", "chr1", 400, "C", "G", "SNV", "NA", 0.02,
             "no_evidence", "TP53", ".", "intron_variant", "other", "MODIFIER", "."],
        ],
    )  # fmt: skip
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
        ["truth", "query", "status", *report_content.CLASS_ORDER],
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
    return [
        "--per-pair", per_pair, "--metrics", metrics, "--variants", variants,
        "--truth-variants", truth, "--classes", classes, "--accounting", accounting,
        "--removed", removed, "--region-steps", steps, "--checks", checks,
        "--samplesheet", str(sheet),
    ]  # fmt: skip


def run(tmp_path, *args):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args], capture_output=True, text=True, cwd=tmp_path
    )


def test_multiqc_sections(tmp_path):
    result = run(tmp_path, *inputs(tmp_path), "--multiqc", "mqc")
    assert result.returncode == 0, result.stderr
    files = {f.name for f in (tmp_path / "mqc").iterdir()}
    for section in (
        "summary",
        "introduction",
        "truth_found",
        "why",
        "characterization",
        "metrics",
        "regions",
        "filtering",
        "checks",
    ):
        assert f"vc_{section}_mqc.html" in files, section
    summary = (tmp_path / "mqc" / "vc_summary_mqc.html").read_text()
    # banner, one card set per matching method (toggle), the plot library
    assert "All 1 consistency checks passed" in summary
    assert summary.count("data-vc-group='summary_method'") == 3
    assert "window.vcLib" in summary and "mqc_hidesamples" in summary
    # recall vs precision: native MultiQC scatter plot, axes 0-1, one dataset per method
    scatter = json.loads((tmp_path / "mqc" / "vc_recall_precision_mqc.json").read_text())
    assert scatter["plot_type"] == "scatter"
    assert (scatter["pconfig"]["xmin"], scatter["pconfig"]["xmax"]) == (0, 1)
    assert [d["name"] for d in scatter["pconfig"]["data_labels"]] == [
        "RTG vcfeval",
        "isec: same allele",
        "isec: same position",
    ]
    assert scatter["data"][2]["S1 (wes)"] == {
        "x": 1.0,
        "y": 1.0,
        "color": report_content.PALETTE[0],
    }
    # truth variants found: vertical bars with the label found/total
    assert '"2/3"' in (tmp_path / "mqc" / "vc_truth_found_mqc.html").read_text()
    # discordant variants in detail: both directions, several reasons, protein-altering
    detail = (tmp_path / "mqc" / "vc_characterization_mqc.html").read_text()
    assert (
        "wgs only</b> (missed by wes, FN)" in detail
        and "1 of 1 have more than one reason" in detail
    )
    assert "TP53 p.Arg175fs" in detail
    # region steps in Mbp/kbp, not bp
    assert "10.0 kbp" in (tmp_path / "mqc" / "vc_regions_mqc.html").read_text()  # fmt: skip


def test_failed_check_banner(tmp_path):
    result = run(tmp_path, *inputs(tmp_path, check_status="FAIL"), "--multiqc", "mqc")
    assert result.returncode == 0, result.stderr
    assert (
        "1 of 1 consistency checks FAILED" in (tmp_path / "mqc" / "vc_summary_mqc.html").read_text()
    )


def test_datavzrd_tables_and_links(tmp_path):
    result = run(tmp_path, *inputs(tmp_path), "--datavzrd", "dv")
    assert result.returncode == 0, result.stderr
    config = json.loads((tmp_path / "dv" / "datavzrd.yaml").read_text())
    assert config["default-view"] == "multiqc-report"
    views = config["views"]
    # the MultiQC report as a page of the table browser
    assert views["multiqc-report"]["render-html"] == {"script-path": "multiqc_page.js"}
    assert (
        "../../../multiqc/multiqc_report.html" in (tmp_path / "dv" / "multiqc_page.js").read_text()
    )
    # one table of all discordant variants, linked to the truth variants and the genes
    links = config["datasets"]["discordant-variants"]["links"]
    assert links["truth variant"] == {
        "column": "variant_id",
        "table-row": "truth-variants/variant_id",
        "optional": True,
    }
    assert links["gene table"]["table-row"] == "genes/gene"
    # every view starts with a description that links back to the MultiQC report
    assert all("[MultiQC report](" in v["desc"] for k, v in views.items() if k != "multiqc-report")
    # columns: variant name first, coloured status/class, links to dbSNP and Ensembl
    columns = views["discordant-variants"]["render-table"]["columns"]
    assert columns["status"]["plot"]["heatmap"]["domain"] == ["TP", "FN", "FP"]
    assert (
        "dbSNP" in columns["dbsnp"]["link-to-url"] and "Ensembl" in columns["locus"]["link-to-url"]
    )
    assert columns["reasons"]["plot"]["pills"]["separator"] == ","
    for dataset in config["datasets"].values():
        assert (tmp_path / "dv" / dataset["path"]).exists()
    header = (tmp_path / "dv" / "discordant-variants.tsv").read_text().splitlines()[0].split("\t")
    assert header[:4] == ["variant", "gene", "hgvsp", "consequence"] and header[-1] == "variant_id"
    # the main reason is the first entry of reasons: no column class in the tables
    assert "class" not in header and "reasons" in header
    rows = (tmp_path / "dv" / "discordant-variants.tsv").read_text().splitlines()
    assert rows[1].startswith("TP53 p.Arg175fs\tTP53\t")  # fmt: skip


def test_variant_label():
    row = {"chrom": "chr1", "pos": 5, "ref": "A", "alt": "C"}
    assert report_content.variant_label({**row, "gene": "KRAS", "hgvsp": "p.G12V"}) == "KRAS p.G12V"
    assert (
        report_content.variant_label({**row, "gene": "KRAS", "hgvsp": ".", "hgvsc": "c.35G>T"})
        == "KRAS c.35G>T"
    )
    assert (
        report_content.variant_label({**row, "gene": "KRAS", "consequence": "intron_variant&x"})
        == "KRAS intron_variant"
    )
    assert report_content.variant_label({**row, "gene": "."}) == "chr1:5 A>C"


def test_fmt_bp():
    assert report_content.fmt_bp(3_100_000_000) == "3,100.00 Mbp"
    assert report_content.fmt_bp(1_987_915) == "1.99 Mbp"
    assert report_content.fmt_bp(125_691) == "125.7 kbp"
    assert report_content.fmt_bp(500) == "500 bp"


def test_html_table_shows_at_most_max_rows():
    import pandas as pd

    df = pd.DataFrame({"a": range(40)})
    text = report_content.html_table(df, link="x.html")
    assert text.count("<tr>") == report_content.MAX_ROWS + 1  # header + rows
    assert "Showing 15 of 40 rows" in text


def test_needs_an_output(tmp_path):
    result = run(tmp_path, *inputs(tmp_path))
    assert result.returncode == 1 and "nothing to write" in result.stderr


def test_venn2_shows_truth_only_both_query_only():
    svg = report_content.venn2_svg("wgs", "wes", 3, 10, 2)
    assert ">3<title>only in wgs (FN): 3</title>" in svg
    assert ">10<title>in both: 10 truth variants found by wes</title>" in svg
    assert ">2<title>only in wes (FP): 2</title>" in svg
    # one truth record matched by two query records (vcfeval): calls in brackets
    assert ">10 (11)<" in report_content.venn2_svg("wgs", "wes", 3, 10, 2, 11)


def test_overlap_section_has_one_venn_per_sample_and_method(tmp_path):
    result = run(tmp_path, *inputs(tmp_path), "--multiqc", "mqc")
    assert result.returncode == 0, result.stderr
    text = (tmp_path / "mqc" / "vc_overlap_mqc.html").read_text()
    # three methods (toggle), one diagram each for the sample
    assert text.count("data-vc-group='overlap_method'") == 3
    assert text.count("<svg") == 3
    # isec position: 3 TP, 0 FN, 0 FP in the fixture
    assert ">0<title>only in wgs (FN): 0</title>" in text
