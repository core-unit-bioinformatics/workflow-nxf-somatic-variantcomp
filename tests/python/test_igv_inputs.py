"""Unit tests for bin/igv_inputs.py (inputs of the IGV views, --igv_report)."""

import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "bin" / "igv_inputs.py"
TEMPLATE = ROOT / "assets" / "igv" / "variant_template.html"

sys.path.insert(0, str(ROOT / "bin"))
import igv_inputs  # noqa: E402
import report_content  # noqa: E402

VARIANT_COLUMNS = ["sample", "truth", "query", "status", "chrom", "pos", "ref", "alt", "type"]
VARIANT_COLUMNS += ["reasons", "gene"]


def run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)


def tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return str(path)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fai(tmp_path):
    """Reference with chr2 before chr10 (sort by reference order, not by name)."""
    path = tmp_path / "ref.fa.fai"
    path.write_text("chr2\t1000\t6\t60\t61\nchr10\t1000\t1030\t60\t61\n")
    return str(path)


def variants(tmp_path):
    """Two queries (wes, panel) of truth wgs; S1 has one TP in both, an FN in wes only, an FP
    in panel only; S2 has one TP."""
    return tsv(
        tmp_path / "variants.tsv",
        VARIANT_COLUMNS,
        [
            ["S1", "wgs", "wes", "TP", "chr10", 50, "A", "C", "SNV", "concordant", "GENE1"],
            ["S1", "wgs", "panel", "TP", "chr10", 50, "A", "C", "SNV", "concordant", "GENE1"],
            ["S1", "wgs", "wes", "FN", "chr2", 300, "AT", "A", "DEL", "low_coverage,indel", "."],
            ["S1", "wgs", "panel", "TP", "chr2", 300, "AT", "A", "DEL", "concordant", "."],
            ["S1", "wgs", "panel", "FP", "chr2", 100, "G", "T", "SNV", "no_evidence", "GENE2"],
            ["S2", "wgs", "wes", "TP", "chr2", 10, "C", "G", "SNV", "concordant", "NA"],
        ],
    )


def required(tmp_path):
    """A required variant of S1, compared (chr10:50) in wgs and wes; one not compared."""
    return tsv(
        tmp_path / "required.tsv",
        ["list", "sample", "technology", "label", "gene", "chrom", "pos", "ref", "alt", "status"],
        [
            ["clin", "S1", "wgs", "v1", "GENE1", "chr10", 50, "A", "C", "PASS"],
            ["clin", "S1", "wes", "v1", "GENE1", "chr10", 50, "A", "C", "PASS"],
            ["clin", "S1", "wgs", "v2", "GENE3", "chr2", 500, "C", "CTT", "FAIL_not_called"],
            ["clin", "S1", "wes", "v2", "GENE3", "chr2", 500, "C", "CTT", "FAIL_not_called"],
        ],
    )


def test_sites_one_row_per_variant_sorted_and_split(tmp_path):
    out, links = tmp_path / "sites", tmp_path / "links.tsv"
    proc = run(
        "sites",
        "--variants",
        variants(tmp_path),
        "--required",
        required(tmp_path),
        "--fai",
        fai(tmp_path),
        "--max-sites",
        "2",
        "--outdir",
        str(out),
        "--links",
        str(links),
    )
    assert proc.returncode == 0, proc.stderr
    # S1: 4 sites (3 compared, 1 required only) in 2 parts of at most 2 sites; S2: 1 part
    assert sorted(p.name for p in out.rglob("*.sites.tsv")) == [
        "S1_001.sites.tsv",
        "S1_002.sites.tsv",
        "S2_001.sites.tsv",
    ]
    part1 = read_tsv(out / "S1" / "S1_001.sites.tsv")
    part2 = read_tsv(out / "S1" / "S1_002.sites.tsv")
    # reference order: chr2 before chr10; within a contig by position
    assert [r["variant_id"] for r in part1 + part2] == [
        "S1:chr2:100:G>T",
        "S1:chr2:300:AT>A",
        "S1:chr2:500:C>CTT",
        "S1:chr10:50:A>C",
    ]
    header = list(part1[0])
    assert header[:6] == ["variant_id", "chrom", "pos", "end", "gene", "type"]
    assert header[6:] == [
        "status_wes",
        "reasons_wes",
        "status_panel",
        "reasons_panel",
        "known_required",
        "required_wgs",
        "required_wes",
    ]
    fp, deletion = part1
    # status per comparison, '.' where the variant is not in the comparison
    assert (fp["status_wes"], fp["status_panel"], fp["reasons_panel"]) == (".", "FP", "no_evidence")
    assert (deletion["status_wes"], deletion["reasons_wes"]) == ("FN", "low_coverage,indel")
    # end: last reference base (1-based, inclusive)
    assert (deletion["pos"], deletion["end"]) == ("300", "301")
    required_only, compared_required = part2
    assert required_only["type"] == "INS" and required_only["status_wes"] == "."
    assert required_only["known_required"] == "yes"
    assert required_only["required_wes"] == "FAIL_not_called"
    assert compared_required["known_required"] == "yes" and compared_required["gene"] == "GENE1"
    assert fp["known_required"] == "no" and fp["required_wgs"] == "."
    # links: every site once, file name of the IGVREPORTS output
    rows = read_tsv(links)
    assert len(rows) == 5 and len({r["variant_id"] for r in rows}) == 5
    assert {r["igv_file"] for r in rows if r["sample"] == "S1"} == {
        "S1_001_report.html",
        "S1_002_report.html",
    }
    assert rows[-1] == {
        "variant_id": "S2:chr2:10:C>G",
        "sample": "S2",
        "part": "001",
        "igv_file": "S2_001_report.html",
    }


def test_variant_ids_match_the_report_tables(tmp_path):
    """The links use exactly the variant_id of the datavzrd tables."""
    path = variants(tmp_path)
    expected = set(report_content.add_variant_columns(report_content.read(path))["variant_id"])
    links = tmp_path / "links.tsv"
    proc = run(
        "sites",
        "--variants",
        path,
        "--fai",
        fai(tmp_path),
        "--outdir",
        str(tmp_path / "s"),
        "--links",
        str(links),
    )
    assert proc.returncode == 0, proc.stderr
    assert {r["variant_id"] for r in read_tsv(links)} == expected


def test_sites_errors(tmp_path):
    base = ["sites", "--fai", fai(tmp_path), "--outdir", str(tmp_path / "s")]
    base += ["--links", str(tmp_path / "l.tsv")]
    # a contig of no reference
    bad = tsv(
        tmp_path / "bad.tsv",
        VARIANT_COLUMNS,
        [["S1", "wgs", "wes", "TP", "chrX", 5, "A", "C", "SNV", "concordant", "."]],
    )
    proc = run(*base, "--variants", bad)
    assert proc.returncode == 1 and "contigs not in the reference index: chrX" in proc.stderr
    # a missing column
    bad = tsv(tmp_path / "bad2.tsv", ["sample", "chrom"], [["S1", "chr2"]])
    proc = run(*base, "--variants", bad)
    assert proc.returncode == 1 and "missing column(s) truth" in proc.stderr
    # the initial view larger than the embedded region
    proc = run(*base, "--variants", variants(tmp_path), "--flanking", "50", "--window", "100")
    assert proc.returncode == 1 and "set --igv_window <= --igv_flanking" in proc.stderr
    # the same variant twice in one comparison
    dup = tsv(
        tmp_path / "dup.tsv",
        VARIANT_COLUMNS,
        [["S1", "wgs", "wes", "TP", "chr2", 5, "A", "C", "SNV", "concordant", "."]] * 2,
    )
    proc = run(*base, "--variants", dup)
    assert proc.returncode == 1 and "occurs twice in comparison wgs vs wes" in proc.stderr


def test_sites_without_variants_writes_an_empty_links_table(tmp_path):
    empty = tsv(tmp_path / "empty.tsv", VARIANT_COLUMNS, [])
    links = tmp_path / "links.tsv"
    proc = run(
        "sites",
        "--variants",
        empty,
        "--fai",
        fai(tmp_path),
        "--outdir",
        str(tmp_path / "s"),
        "--links",
        str(links),
    )
    assert proc.returncode == 0, proc.stderr
    assert links.read_text() == "variant_id\tsample\tpart\tigv_file\n"
    assert not list(tmp_path.rglob("*.sites.tsv"))


def test_all_pairs_columns_name_truth_and_query(tmp_path):
    rows = [
        ["S1", "wgs", "wes", "TP", "chr2", 5, "A", "C", "SNV", "concordant", "."],
        ["S1", "wes", "panel", "FN", "chr2", 5, "A", "C", "SNV", "absent", "."],
    ]
    path = tsv(tmp_path / "v.tsv", VARIANT_COLUMNS, rows)
    out = tmp_path / "s"
    proc = run(
        "sites",
        "--variants",
        path,
        "--fai",
        fai(tmp_path),
        "--outdir",
        str(out),
        "--links",
        str(tmp_path / "l.tsv"),
    )
    assert proc.returncode == 0, proc.stderr
    row = read_tsv(out / "S1" / "S1_001.sites.tsv")[0]
    assert row["status_wgs_vs_wes"] == "TP" and row["status_wes_vs_panel"] == "FN"


def sheet(tmp_path, wes_library="capture"):
    path = tmp_path / "sheet.csv"
    path.write_text(
        "sample,technology,library_type,tumor_id,normal_id,caller,is_truth,vcf,tumor_bam,"
        "normal_bam,target_bed\n"
        f"S1,wes,{wes_library},T1,N1,x,false,a.vcf.gz,t.bam,n.bam,\n"
        "S1,wgs,wgs,T2,N2,x,true,b.vcf.gz,t.bam,n.bam,\n"
        "S2,wgs,wgs,T,N,x,true,c.vcf.gz,t.bam,n.bam,\n"
    )
    return str(path)


TRACK_ARGS = [
    "--labels",
    "wes:calls",
    "wes:normal",
    "wes:tumor",
    "wgs:calls",
    "wgs:normal",
    "wgs:tumor",
    "--files",
    "tracks/01/wes.vcf.gz",
    "tracks/02/N.bam",
    "tracks/03/T.bam",
    "tracks/04/wgs.vcf.gz",
    "tracks/05/N.cram",
    "tracks/06/T.bam",
]


def test_tracks_truth_first_tumor_before_normal(tmp_path):
    out, flags = tmp_path / "t.json", tmp_path / "flags.txt"
    proc = run(
        "tracks",
        "--sample",
        "S1",
        "--samplesheet",
        sheet(tmp_path),
        *TRACK_ARGS,
        "--sampling-depth",
        "300",
        "--output",
        str(out),
        "--flags-output",
        str(flags),
    )
    assert proc.returncode == 0, proc.stderr
    config = json.loads(out.read_text())
    assert [(t["name"], t["url"], t["order"]) for t in config] == [
        ("WGS tumor (T2)", "tracks/06/T.bam", 1),
        ("WGS normal (N2)", "tracks/05/N.cram", 2),
        ("WES tumor (T1)", "tracks/03/T.bam", 3),
        ("WES normal (N1)", "tracks/02/N.bam", 4),
        ("WGS calls (normalized)", "tracks/04/wgs.vcf.gz", 5),
        ("WES calls (normalized)", "tracks/01/wes.vcf.gz", 6),
    ]
    assert config[1]["format"] == "cram" and config[4]["type"] == "variant"
    assert all(t["samplingDepth"] == 300 for t in config[:4])
    # no amplicon library: duplicates are hidden and not embedded (1024 + 512)
    assert all(t["filter"]["duplicate"] for t in config[:4])
    assert flags.read_text() == "1536\n"


def test_tracks_amplicon_keeps_duplicates(tmp_path):
    out, flags = tmp_path / "t.json", tmp_path / "flags.txt"
    args = ["tracks", "--sample", "S1", "--samplesheet", sheet(tmp_path, "amplicon"), *TRACK_ARGS]
    proc = run(*args, "--output", str(out), "--flags-output", str(flags))
    assert proc.returncode == 0, proc.stderr
    filters = {t["name"]: t["filter"]["duplicate"] for t in json.loads(out.read_text())[:4]}
    # duplicates are embedded (only QC-failed reads removed) and shown for the amplicon
    # library only
    assert filters == {
        "WGS tumor (T2)": True,
        "WGS normal (N2)": True,
        "WES tumor (T1)": False,
        "WES normal (N1)": False,
    }
    assert flags.read_text() == "512\n"
    # --amplicon-keep-duplicates false: duplicates hidden everywhere
    proc = run(
        *args,
        "--amplicon-keep-duplicates",
        "false",
        "--output",
        str(out),
        "--flags-output",
        str(flags),
    )
    assert proc.returncode == 0, proc.stderr
    assert flags.read_text() == "1536\n"


def test_tracks_errors(tmp_path):
    out = ["--output", str(tmp_path / "t.json"), "--flags-output", str(tmp_path / "f.txt")]
    # a missing normal alignment
    labels = ["--labels", "wgs:tumor", "wgs:calls", "--files", "T.bam", "c.vcf.gz"]
    proc = run("tracks", "--sample", "S2", "--samplesheet", sheet(tmp_path), *labels, *out)
    assert proc.returncode == 1 and "no normal alignment of technology wgs" in proc.stderr
    # unknown sample
    proc = run("tracks", "--sample", "S9", "--samplesheet", sheet(tmp_path), *labels, *out)
    assert proc.returncode == 1 and "sample 'S9' not found" in proc.stderr
    # labels and files differ in number
    proc = run(
        "tracks",
        "--sample",
        "S2",
        "--samplesheet",
        sheet(tmp_path),
        "--labels",
        "wgs:tumor",
        "--files",
        "a",
        "b",
        *out,
    )
    assert proc.returncode == 1 and "1 labels for 2 track files" in proc.stderr
    # a bad role
    proc = run(
        "tracks",
        "--sample",
        "S2",
        "--samplesheet",
        sheet(tmp_path),
        "--labels",
        "wgs:germline",
        "--files",
        "a",
        *out,
    )
    assert proc.returncode == 1 and "expected <technology>:<tumor|normal|calls>" in proc.stderr


def test_template_inlines_igv_js(tmp_path):
    js = tmp_path / "igv.min.js"
    js.write_text("var igv = {};\n//# sourceMappingURL=igv.min.js.map\n")
    out = tmp_path / "template.html"
    proc = run("template", "--template", str(TEMPLATE), "--igv-js", str(js), "--output", str(out))
    assert proc.returncode == 0, proc.stderr
    text = out.read_text()
    # no CDN script, the code inlined, the source-map comment removed
    assert "cdn.jsdelivr.net" not in text and "var igv = {};" in text
    assert "sourceMappingURL" not in text
    # the placeholders of create_report are kept exactly once each (none in comments)
    assert text.count('"@TABLE_JSON@"') == 1 and text.count('"@SESSION_DICTIONARY@"') == 1
    # deep link to the row of #<variant_id>
    assert "selectFromHash" in text and "history.replaceState" in text


def test_template_errors(tmp_path):
    js = tmp_path / "igv.min.js"
    js.write_text("var a = '</script>';\n")
    out = ["--output", str(tmp_path / "o.html")]
    proc = run("template", "--template", str(TEMPLATE), "--igv-js", str(js), *out)
    assert proc.returncode == 1 and "contains '</script'" in proc.stderr
    js.write_text("var igv = {};\n")
    no_cdn = tmp_path / "t.html"
    no_cdn.write_text("<html></html>\n")
    proc = run("template", "--template", str(no_cdn), "--igv-js", str(js), *out)
    assert proc.returncode == 1 and "found 0" in proc.stderr


def test_index_lists_every_part(tmp_path):
    links = tsv(
        tmp_path / "links.tsv",
        ["variant_id", "sample", "part", "igv_file"],
        [
            ["S1:chr1:1:A>C", "S1", "001", "S1_001_report.html"],
            ["S1:chr1:2:A>C", "S1", "001", "S1_001_report.html"],
            ["S1:chr1:3:A>C", "S1", "002", "S1_002_report.html"],
        ],
    )
    out = tmp_path / "index.html"
    proc = run("index", "--links", links, "--max-sites", "2", "--output", str(out))
    assert proc.returncode == 0, proc.stderr
    text = out.read_text()
    assert "<td>S1</td><td>001</td><td class='n'>2</td>" in text
    assert "href='S1_002_report.html'" in text and "../../multiqc/multiqc_report.html" in text


def test_variant_type():
    assert igv_inputs.variant_type("A", "C") == "SNV"
    assert igv_inputs.variant_type("AC", "GT") == "MNV"
    assert igv_inputs.variant_type("A", "AT") == "INS"
    assert igv_inputs.variant_type("AT", "A") == "DEL"
    assert igv_inputs.variant_type("AT", "C") == "COMPLEX"
