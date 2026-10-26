"""Unit tests for bin/report_bundle.py (start page and report zip)."""

import subprocess
import sys
import zipfile
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "report_bundle.py"


def run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)


def reports(tmp_path):
    """A MultiQC report and a datavzrd folder with a MultiQC page and one table view."""
    multiqc = tmp_path / "multiqc_report.html"
    multiqc.write_text("<html>mqc</html>")
    tables = tmp_path / "tables"
    (tables / "multiqc-report").mkdir(parents=True)
    (tables / "truth-variants").mkdir()
    (tables / "index.html").write_text("index")
    (tables / "multiqc-report" / "index_1.html").write_text("page")
    (tables / "truth-variants" / "index_1.html").write_text("table")
    return str(multiqc), str(tables)


def test_start_page_and_zip_keep_the_output_layout(tmp_path):
    multiqc, tables = reports(tmp_path)
    zip_path, main = tmp_path / "r.zip", tmp_path / "main.html"
    proc = run(
        "--multiqc",
        multiqc,
        "--tables",
        tables,
        "--output",
        str(zip_path),
        "--start-page",
        str(main),
    )
    assert proc.returncode == 0, proc.stderr
    # the start page opens the MultiQC page of the table browser
    assert 'url=report/tables/multiqc-report/index_1.html"' in main.read_text()
    names = set(zipfile.ZipFile(zip_path).namelist())
    assert {
        "variantcomp_report/main.html",
        "variantcomp_report/multiqc/multiqc_report.html",
        "variantcomp_report/report/tables/index.html",
        "variantcomp_report/report/tables/multiqc-report/index_1.html",
        "variantcomp_report/report/tables/truth-variants/index_1.html",
    } == names
    # the page the start page points to is in the zip, relative to the start page
    assert (
        zipfile.ZipFile(zip_path).read("variantcomp_report/main.html").decode() == main.read_text()
    )


def test_missing_tables_folder_is_an_error(tmp_path):
    multiqc, _ = reports(tmp_path)
    proc = run(
        "--multiqc",
        multiqc,
        "--tables",
        str(tmp_path / "none"),
        "--output",
        str(tmp_path / "r.zip"),
        "--start-page",
        str(tmp_path / "m.html"),
    )
    assert proc.returncode == 1 and "no index.html" in proc.stderr
