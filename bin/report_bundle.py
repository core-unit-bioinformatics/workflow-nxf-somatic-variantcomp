#!/usr/bin/env python3
"""The start page of the reports and one zip file that can be sent and opened offline.

The MultiQC report is one self-contained HTML file; the datavzrd table browser is a folder
of HTML, JavaScript and data files whose first page shows the MultiQC report. Both link to
each other with relative paths that assume the layout of the pipeline output folder
(multiqc/ next to report/tables/). The start page main.html sits next to these folders and
opens the table browser on its MultiQC page, so the MultiQC report and every table are in
one window. The zip keeps exactly this layout under one top folder, with main.html at its
top, so the links work after unzipping.

INPUT
  --multiqc   the MultiQC report (multiqc_report.html)
  --tables    the datavzrd output folder (with index.html)
  --name      name of the top folder in the zip (default variantcomp_report)
OUTPUT
  --start-page  the start page, e.g. main.html (for the top of the output folder)
  --output      zip file, e.g. variantcomp_report.zip, with
                  <name>/main.html                      start page
                  <name>/multiqc/multiqc_report.html    the MultiQC report
                  <name>/report/tables/...              the datavzrd table browser

Why a script: zip is not available in every container and its options differ between
systems; Python's zipfile writes the same archive everywhere (no tool packs this layout
with a start page).
"""

import argparse
import os
import sys
import zipfile

# path of the MultiQC page of the table browser, relative to the start page
START_TARGET = "report/tables/multiqc-report/index_1.html"
MAIN_HTML = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta http-equiv="refresh" content="0; url={START_TARGET}">
<title>variantcomp report</title></head>
<body style="font-family: sans-serif; margin: 40px">
<h1>variantcomp report</h1>
<p>Opening the report ... If nothing happens, use a link:</p>
<ul>
<li><a href="{START_TARGET}">Report with the variant tables</a>: the MultiQC report inside the
table browser (datavzrd); the menu at the top switches to the tables</li>
<li><a href="multiqc/multiqc_report.html">MultiQC report</a> on its own</li>
<li><a href="report/tables/index.html">Variant tables</a> on their own</li>
</ul>
<p>Keep the folders multiqc/ and report/ next to this page: the reports link to each other.</p>
</body></html>
"""


def check_inputs(multiqc, tables):
    """Stop with an error if the MultiQC report or the datavzrd folder is missing."""
    if not os.path.isfile(multiqc):
        sys.exit(f"ERROR: --multiqc {multiqc}: file not found; give the MultiQC report")
    if not os.path.isfile(os.path.join(tables, "index.html")):
        sys.exit(f"ERROR: --tables {tables}: no index.html; give the datavzrd output folder")


def write_zip(multiqc, tables, name, output):
    """Write the zip; returns the number of files packed."""
    n = 0
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        # start page and MultiQC report
        zf.writestr(f"{name}/main.html", MAIN_HTML)
        zf.write(multiqc, f"{name}/multiqc/multiqc_report.html")
        n += 2
        # every file of the datavzrd folder, with its path relative to the folder
        for root, _dirs, files in os.walk(tables):
            for file in sorted(files):
                path = os.path.join(root, file)
                rel = os.path.relpath(path, tables)
                zf.write(path, f"{name}/report/tables/{rel}")
                n += 1
    return n


def main():
    """Parse the command line and write the zip."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--multiqc", required=True, help="MultiQC report (HTML)")
    parser.add_argument("--tables", required=True, help="datavzrd output folder")
    parser.add_argument("--name", default="variantcomp_report", help="top folder in the zip")
    parser.add_argument("--output", required=True, help="output zip file")
    parser.add_argument("--start-page", required=True, help="output: start page (main.html)")
    args = parser.parse_args()
    check_inputs(args.multiqc, args.tables)
    # the start page for the output folder, then the zip with the same page at its top
    with open(args.start_page, "w") as fh:
        fh.write(MAIN_HTML)
    n = write_zip(args.multiqc, args.tables, args.name, args.output)
    size = os.path.getsize(args.output) / 1e6
    print(f"INFO: {args.output}: {n} files, {size:.1f} MB", file=sys.stderr)


if __name__ == "__main__":
    main()
