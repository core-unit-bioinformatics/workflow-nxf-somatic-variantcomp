#!/usr/bin/env python3
"""Pack the reports of a run into one zip file that can be sent and opened offline.

The MultiQC report is one self-contained HTML file; the datavzrd table browser is a folder
of HTML, JavaScript and data files. Both link to each other with relative paths that
assume the layout of the pipeline output folder (multiqc/ next to report/tables/). The
zip keeps exactly this layout under one top folder and adds an index.html that opens the
MultiQC report, so the links work after unzipping.

INPUT
  --multiqc   the MultiQC report (multiqc_report.html)
  --tables    the datavzrd output folder (with index.html)
  --name      name of the top folder in the zip (default variantcomp_report)
OUTPUT
  --output    zip file, e.g. variantcomp_report.zip, with
                <name>/index.html                     start page (opens the MultiQC report)
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

INDEX_HTML = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta http-equiv="refresh" content="0; url=multiqc/multiqc_report.html">
<title>variantcomp report</title></head>
<body style="font-family: sans-serif; margin: 40px">
<h1>variantcomp report</h1>
<ul>
<li><a href="multiqc/multiqc_report.html">MultiQC report</a>: summary, results and quality control</li>
<li><a href="report/tables/index.html">Variant tables (datavzrd)</a>: every variant, linked tables</li>
</ul>
<p>Keep the folders together: the reports link to each other.</p>
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
        zf.writestr(f"{name}/index.html", INDEX_HTML)
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
    args = parser.parse_args()
    check_inputs(args.multiqc, args.tables)
    n = write_zip(args.multiqc, args.tables, args.name, args.output)
    size = os.path.getsize(args.output) / 1e6
    print(f"INFO: {args.output}: {n} files, {size:.1f} MB", file=sys.stderr)


if __name__ == "__main__":
    main()
