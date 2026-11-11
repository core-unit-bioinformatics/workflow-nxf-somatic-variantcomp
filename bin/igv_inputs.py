#!/usr/bin/env python3
"""Inputs of the IGV views (--igv_report, decision D40): sites, tracks, template, index.

igv-reports (create_report) writes one self-contained HTML file per sites table: a table
of the sites and an igv.js browser with the reads of every alignment around each site.
All read data are embedded, so the reports are split: one report per sample and part of
at most --max-sites sites. This script prepares the inputs of create_report and the
overview page; create_report itself runs in the nf-core module IGVREPORTS.

Subcommands
  sites     the sites of every sample: every compared variant (TP, FN, FP of every
            comparison, the rows of the datavzrd table compared-variants) and every
            required known variant, one row per variant_id; split into parts; the links
            table variant_id -> IGV report file (used by the datavzrd tables, the index and
            the consistency checks)
  tracks    the track configuration (igv.js track objects, create_report --track-config)
            of one sample, and the SAM flags excluded when the reads are embedded
  template  the report template with igv.js inlined (works offline, no CDN)
  index     the overview page report/igv/index.html (every sample and part)

INPUT
  sites     --variants  discordance/variants.annotated.tsv (TSV; columns sample, truth,
                        query, status, chrom, pos (1-based), ref, alt, type, reasons;
                        optional gene)
            --required  known_variants/required_variants.tsv (optional; columns sample,
                        technology, chrom, pos (1-based), ref, alt, status; optional gene)
            --fai       reference .fai (contig order of the sort)
  tracks    --samplesheet  pipeline_info/samplesheet.validated.csv (CSV; sample,
                           technology, library_type, tumor_id, normal_id, is_truth)
            --labels       one label per track file, <technology>:<tumor|normal|calls>
            --files        the track files as staged for create_report (same order)
  template  --template  assets/igv/variant_template.html (pipeline template, from
                        igv-reports 1.16.0 templates/variant_template.html, MIT)
            --igv-js    igv.js (igv.min.js) to inline
  index     --links     igv_links.tsv of `sites`
OUTPUT
  sites     --outdir/<sample>/<sample>_<part>.sites.tsv  TSV, one row per site, sorted by
              contig (reference order), pos, ref, alt; columns:
                variant_id   sample:chrom:pos:ref>alt (as in the datavzrd tables)
                chrom        contig (reference naming)
                pos          first base of REF (1-based)
                end          last base of REF (1-based, inclusive: pos + len(REF) - 1)
                gene         gene symbol (VEP; '.' if none) - if the variant table has it
                type         SNV, MNV, INS, DEL or COMPLEX
                status_<c>   TP, FN, FP of comparison <c> ('.': not in that comparison)
                reasons_<c>  reasons of comparison <c> (main reason first; concordant for TP)
                known_required  yes if the variant is a required known variant, else no
                required_<technology>  status of the required variant (PASS, FAIL_*; only
                                       with required variants; '.' for other sites)
              <c> is the query technology; with several truth technologies
              (--comparison_mode all_pairs) <truth>_vs_<query>.
              <part> is 1-based and zero-padded (001, 002, ...). A sample without sites
              gets no file (INFO message).
            --links  TSV: variant_id, sample, part, igv_file (<sample>_<part>_report.html,
                     the name of the IGVREPORTS output)
  tracks    --output         JSON array of igv.js track objects (alignments: tumor, then
                             normal, of every technology, truth technology first, then the
                             samplesheet order; then the normalized calls of every
                             technology in the same order)
            --flags-output   one number: the SAM flags excluded by samtools when the reads
                             are embedded (create_report --exclude-flags)
  template  --output   HTML template for create_report --template
  index     --output   HTML page

Why a script: create_report needs a sites table, a track configuration and a template;
no tool builds them from the pipeline tables. igv-reports' own --standalone downloads
igv.js from a CDN at run time (not possible offline), and its template has no deep link.
"""

import argparse
import csv
import html
import json
import os
import sys
from pathlib import Path

# readers of the reference index and the variant-id logic of the report tables (the IGV
# links must use exactly the variant_id of the datavzrd tables)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_reference import read_fai  # noqa: E402
from report_content import add_variant_columns, read  # noqa: E402

# SAM flags: 0x400 PCR/optical duplicate, 0x200 not passing quality controls
FLAG_DUPLICATE = 1024
FLAG_QC_FAIL = 512
# the line of the template that loads igv.js from a CDN; replaced by the inlined script
IGV_SCRIPT_MARK = "igv.min.js"
ROLES = ("tumor", "normal", "calls")


def log(level, message):
    """Print a message with level prefix (INFO, WARNING) to stderr."""
    print(f"{level}: {message}", file=sys.stderr)


# ----------------------------------------------------------------------------- sites


def variant_type(ref, alt):
    """SNV, MNV, INS, DEL or COMPLEX of one allele (as harmonize_vcf.py names them)."""
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    # an insertion or deletion keeps the first (anchor) base
    if len(ref) == 1 and alt.startswith(ref):
        return "INS"
    if len(alt) == 1 and ref.startswith(alt):
        return "DEL"
    return "COMPLEX"


def comparison_keys(variants):
    """{(truth, query): column suffix}: the query technology, or truth_vs_query if the
    comparisons have more than one truth technology (--comparison_mode all_pairs)."""
    pairs = list(dict.fromkeys(zip(variants["truth"], variants["query"], strict=True)))
    truths = {t for t, _ in pairs}
    return {(t, q): (q if len(truths) == 1 else f"{t}_vs_{q}") for t, q in pairs}


def require_columns(df, columns, path):
    """Stop with an error naming the file and the missing columns."""
    missing = [c for c in columns if c not in df.columns]
    if missing:
        sys.exit(
            f"ERROR: {path}: missing column(s) {', '.join(missing)} "
            f"(columns: {', '.join(map(str, df.columns))})"
        )


def compared_sites(variants, path):
    """One row per variant_id of the compared variants: chrom, pos, ref, alt, sample, gene,
    type, and status/reasons per comparison."""
    require_columns(
        variants,
        ["sample", "truth", "query", "status", "chrom", "pos", "ref", "alt", "type", "reasons"],
        path,
    )
    keys = comparison_keys(variants)
    # variant_id exactly as in the datavzrd tables (report_content.add_variant_columns)
    variants = add_variant_columns(variants)
    sites = {}
    for row in variants.to_dict("records"):
        site = sites.setdefault(
            row["variant_id"],
            {
                "variant_id": row["variant_id"],
                "sample": str(row["sample"]),
                "chrom": str(row["chrom"]),
                "pos": int(row["pos"]),
                "ref": str(row["ref"]),
                "alt": str(row["alt"]),
                "gene": ".",
                "type": str(row["type"]),
            },
        )
        # the gene of the first row that has one (VEP annotates the allele, not the row)
        gene = str(row.get("gene", "."))
        if site["gene"] == "." and gene not in ("", ".", "nan"):
            site["gene"] = gene
        key = keys[(row["truth"], row["query"])]
        # a variant_id occurs at most once per comparison (one record per allele after
        # normalization); a second row would hide a status
        if f"status_{key}" in site:
            sys.exit(
                f"ERROR: {path}: variant {row['variant_id']} occurs twice in comparison "
                f"{row['truth']} vs {row['query']}; the variant table must have one row per "
                "variant and comparison"
            )
        site[f"status_{key}"] = str(row["status"])
        site[f"reasons_{key}"] = str(row["reasons"])
    return sites, list(keys.values())


def add_required(sites, required, path):
    """Add the required known variants (new sites for the ones not compared) with their
    status per technology; returns the technologies of the required table."""
    require_columns(
        required, ["sample", "technology", "chrom", "pos", "ref", "alt", "status"], path
    )
    required = add_variant_columns(required)
    technologies = list(dict.fromkeys(required["technology"].astype(str)))
    for row in required.to_dict("records"):
        site = sites.setdefault(
            row["variant_id"],
            {
                "variant_id": row["variant_id"],
                "sample": str(row["sample"]),
                "chrom": str(row["chrom"]),
                "pos": int(row["pos"]),
                "ref": str(row["ref"]),
                "alt": str(row["alt"]),
                "gene": ".",
                # a required variant that was not compared: type from its alleles
                "type": variant_type(str(row["ref"]), str(row["alt"])),
            },
        )
        gene = str(row.get("gene", "."))
        if site["gene"] == "." and gene not in ("", ".", "nan"):
            site["gene"] = gene
        site["known_required"] = "yes"
        site[f"required_{row['technology']}"] = str(row["status"])
    return technologies


def site_table(sites, contigs, comparisons, required_techs, with_gene):
    """Rows (lists) and header of the sites table, sorted by reference contig order, pos,
    ref, alt; stops if a site is on a contig of no reference."""
    order = {c: i for i, c in enumerate(contigs)}
    unknown = sorted({s["chrom"] for s in sites.values() if s["chrom"] not in order})
    if unknown:
        sys.exit(
            f"ERROR: sites on contigs not in the reference index: {', '.join(unknown)}; the "
            "variant tables must use the reference naming"
        )
    header = ["variant_id", "chrom", "pos", "end"]
    header += ["gene"] if with_gene else []
    header += ["type"]
    header += [f"{p}_{c}" for c in comparisons for p in ("status", "reasons")]
    header += ["known_required"]
    header += [f"required_{t}" for t in required_techs]
    rows = []
    for s in sorted(
        sites.values(), key=lambda s: (order[s["chrom"]], s["pos"], s["ref"], s["alt"])
    ):
        # end: last reference base of the allele (1-based, inclusive)
        s = {**s, "end": s["pos"] + len(s["ref"]) - 1}
        s.setdefault("known_required", "no")
        rows.append((s["sample"], [str(s.get(col, ".")) for col in header]))
    return header, rows


def write_sites(args):
    """Sites tables per sample and part, and the links table."""
    if args.max_sites < 1:
        sys.exit(f"ERROR: --max-sites {args.max_sites}: must be at least 1 (--igv_max_sites)")
    # the initial view must lie inside the embedded data (flanking/2 on either side)
    if args.window is not None and args.window > args.flanking:
        sys.exit(
            f"ERROR: --window {args.window} is larger than --flanking {args.flanking}: the "
            "initial view would show bases without data; set --igv_window <= --igv_flanking"
        )
    contigs = list(read_fai(args.fai))
    variants = read(args.variants)
    log("INFO", f"{args.variants}: {len(variants)} compared variant rows")
    sites, comparisons = ({}, []) if variants.empty else compared_sites(variants, args.variants)
    n_compared = len(sites)
    required_techs = []
    if args.required:
        required = read(args.required)
        log("INFO", f"{args.required}: {len(required)} required variant rows")
        if not required.empty:
            required_techs = add_required(sites, required, args.required)
    log("INFO", f"{n_compared} compared variants, {len(sites) - n_compared} more required ones")
    with_gene = "gene" in variants.columns or any(s["gene"] != "." for s in sites.values())
    header, rows = site_table(sites, contigs, comparisons, required_techs, with_gene)
    # samples by name; parts of at most max_sites sites
    by_sample = {}
    for sample, row in rows:
        by_sample.setdefault(sample, []).append(row)
    links = []
    for sample, sample_rows in sorted(by_sample.items()):
        n_parts = (len(sample_rows) + args.max_sites - 1) // args.max_sites
        width = max(3, len(str(n_parts)))
        os.makedirs(os.path.join(args.outdir, sample), exist_ok=True)
        for i in range(n_parts):
            part = str(i + 1).zfill(width)
            chunk = sample_rows[i * args.max_sites : (i + 1) * args.max_sites]
            path = os.path.join(args.outdir, sample, f"{sample}_{part}.sites.tsv")
            with open(path, "w", newline="") as fh:
                writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                writer.writerows(chunk)
            # the report of this part: IGVREPORTS writes <prefix>_report.html
            links += [(r[0], sample, part, f"{sample}_{part}_report.html") for r in chunk]
        log("INFO", f"{sample}: {len(sample_rows)} sites in {n_parts} part(s)")
    # samples of the variant table without any site (e.g. no variant in the region)
    samples = set(variants["sample"].astype(str)) if "sample" in variants.columns else set()
    for sample in sorted(samples - set(by_sample)):
        log("INFO", f"{sample}: no sites, no IGV report")
    with open(args.links, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["variant_id", "sample", "part", "igv_file"])
        writer.writerows(links)
    log("INFO", f"{args.links}: {len(links)} links")


# ----------------------------------------------------------------------------- tracks


def read_samplesheet(path, sample):
    """The rows of one sample of the validated samplesheet, truth technology first, then
    the samplesheet order; stops if the sample has none."""
    with open(path, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["sample"] == sample]
    if not rows:
        sys.exit(f"ERROR: {path}: sample '{sample}' not found; check the sample names")
    # sorted() is stable: the samplesheet order is kept within truth / other technologies
    return sorted(rows, key=lambda r: r["is_truth"] != "true")


def track_files(labels, files):
    """{(technology, role): file} from the labels and files given in the same order."""
    if len(labels) != len(files):
        sys.exit(f"ERROR: {len(labels)} labels for {len(files)} track files; give one per file")
    tracks = {}
    for label, file in zip(labels, files, strict=True):
        technology, _, role = label.rpartition(":")
        if role not in ROLES or not technology:
            sys.exit(f"ERROR: label '{label}': expected <technology>:<{'|'.join(ROLES)}>")
        if (technology, role) in tracks:
            sys.exit(f"ERROR: label '{label}' given twice")
        tracks[(technology, role)] = file
    return tracks


def track_config(rows, tracks, args):
    """igv.js track objects of one sample and the excluded SAM flags.

    Duplicate-flagged reads are hidden, except in amplicon libraries with
    --amplicon-keep-duplicates (every read of an amplicon has the same start and end, so
    duplicate marking flags real reads; the coverage step keeps them as well). samtools
    (create_report --exclude-flags) removes reads before they are embedded: duplicates are
    only embedded if an alignment of the sample needs them (flags 512 instead of 1536);
    the igv.js filter of each alignment track then hides them in the other tracks.
    """
    keep_dups = {
        r["technology"]: args.amplicon_keep_duplicates and r["library_type"] == "amplicon"
        for r in rows
    }
    known = {r["technology"] for r in rows}
    extra = sorted({t for t, _ in tracks if t not in known})
    if extra:
        sys.exit(
            f"ERROR: tracks of technologies {', '.join(extra)} given for sample "
            f"{args.sample}, which has no such rows in {args.samplesheet}"
        )
    config = []
    # alignments: tumor before normal, truth technology first
    for r in rows:
        tech = r["technology"]
        for role, sample_id in (("tumor", r["tumor_id"]), ("normal", r["normal_id"])):
            if not sample_id or sample_id == "n/a":
                # tumor-only row: no normal alignment
                continue
            if (tech, role) not in tracks:
                sys.exit(
                    f"ERROR: sample {args.sample}: no {role} alignment of technology {tech} "
                    f"({role}_id {sample_id}); --igv_report needs tumor and normal "
                    "alignments of every row"
                )
            path = tracks[(tech, role)]
            config.append(
                {
                    "name": f"{tech.upper()} {role} ({sample_id})",
                    "url": path,
                    "format": "cram" if path.endswith(".cram") else "bam",
                    "type": "alignment",
                    "height": args.height,
                    "samplingDepth": args.sampling_depth,
                    # igv.js display filter: duplicate = hide duplicate-flagged reads,
                    # vendorFailed = hide reads failing quality controls
                    "filter": {"duplicate": not keep_dups[tech], "vendorFailed": True},
                }
            )
    # normalized calls of every technology (PASS and non-PASS records), same order
    for r in rows:
        tech = r["technology"]
        if (tech, "calls") not in tracks:
            sys.exit(f"ERROR: sample {args.sample}: no normalized VCF of technology {tech}")
        config.append(
            {
                "name": f"{tech.upper()} calls (normalized)",
                "url": tracks[(tech, "calls")],
                "format": "vcf",
                "type": "variant",
                # one row per sample column of the VCF (tumor, normal)
                "displayMode": "EXPANDED",
            }
        )
    # explicit order: as listed (create_report would number them in the same way)
    for i, track in enumerate(config, start=1):
        track["order"] = i
    flags = FLAG_QC_FAIL if any(keep_dups.values()) else FLAG_QC_FAIL + FLAG_DUPLICATE
    return config, flags


def write_tracks(args):
    """Track configuration and excluded SAM flags of one sample."""
    rows = read_samplesheet(args.samplesheet, args.sample)
    tracks = track_files(args.labels, args.files)
    config, flags = track_config(rows, tracks, args)
    with open(args.output, "w") as fh:
        json.dump(config, fh, indent=2)
        fh.write("\n")
    with open(args.flags_output, "w") as fh:
        fh.write(f"{flags}\n")
    log(
        "INFO",
        f"{args.sample}: {len(config)} tracks "
        f"({sum(t['type'] == 'alignment' for t in config)} alignments); "
        f"excluded SAM flags {flags}",
    )


# ----------------------------------------------------------------------------- template


def write_template(args):
    """The template with igv.js inlined instead of the CDN script tag."""
    with open(args.igv_js) as fh:
        # the source-map comment would make the browser look for a file that is not there
        script = "".join(
            line for line in fh if not line.startswith("//# sourceMappingURL=")
        ).rstrip()
    # a closing script tag inside the code would end the inlined script early
    if "</script" in script.lower():
        sys.exit(f"ERROR: {args.igv_js}: contains '</script'; it can't be inlined")
    with open(args.template) as fh:
        lines = fh.readlines()
    # the line(s) that load igv.js: <script src="...igv.min.js"></script>
    marks = [
        i for i, line in enumerate(lines) if "<script src=" in line and IGV_SCRIPT_MARK in line
    ]
    if len(marks) != 1:
        sys.exit(
            f"ERROR: {args.template}: expected one '<script src=...{IGV_SCRIPT_MARK}>' line, "
            f"found {len(marks)}"
        )
    lines[marks[0]] = f'<script type="text/javascript">\n{script}\n</script>\n'
    with open(args.output, "w") as fh:
        fh.writelines(lines)
    log("INFO", f"{args.output}: igv.js inlined ({len(script) / 1e6:.1f} MB)")


# ----------------------------------------------------------------------------- index

INDEX_HTML = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>variantcomp: IGV views</title>
<style>
body {{ font-family: sans-serif; margin: 32px; max-width: 960px; }}
table {{ border-collapse: collapse; }}
th, td {{ border-bottom: 1px solid #ddd; padding: 4px 12px; text-align: left; }}
th {{ background: #f2f2f2; }}
td.n {{ text-align: right; }}
</style></head>
<body>
<h1>IGV views of the compared variants</h1>
<p>Every compared variant (TP, FN and FP of every comparison) and every required known
variant, with the reads of the tumor and normal alignments of every technology and the
normalized calls of every technology. One page per sample and part of at most
{max_sites} sites (all read data are embedded in the pages). Click a row of a page to see
the site; the address of the page then ends with <code>#&lt;variant_id&gt;</code> and can be
shared. Duplicate-flagged reads are hidden (except in amplicon libraries).</p>
<p>See also: <a href="{multiqc}">MultiQC report</a> &middot;
<a href="{tables}">variant tables</a> (column <em>igv</em> links to the page of a variant)</p>
{body}
</body></html>
"""


def write_index(args):
    """Overview page: every sample and part with its number of sites."""
    with open(args.links, newline="") as fh:
        links = list(csv.DictReader(fh, delimiter="\t"))
    # parts in the order of the links table (samples, then parts)
    parts = {}
    for r in links:
        key = (r["sample"], r["part"], r["igv_file"])
        parts[key] = parts.get(key, 0) + 1
    if parts:
        rows = "\n".join(
            f"<tr><td>{html.escape(s)}</td><td>{html.escape(p)}</td><td class='n'>{n}</td>"
            f"<td><a href='{html.escape(f)}'>{html.escape(f)}</a></td></tr>"
            for (s, p, f), n in parts.items()
        )
        body = (
            "<table><thead><tr><th>sample</th><th>part</th><th>sites</th><th>page</th></tr>"
            f"</thead><tbody>\n{rows}\n</tbody></table>"
        )
    else:
        body = "<p>No variant was compared, so there are no IGV views.</p>"
    with open(args.output, "w") as fh:
        fh.write(
            INDEX_HTML.format(
                max_sites=args.max_sites,
                multiqc=html.escape(args.multiqc_link),
                tables=html.escape(args.tables_link),
                body=body,
            )
        )
    log("INFO", f"{args.output}: {len(parts)} pages, {len(links)} sites")


# ----------------------------------------------------------------------------- main


def main():
    """Parse the command line and run the subcommand."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    s = sub.add_parser("sites", help="sites tables per sample and part, links table")
    s.add_argument("--variants", required=True, help="variants.annotated.tsv")
    s.add_argument("--required", help="required_variants.tsv (optional)")
    s.add_argument("--fai", required=True, help="reference .fai")
    s.add_argument("--max-sites", type=int, default=100, help="--igv_max_sites")
    s.add_argument("--flanking", type=int, default=200, help="--igv_flanking (checked)")
    s.add_argument("--window", type=int, help="--igv_window (checked against --flanking)")
    s.add_argument("--outdir", required=True, help="output folder of the sites tables")
    s.add_argument("--links", required=True, help="output: links table")
    t = sub.add_parser("tracks", help="track configuration of one sample")
    t.add_argument("--sample", required=True)
    t.add_argument("--samplesheet", required=True, help="samplesheet.validated.csv")
    t.add_argument("--labels", nargs="+", required=True, help="<technology>:<role> per file")
    t.add_argument("--files", nargs="+", required=True, help="track files (staged names)")
    t.add_argument("--height", type=int, default=250, help="height of an alignment track (px)")
    t.add_argument("--sampling-depth", type=int, default=500, help="--igv_sampling_depth (igv.js)")
    t.add_argument(
        "--amplicon-keep-duplicates",
        type=lambda v: v.lower() == "true",
        default=True,
        help="--amplicon_keep_duplicates (true/false)",
    )
    t.add_argument("--output", required=True, help="output: track configuration (JSON)")
    t.add_argument("--flags-output", required=True, help="output: excluded SAM flags")
    m = sub.add_parser("template", help="report template with igv.js inlined")
    m.add_argument("--template", required=True)
    m.add_argument("--igv-js", required=True, help="igv.min.js (--igv_js)")
    m.add_argument("--output", required=True)
    i = sub.add_parser("index", help="overview page of the IGV views")
    i.add_argument("--links", required=True, help="igv_links.tsv")
    i.add_argument("--max-sites", type=int, default=100, help="--igv_max_sites (shown)")
    i.add_argument("--multiqc-link", default="../../multiqc/multiqc_report.html")
    i.add_argument("--tables-link", default="../tables/index.html")
    i.add_argument("--output", required=True)
    args = parser.parse_args()
    {
        "sites": write_sites,
        "tracks": write_tracks,
        "template": write_template,
        "index": write_index,
    }[args.command](args)


if __name__ == "__main__":
    main()
