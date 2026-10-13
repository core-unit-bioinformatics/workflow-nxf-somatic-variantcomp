#!/usr/bin/env python3
"""Draw the nf-core-style metro map of the pipeline and its step table.

Layout: the pipeline runs from top to bottom. Every subworkflow is a box; every row
in a box is one step that a sample goes through, with its number, what happens and
the tools. Each coloured vertical line is one use case or option (see LINES): a
line has a station in every step that runs in that use case. A station spanning
several lines is a step shared by these use cases.

The style follows the nf-core metro-map components (CC0; James A. Fellows Yates,
Maxime U. Garcia, Louis Le Nezet; https://nf-co.re/docs/community/brand/workflow-schematics):
nf-core line colours, white stations with black borders, file icons at the line
starts.

Usage:
  scripts/make_metro_map.py --outdir docs/images     SVG (+ PNG with cairosvg)
  scripts/make_metro_map.py --markdown               step table for the README

Edit BANDS when the workflow changes; the map and the README table come from it.
"""

import argparse
import sys
from pathlib import Path
from xml.sax.saxutils import escape

FONT = "Maven Pro, Helvetica, Arial, sans-serif"
MONO = "DejaVu Sans Mono, Menlo, Consolas, monospace"

# use cases / options: key -> (colour, start icon label, legend text); the order is
# the left-to-right order of the vertical lines
LINES = {
    "regions": (
        "#3d95fd",
        "bed",
        "Targeted technologies: target BEDs (+ --regions_bed / --exclude_bed)",
    ),
    "core": (
        "#24B064",
        "vcf",
        "Always: compare the VCF of every technology with the truth technology",
    ),
    "bam": ("#F18046", "bam", "With alignments (BAM/CRAM): callable regions and read evidence"),
    "required": (
        "#ECDC86",
        "csv",
        "--known_variants, mode required: clinical checklist (PASS/FAIL)",
    ),
    "truthset": (
        "#FA6863",
        "vcf",
        "--known_variants, mode benchmark: external truth set (e.g. SEQC2)",
    ),
    "annot": ("#E879E3", "bed", "With --gene_bed / --stratification: genes and genomic context"),
}
ALL = list(LINES)

# subworkflow boxes: (subworkflow name(s), purpose, steps)
# step: (number, what happens, tools, lines that run this step)
BANDS = [
    (
        "PIPELINE_INITIALISATION, CHECK_SAMPLESHEET, PREPARE_ALIGNMENTS",
        "read and check the samplesheet",
        [
            (
                "1.1",
                "Read the samplesheet, check every column",
                "nf-schema",
                ["regions", "core", "bam"],
            ),
            (
                "1.2",
                "Check rows: truth technology per sample, one BED per technology",
                "check_sheets.py",
                ["regions", "core", "bam"],
            ),
            ("1.3", "Index alignments without an index", "samtools index", ["bam"]),
        ],
    ),
    (
        "PREPARE_REFERENCE",
        "index the reference genome",
        [
            ("2.1", "Unpack and index the reference FASTA", "gunzip, samtools faidx", ["core"]),
            ("2.2", "Genome BED: one interval per reference contig", "gawk", ["core"]),
        ],
    ),
    (
        "PREPARE_REGIONS",
        "build the common comparison region",
        [
            (
                "3.1",
                "Check target BEDs against the reference, rename contigs",
                "check_reference.py",
                ["regions"],
            ),
            ("3.2", "Sort the BEDs in reference order", "bedtools sort", ["regions"]),
            (
                "3.3",
                "Keep bases targeted by every technology",
                "bedtools multiinter, gawk, merge",
                ["regions"],
            ),
            (
                "3.4",
                "Restrict to --regions_bed, remove --exclude_bed (optional)",
                "bedtools intersect, subtract",
                ["regions"],
            ),
            (
                "3.5",
                "Common region (genome if no target BEDs) and its size",
                "bedtools merge, region_summary.py",
                ["regions", "core"],
            ),
        ],
    ),
    (
        "COVERAGE",
        "find the bases callable in all technologies of a sample",
        [
            (
                "4.1",
                "Check alignment contigs against the reference",
                "samtools idxstats, check_reference.py",
                ["bam"],
            ),
            ("4.2", "Depth per base and per common-region interval", "mosdepth", ["bam"]),
            (
                "4.3",
                "Callable bases per alignment (depth >= --min_coverage)",
                "mosdepth_callable.py, bedtools sort",
                ["bam"],
            ),
            (
                "4.4",
                "Bases callable in all alignments of a sample (in common region)",
                "bedtools multiinter, gawk, intersect",
                ["bam"],
            ),
        ],
    ),
    (
        "PREPROCESS_VCF",
        "normalize every VCF in the same way",
        [
            (
                "5.1",
                "Check VCF contigs, contig lengths and sample names",
                "bcftools view -h, query, check_reference.py",
                ["core"],
            ),
            (
                "5.2",
                "Keep tumor (+ normal) columns, recompress with BGZF",
                "bcftools view",
                ["core"],
            ),
            (
                "5.3",
                "Rename contigs to the reference (chr prefix, MT/chrM)",
                "bcftools annotate",
                ["core"],
            ),
            ("5.4", "Remove records on other contigs and gVCF blocks", "bcftools view", ["core"]),
            ("5.5", "Rewrite the contig header from the reference", "bcftools reheader", ["core"]),
            (
                "5.6",
                "Split multiallelics, left-align, split MNVs, deduplicate",
                "bcftools norm",
                ["core"],
            ),
            (
                "5.7",
                "Stop if too many REF alleles differ (other genome build)",
                "check_reference.py",
                ["core"],
            ),
            (
                "5.8",
                "Same VAF, depth and alt reads for every caller (from AD)",
                "harmonize_vcf.py",
                ["core"],
            ),
            (
                "5.9",
                "Uniform filters; keep the common (or callable) region",
                "bcftools view",
                ["core", "bam"],
            ),
        ],
    ),
    (
        "BENCHMARK",
        "compare every technology with the truth technology",
        [
            ("6.1", "Match calls haplotype-aware: TP, FN, FP", "rtg vcfeval", ["core"]),
            ("6.2", "Exact-allele and same-position matches", "bcftools isec", ["core"]),
            (
                "6.3",
                "Recall, precision, F1 per sample (overall, type, VAF)",
                "benchmark_metrics.py",
                ["core"],
            ),
            (
                "6.4",
                "Combine samples: summed (micro) and mean (macro)",
                "aggregate_metrics.py",
                ["core"],
            ),
        ],
    ),
    (
        "DISCORDANCE",
        "explain why variants differ",
        [
            ("7.1", "Positions of all compared variants of a sample", "bcftools query", ["bam"]),
            (
                "7.2",
                "Read counts of every alignment at these positions",
                "bcftools mpileup + call, bcftools norm",
                ["bam"],
            ),
            (
                "7.3",
                "Reason for every FN/FP; detection power of the query",
                "classify_discordance.py",
                ["core", "bam"],
            ),
            (
                "7.4",
                "Check the stratification sheet, the gene and stratification BEDs",
                "check_sheets.py, check_reference.py, bedtools sort",
                ["annot"],
            ),
            (
                "7.5",
                "Add genes and genomic context to every variant",
                "bedtools map, intersect -c, annotate_variants.py",
                ["annot"],
            ),
        ],
    ),
    (
        "PREPARE_KNOWN_VARIANTS, KNOWN_VARIANTS",
        "check known variants",
        [
            (
                "8.1",
                "Check the known-variant sheet; read the lists, check positions",
                "check_sheets.py, known_variants.py prepare",
                ["required", "truthset"],
            ),
            (
                "8.2",
                "Normalize like the calls; stop on any REF mismatch",
                "bcftools norm, check_reference.py",
                ["required", "truthset"],
            ),
            (
                "8.3",
                "PASS/FAIL of each variant per sample and technology",
                "known_variants.py check",
                ["bam", "required"],
            ),
            ("8.4", "Checklist summary over all samples", "known_summary.py", ["required"]),
            (
                "8.5",
                "Recall/precision of every technology vs the truth set",
                "rtg vcfeval, bcftools isec, benchmark_metrics.py",
                ["truthset"],
            ),
        ],
    ),
    (
        "REPORT",
        "summarize",
        [
            (
                "9.1",
                "Enrichment statistics, gene table and plots",
                "variantcomp_report.py",
                ["core", "annot"],
            ),
            (
                "9.2",
                "Interactive report of all results",
                "MultiQC",
                [k for k in ALL if k != "regions"],
            ),
        ],
    ),
]

# longer explanation of every step for the README table (the map shows the short text)
DETAILS = {
    "1.1": "The samplesheet is read as in the nf-core template: nf-schema checks every column (names, library types, file extensions). Relative paths are resolved against the samplesheet's directory; VCF and alignment indices are looked up next to the files.",
    "1.2": "Rules across rows, before anything else runs: every sample has the truth technology and another technology, no duplicate rows, one target BED and library type per technology, existing files. All problems are reported at once with line numbers.",
    "1.3": "Alignments without a `.bai`/`.csi`/`.crai` next to them are indexed once; they are needed for coverage and read evidence.",
    "2.1": "A gzipped FASTA is unpacked and indexed; the `.fai` (contig names and lengths) is the reference every input is checked against.",
    "2.2": "One interval per contig, used as comparison region when no technology has a target BED (WGS only) and to drop records on other contigs.",
    "3.1": "Contig names that differ only by the `chr` prefix or `MT`/`chrM` are renamed, intervals on unknown contigs are dropped and reported; intervals beyond a contig end stop the run (other genome build).",
    "3.2": "All BEDs are sorted in the contig order of the reference, which the following bedtools steps need.",
    "3.3": "Only bases covered by the target BED of every targeted technology are kept (WGS counts as genome-wide), so every technology could have called a variant there.",
    "3.4": "Optional: the region is restricted to `--regions_bed` and `--exclude_bed` is removed (e.g. difficult regions).",
    "3.5": "The final common comparison region (`regions/common_targets.bed`), merged; its size and the share of each technology's targets are summarized.",
    "4.1": "Contig names and lengths of every alignment are compared with the reference (same rules as for BEDs); the common region is translated to the alignment's contig names.",
    "4.2": "mosdepth classifies every base as below or above `--min_coverage` and computes the mean depth per common-region interval (duplicates are kept for amplicon libraries).",
    "4.3": "The callable bases of each alignment (base mode; region mode uses the mean depth per interval), renamed to reference contigs and sorted.",
    "4.4": "Per sample: bases callable in the tumor (and normal) alignment of every technology, inside the common region. Variants are only compared there, so missing coverage in one technology does not count as a miss.",
    "5.1": "Contig names and lengths of the VCF header and records are checked against the reference, and the tumor/normal sample columns must exist. Safe fixes are prepared, a different build stops the run.",
    "5.2": "Only the tumor (and normal) columns are kept, in this order, and the VCF is written as bgzipped VCF with index (inputs may be plain gzip).",
    "5.3": "Contigs are renamed to the reference names (`1` -> `chr1`, `MT` -> `chrM`).",
    "5.4": "Records on contigs absent from the reference and gVCF reference blocks (only ALT `<NON_REF>`/`<*>`) are removed.",
    "5.5": "Contig header lines are rewritten from the reference index (fixes wrong lengths, e.g. chrM, and missing contig lines).",
    "5.6": "Multiallelic records are split, indels left-aligned, MNVs split into SNVs, records whose REF differs from the FASTA removed and counted, duplicates removed. Same representation for every caller.",
    "5.7": "The fraction of REF mismatches from step 5.6 is checked: above `--max_ref_mismatch` the VCF is on another reference and the run stops.",
    "5.8": "VAF, depth and alt reads are computed the same way for every caller (from FORMAT/AD, else AF + DP; documented in docs/usage.md) and stored as `VC_*` INFO tags; GT is set to tumor 0/1, normal 0/0. Every tumor record needs a VAF.",
    "5.9": "The same filters for every technology (`--filter_pass_only`, `--min_vaf`, `--min_depth`, `--min_alt_reads`) and restriction to the comparison region: the common region, or with alignments the sample's callable region. The unfiltered VCF is kept for the discordance step.",
    "6.1": "Each technology's calls are matched with the truth technology's calls of the same sample, haplotype-aware, so different representations of one variant still match (TP, FN = missed, FP = only in the query).",
    "6.2": "Record-by-record comparison by CHROM/POS/REF/ALT and by position only, to see how much discordance is representation or another allele.",
    "6.3": "Recall, precision and F1 per sample and technology: overall, SNV vs. indel, by truth VAF (detection limit) and query VAF (low-VAF false positives), exact and position matches.",
    "6.4": "Over all samples: from the summed counts (micro) and as mean +- sd of the per-sample values (macro). Main result: `benchmark/metrics.aggregated.tsv`.",
    "7.1": "All positions of compared variants of a sample are collected, so every alignment is looked at the same sites.",
    "7.2": "Allele counts of every tumor and normal alignment at these positions (MQ >= 20, BQ >= 13), normalized like the calls.",
    "7.3": "Each FN/FP gets a reason from the other technology's view: called but filtered, germline (alt reads in the normal), low coverage, alt reads below the caller's threshold, no alt reads, or absent (no alignment). For truth variants, the binomial detection power of the query gives the power-adjusted recall.",
    "7.4": "The `--stratification` sheet is checked (columns, unique names, existing BEDs); gene and stratification BEDs (repeats, homopolymers, ...) are checked against the reference and sorted.",
    "7.5": "Every compared variant gets its gene(s) and a yes/no per stratification; result: `discordance/variants.annotated.tsv`.",
    "8.1": "The `--known_variants` sheet is checked (columns, unique names, modes, samples of the samplesheet, existing files). Known variants from CSV, TSV or VCF lists are mapped to reference contigs, positions and alleles checked, and written as a sorted VCF.",
    "8.2": "The list is normalized like the calls (split, left-align, atomize); any REF allele that differs from the FASTA stops the run.",
    "8.3": "Required mode: each listed variant must be found in every technology of the sample: PASS, or FAIL as filtered, outside the comparison region or not called (with read evidence if alignments are given).",
    "8.4": "Required mode: one table over all lists, samples and technologies, a warning per failed variant, and optionally a stop (`--fail_on_missing_required`).",
    "8.5": "Benchmark mode: every technology is compared with the external truth set inside its confident regions, with the same metrics as the main benchmark.",
    "9.1": "Are discordant variants enriched for indels, low VAF, stratifications or genes compared with concordant ones (Fisher's exact test, Benjamini-Hochberg)? Plus gene table, class counts and plots.",
    "9.2": "One report with metrics, discordance classes, recall by VAF, coverage (mosdepth), required variants and software versions.",
}

# geometry
TRACK_X0, TRACK_DX = 90, 62  # x of the first vertical line, distance between lines
LABEL_X = 470  # step number and text
TOOLS_X = 1110  # tools column
WIDTH = 1620
ROW = 38  # height of a step row
BAND_HEAD = 44  # height of a box header
BAND_GAP = 14  # space between boxes
TOP = 330  # y of the first box


def track_x(line):
    """x coordinate of a line's vertical track."""
    return TRACK_X0 + TRACK_DX * list(LINES).index(line)


def layout():
    """Positions: boxes [(y0, y1, band)], stations [(y, step)], total height."""
    boxes, stations = [], []
    y = TOP
    for band in BANDS:
        y0 = y
        y += BAND_HEAD
        for step in band[2]:
            stations.append((y + ROW / 2, step))
            y += ROW
        y += 10
        boxes.append((y0, y, band))
        y += BAND_GAP
    return boxes, stations, y + 40


def check_data():
    """Every step uses known lines; every line has at least one station."""
    used = set()
    for _name, _purpose, steps in BANDS:
        for number, _text, _tools, lines in steps:
            unknown = set(lines) - set(LINES)
            if unknown:
                sys.exit(f"ERROR: step {number}: unknown line(s) {', '.join(sorted(unknown))}")
            used.update(lines)
    unused = set(LINES) - used
    if unused:
        sys.exit(f"ERROR: line(s) without station: {', '.join(sorted(unused))}")


def file_icon(cx, y, label, colour):
    """nf-core style file icon centred on cx, top at y, with a coloured label box."""
    w, h, fold = 40, 50, 12
    x = cx - w / 2
    return (
        f'<g stroke="#000" stroke-width="4" fill="#fff" stroke-linejoin="round">'
        f'<path d="M {x} {y} L {x + w - fold} {y} L {x + w} {y + fold} L {x + w} {y + h} L {x} {y + h} Z"/>'
        f'<path d="M {x + w - fold} {y} L {x + w - fold} {y + fold} L {x + w} {y + fold}" fill="#000"/></g>'
        f'<rect x="{x - 6}" y="{y + 16}" width="{w + 12}" height="22" rx="4" fill="{colour}" stroke="#000" stroke-width="2"/>'
        f'<text x="{cx}" y="{y + 32}" font-family="{FONT}" font-size="15" font-weight="bold" '
        f'text-anchor="middle" fill="#000">{escape(label)}</text>'
    )


def build_svg():
    boxes, stations, height = layout()
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{height}" '
        f'viewBox="0 0 {WIDTH} {height}">',
        f'<rect width="{WIDTH}" height="{height}" fill="#fff"/>',
    ]
    # title
    out.append(
        f'<text x="40" y="56" font-family="{FONT}" font-size="30" font-weight="bold">'
        "core-unit-bioinformatics/workflow-nxf-somatic-variantcomp</text>"
    )
    out.append(
        f'<text x="40" y="88" font-family="{FONT}" font-size="17" fill="#444">'
        "Each box is a subworkflow, each row a step every sample goes through (top to bottom). "
        "Each line is a use case: it stops at the steps that run in it.</text>"
    )
    # legend: one row per line, right of the file icons
    for i, (colour, _icon, text) in enumerate(LINES.values()):
        y = 136 + 30 * i
        out.append(
            f'<line x1="{LABEL_X}" y1="{y}" x2="{LABEL_X + 50}" y2="{y}" stroke="{colour}" '
            'stroke-width="10" stroke-linecap="round"/>'
        )
        out.append(
            f'<text x="{LABEL_X + 66}" y="{y + 6}" font-family="{FONT}" font-size="17">{escape(text)}</text>'
        )

    # subworkflow boxes with their header
    for y0, y1, (name, purpose, _steps) in boxes:
        out.append(
            f'<rect x="30" y="{y0}" width="{WIDTH - 60}" height="{y1 - y0}" rx="12" '
            'fill="#f7f7f7" stroke="#9a9a9a" stroke-width="2" stroke-dasharray="7 5"/>'
        )
        out.append(
            f'<text x="{LABEL_X}" y="{y0 + 29}" font-family="{MONO}" font-size="17" font-weight="bold" '
            f'fill="#333">{escape(name)}</text>'
        )
        out.append(
            f'<text x="{WIDTH - 50}" y="{y0 + 29}" font-family="{FONT}" font-size="16" '
            f'font-style="italic" fill="#555" text-anchor="end">{escape(purpose)}</text>'
        )

    # vertical lines: from the file icon to the last station of the line
    last = {line: max(y for y, step in stations if line in step[3]) for line in LINES}
    icon_y = 210
    for line, (colour, icon, _text) in LINES.items():
        x = track_x(line)
        out.append(
            f'<line x1="{x}" y1="{icon_y + 50}" x2="{x}" y2="{last[line]}" stroke="{colour}" '
            'stroke-width="10" stroke-linecap="round"/>'
        )
        out.append(file_icon(x, icon_y, icon, colour))

    # stations and step labels
    for y, (number, text, tools, lines) in stations:
        xs = [track_x(line) for line in lines]
        if len(xs) == 1:
            out.append(
                f'<circle cx="{xs[0]}" cy="{y}" r="11" fill="#fff" stroke="#000" stroke-width="5"/>'
            )
        elif max(xs) - min(xs) == TRACK_DX * (len(xs) - 1):
            # a step shared by use cases on neighbouring lines: one station across them
            x0, x1 = min(xs) - 14, max(xs) + 14
            out.append(
                f'<rect x="{x0}" y="{y - 13}" width="{x1 - x0}" height="26" rx="13" fill="#fff" '
                'stroke="#000" stroke-width="5"/>'
            )
        else:
            # lines in between don't stop here: one station per line, joined by a bar
            out.append(
                f'<line x1="{min(xs)}" y1="{y}" x2="{max(xs)}" y2="{y}" stroke="#000" stroke-width="7"/>'
            )
            for x in xs:
                out.append(
                    f'<circle cx="{x}" cy="{y}" r="11" fill="#fff" stroke="#000" stroke-width="5"/>'
                )
        # thin guide from the right-most station to the label
        out.append(
            f'<line x1="{max(xs) + 16}" y1="{y}" x2="{LABEL_X - 10}" y2="{y}" stroke="#d0d0d0" '
            'stroke-width="1.5" stroke-dasharray="2 4"/>'
        )
        out.append(
            f'<text x="{LABEL_X}" y="{y + 6}" font-family="{FONT}" font-size="17">'
            f'<tspan font-weight="bold">{number}</tspan>  {escape(text)}</text>'
        )
        out.append(
            f'<text x="{TOOLS_X}" y="{y + 6}" font-family="{MONO}" font-size="14" fill="#555">'
            f"{escape(tools)}</text>"
        )

    out.append(
        f'<text x="40" y="{height - 16}" font-family="{FONT}" font-size="13" fill="#666">'
        "Metro-map style: nf-core components (CC0), J. A. Fellows Yates, M. U. Garcia, L. Le Nezet. "
        "Step details: README.md</text>"
    )
    out.append("</svg>")
    return "\n".join(out)


def markdown():
    """README table of all steps, generated from the same data as the map."""
    names = {key: text.split(":")[0] for key, (_c, _i, text) in LINES.items()}
    rows = ["| Step | What happens | Tools | Runs for |", "| --- | --- | --- | --- |"]
    for name, purpose, steps in BANDS:
        rows.append(f"| **{name}** | _{purpose}_ | | |")
        for number, text, tools, lines in steps:
            # every step needs a longer explanation
            if number not in DETAILS:
                sys.exit(f"ERROR: step {number} has no entry in DETAILS")
            runs = (
                "all use cases"
                if set(lines) >= set(ALL) - {"regions"}
                else ", ".join(names[x] for x in lines)
            )
            rows.append(f"| {number} | **{text}.** {DETAILS[number]} | `{tools}` | {runs} |")
    return "\n".join(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--outdir", default="docs/images", type=Path, help="output directory")
    parser.add_argument("--markdown", action="store_true", help="print the README step table")
    args = parser.parse_args()
    check_data()
    if args.markdown:
        print(markdown())
        return
    args.outdir.mkdir(parents=True, exist_ok=True)
    svg = args.outdir / "variantcomp_metro_map.svg"
    svg.write_text(build_svg())
    print(f"INFO: wrote {svg}", file=sys.stderr)
    # the PNG is optional (cairosvg is in the dev environment)
    try:
        import cairosvg
    except ImportError:
        print("INFO: cairosvg not installed; only the SVG was written", file=sys.stderr)
        return
    png = svg.with_suffix(".png")
    cairosvg.svg2png(url=str(svg), write_to=str(png), output_width=WIDTH)
    print(f"INFO: wrote {png}", file=sys.stderr)


if __name__ == "__main__":
    main()
