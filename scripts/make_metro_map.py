#!/usr/bin/env python3
"""Draw the nf-core-style metro map of the pipeline (docs/images/).

The style follows the nf-core metro-map components (CC0; designed by James A.
Fellows Yates, updated by Maxime Garcia and Louis Le Nezet; see
https://nf-co.re/docs/community/brand/workflow-schematics): 10 px lines in the
nf-core line colours, round white stations with a black border, file icons at
the line starts, 45-degree bends.

Every input type has its own horizontal lane. Lines are drawn through explicit
waypoints (LINES), so that bends never cross station labels; a bend between two
lanes is a 45-degree diagonal (dx == dy).

Usage: scripts/make_metro_map.py --outdir docs/images
       (writes variantcomp_metro_map.svg and, with cairosvg installed, .png)

Edit STATIONS and LINES below when the workflow changes.
"""

import argparse
import sys
from pathlib import Path
from xml.sax.saxutils import escape

# nf-core metro-map line colours
GREEN, ORANGE, BLUE, YELLOW = "#24B064", "#F18046", "#3d95fd", "#ECDC86"
BROWN, PINK = "#3F2B29", "#E879E3"
FONT = "Maven Pro, Helvetica, Arial, sans-serif"
WIDTH, HEIGHT = 1960, 1030

# y coordinate of each lane (one lane per input type); 140 px apart
LANE = {"sheet": 170, "bed": 310, "vcf": 450, "bam": 590, "annot": 730, "known": 870}
# x coordinate of the file icons (line starts)
ICON_X = 70

# station id -> (x, y, label lines, label position 'above' | 'below'[, label x shift])
STATIONS = {
    # samplesheet and reference
    "sheets": (260, LANE["sheet"], ["check sheets", "check_sheets.py"], "above"),
    "index": (400, LANE["sheet"], ["index VCF / BAM", "bcftools, samtools"], "above"),
    "faidx": (780, LANE["sheet"], ["index reference", "samtools faidx"], "above"),
    "genome": (920, LANE["sheet"], ["genome BED", "gawk"], "above"),
    # target BEDs
    "bed_check": (260, LANE["bed"], ["reference check", "check_reference.py"], "above"),
    "bed_sort": (400, LANE["bed"], ["sort", "bedtools sort"], "below"),
    "bed_inter": (540, LANE["bed"], ["intersect all", "bedtools multiinter"], "above"),
    "bed_merge": (680, LANE["bed"], ["merge", "gawk, bedtools merge"], "below"),
    "bed_opt": (820, LANE["bed"], ["regions / exclude", "bedtools intersect"], "above"),
    # somatic VCFs (main line)
    "vcf_check": (260, LANE["vcf"], ["reference check", "bcftools view -h, query", "check_reference.py"], "below"),
    "select": (400, LANE["vcf"], ["select samples", "bcftools view"], "above"),
    "rename": (520, LANE["vcf"], ["rename contigs", "bcftools annotate"], "below"),
    "clean": (640, LANE["vcf"], ["clean records", "bcftools view"], "above"),
    "reheader": (760, LANE["vcf"], ["contig header", "bcftools reheader"], "below"),
    "norm": (880, LANE["vcf"], ["normalize", "bcftools norm"], "above"),
    "ref_alleles": (1000, LANE["vcf"], ["REF check", "check_reference.py"], "below"),
    "harmonize": (1110, LANE["vcf"], ["harmonize", "harmonize_vcf.py"], "above"),
    "filter": (1240, LANE["vcf"], ["filter + regions", "bcftools view"], "below", 52),
    "vcfeval": (1380, LANE["vcf"], ["match", "rtg vcfeval", "bcftools isec"], "above"),
    "metrics": (1520, LANE["vcf"], ["metrics", "benchmark_metrics.py", "aggregate_metrics.py"], "above"),
    # alignments (optional)
    "idxstats": (260, LANE["bam"], ["reference check", "samtools idxstats", "check_reference.py"], "below"),
    "mosdepth": (400, LANE["bam"], ["coverage", "mosdepth"], "below"),
    "callable": (540, LANE["bam"], ["callable bases", "mosdepth_callable.py"], "below"),
    "call_sort": (680, LANE["bam"], ["sort", "bedtools sort"], "below"),
    "call_inter": (820, LANE["bam"], ["callable in all", "multiinter, intersect"], "below"),
    "mpileup": (1380, LANE["bam"], ["read evidence", "bcftools mpileup", "bcftools norm"], "above"),
    "classify": (1520, LANE["bam"], ["discordance", "classify_discordance.py"], "below"),
    # genes and stratifications (optional)
    "annot_check": (260, LANE["annot"], ["reference check", "check_reference.py"], "below"),
    "annot_sort": (400, LANE["annot"], ["sort", "bedtools sort"], "below"),
    "annot_map": (1380, LANE["annot"], ["genes, strata", "bedtools map", "bedtools intersect -c"], "below"),
    "annot_merge": (1520, LANE["annot"], ["annotate", "annotate_variants.py"], "below"),
    # known variants (optional)
    "known_prep": (260, LANE["known"], ["prepare", "known_variants.py"], "below"),
    "known_norm": (400, LANE["known"], ["normalize", "bcftools norm"], "below"),
    "known_check": (540, LANE["known"], ["reference check", "check_reference.py"], "below"),
    "known_match": (1380, LANE["known"], ["match", "bcftools isec", "rtg vcfeval"], "below"),
    "known_sum": (1520, LANE["known"], ["checklist", "known_summary.py"], "below"),
    # outputs
    "report": (1730, LANE["bam"], ["report", "variantcomp_report.py"], "below"),
    "multiqc": (1880, LANE["bam"], ["report", "MultiQC"], "below"),
}

# (colour, waypoints, start icon label or None); a waypoint is a station id or (x, y)
LINES = [
    (BROWN, ["sheets", "index"], "csv"),
    (BROWN, [(620, LANE["sheet"]), "faidx", "genome"], "fasta"),
    (BLUE, ["bed_check", "bed_sort", "bed_inter", "bed_merge", "bed_opt", (1100, LANE["bed"]), "filter"], "bed"),
    (GREEN, ["vcf_check", "select", "rename", "clean", "reheader", "norm", "ref_alleles", "harmonize",
             "filter", "vcfeval", "metrics", (1540, LANE["vcf"]), (1680, LANE["bam"]), "report", "multiqc"],
     "vcf"),  # fmt: skip
    (ORANGE, ["idxstats", "mosdepth", "callable", "call_sort", "call_inter", (1100, LANE["bam"]), "filter"], "bam"),
    (ORANGE, ["call_inter", "mpileup", "classify", "report"], None),
    (PINK, ["annot_check", "annot_sort", "annot_map", "annot_merge", (1540, LANE["annot"]), (1680, LANE["bam"]),
            "report"], "bed"),  # fmt: skip
    (YELLOW, ["known_prep", "known_norm", "known_check", "known_match", "known_sum", (1540, LANE["known"]),
              (1820, LANE["bam"]), "multiqc"], "csv"),  # fmt: skip
]

LEGEND = [
    (GREEN, "somatic VCFs (all technologies)"),
    (BLUE, "target BEDs"),
    (ORANGE, "alignments (optional)"),
    (PINK, "genes, stratifications (optional)"),
    (YELLOW, "known variants (optional)"),
    (BROWN, "samplesheet, reference"),
]


def point(waypoint):
    """Coordinates of a waypoint: a station id or an (x, y) tuple."""
    if isinstance(waypoint, str):
        # station ids must exist; a typo would silently draw a wrong map
        if waypoint not in STATIONS:
            sys.exit(f"ERROR: unknown station '{waypoint}' in LINES")
        return STATIONS[waypoint][:2]
    return waypoint


def check_bends(points):
    """Every segment must be horizontal, vertical or a 45-degree diagonal."""
    for (x0, y0), (x1, y1) in zip(points, points[1:], strict=False):
        if y0 != y1 and x0 != x1 and abs(x1 - x0) != abs(y1 - y0):
            sys.exit(f"ERROR: segment ({x0},{y0}) -> ({x1},{y1}) is not a 45-degree bend")


def file_icon(x, y, label):
    """nf-core style file icon (document with folded corner and label box), centred on y."""
    w, h, fold = 52, 64, 16
    top = y - h / 2
    return (
        f'<g stroke="#000" stroke-width="5" fill="#fff" stroke-linejoin="round">'
        f'<path d="M {x} {top} L {x + w - fold} {top} L {x + w} {top + fold} L {x + w} {top + h} '
        f'L {x} {top + h} Z"/>'
        f'<path d="M {x + w - fold} {top} L {x + w - fold} {top + fold} L {x + w} {top + fold}" fill="#000"/>'
        f"</g>"
        f'<rect x="{x - 10}" y="{top + 21}" width="{w + 20}" height="26" rx="5" fill="#000"/>'
        f'<text x="{x + w / 2}" y="{top + 40}" font-family="{FONT}" font-size="18" font-weight="bold" '
        f'fill="#fff" text-anchor="middle">{escape(label)}</text>'
    )


def station(x, y, labels, position, shift=0):
    """Station circle and its labels (first line bold: the step, then the tools).

    shift moves the labels to the right, away from a line that passes the station.
    """
    out = f'<circle cx="{x}" cy="{y}" r="11" fill="#fff" stroke="#000" stroke-width="6"/>'
    for i, text in enumerate(labels):
        # labels above the station grow upwards, labels below grow downwards
        if position == "above":
            offset = -24 - 17 * (len(labels) - 1 - i)
        else:
            offset = 36 + 17 * i
        weight = "bold" if i == 0 else "normal"
        out += (f'<text x="{x + shift}" y="{y + offset}" font-family="{FONT}" font-size="14" font-weight="{weight}" '
                f'text-anchor="middle" fill="#000">{escape(text)}</text>')  # fmt: skip
    return out


def build_svg():
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">',
        f'<rect width="{WIDTH}" height="{HEIGHT}" fill="#fff"/>',
        f'<text x="40" y="58" font-family="{FONT}" font-size="32" font-weight="bold">'
        "core-unit-bioinformatics/workflow-nxf-somatic-variantcomp</text>",
    ]
    # lines first, so that stations and icons are drawn on top
    icons = []
    for colour, waypoints, icon in LINES:
        points = [point(w) for w in waypoints]
        if icon:
            # the line starts at a file icon: at the first waypoint if that is a coordinate,
            # else at the left margin of the lane
            if not isinstance(waypoints[0], str):
                icons.append(file_icon(points[0][0] - 26, points[0][1], icon))
            else:
                points.insert(0, (ICON_X + 26, points[0][1]))
                icons.append(file_icon(ICON_X, points[0][1], icon))
        check_bends(points)
        d = " ".join(f"{'M' if i == 0 else 'L'} {x} {y}" for i, (x, y) in enumerate(points))
        parts.append(f'<path d="{d}" fill="none" stroke="{colour}" stroke-width="10" '
                     'stroke-linecap="round" stroke-linejoin="round"/>')  # fmt: skip
    parts.extend(icons)
    for x, y, *label in STATIONS.values():
        parts.append(station(x, y, *label))
    # legend: one row below the map
    x = 40
    for colour, text in LEGEND:
        parts.append(f'<line x1="{x}" y1="975" x2="{x + 44}" y2="975" stroke="{colour}" stroke-width="10" '
                     'stroke-linecap="round"/>')  # fmt: skip
        parts.append(f'<text x="{x + 58}" y="981" font-family="{FONT}" font-size="16">{escape(text)}</text>')
        x += 70 + 8.5 * len(text)
    parts.append(f'<text x="40" y="{HEIGHT - 14}" font-family="{FONT}" font-size="12" fill="#555">'
                 "Metro-map style: nf-core components (CC0), J. A. Fellows Yates, M. Garcia, L. Le Nezet</text>")  # fmt: skip
    parts.append("</svg>")
    return "\n".join(parts)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--outdir", default="docs/images", type=Path, help="output directory")
    args = parser.parse_args()
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
