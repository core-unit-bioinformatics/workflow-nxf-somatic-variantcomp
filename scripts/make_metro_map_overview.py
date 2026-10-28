#!/usr/bin/env python3
"""Draw the overview metro map of the pipeline (nf-core/rnaseq style).

The overview complements the detailed map of `make_metro_map.py`: it shows only the
stages and the main tool of each step.

Layout: the VCFs start the pipeline on the left and run through Check and Normalize.
The two comparisons, the key step, sit next to each other in the centre: Benchmark on
top (the line forks into the two matching tools, rtg vcfeval and bcftools isec, and
joins again for the metrics) and Known variants below. On the right, Discordance
explains the benchmark result and Report collects everything.

Lines: one colour. Every other input enters in the stage where it is used: the target
BEDs and the known-variant sheet start their own short line there (solid), files that a
single step reads (reference FASTA, alignments, VEP cache, stratification BEDs) reach
that step with a dashed line.

The style follows the nf-core metro-map components (CC0; James A. Fellows Yates,
Maxime U. Garcia, Louis Le Nezet; https://nf-co.re/docs/community/brand/workflow-schematics)
and the nf-core/rnaseq metro map.

Usage:
  scripts/make_metro_map_overview.py --outdir docs/images   SVG (+ PNG with cairosvg)

Edit LINES, STATIONS, INPUTS and STAGES when the workflow changes.
"""

import argparse
import math
import sys
from pathlib import Path
from xml.sax.saxutils import escape

FONT = "Maven Pro, Helvetica, Arial, sans-serif"
MONO = "DejaVu Sans Mono, Menlo, Consolas, monospace"

WIDTH = 2930
HEIGHT = 1010
LINE_COLOUR = "#24B064"
LINE_WIDTH = 10
RADIUS = 50
DASH_COLOUR = "#3F2B29"
KEY_FILL = "#d6efe0"

# y of the line in: row 1 (Check, Normalize, Benchmark, Discordance), Regions,
# the known-variant sheet, Known variants, Report; y of the benchmark fork
Y1, Y3, YK, Y5, Y7 = 310, 720, 700, 860, 740
Y_RTG, Y_ISEC = 225, 395

# lines as centerline points; a point may carry its own corner radius as third value
LINES = [
    # VCFs: Check and Normalize, up to the fork into the two comparisons
    [(125, Y1), (1420, Y1)],
    # benchmark fork: haplotype-aware (rtg vcfeval) and exact allele / position (isec)
    [(1420, Y1), (1505, Y_RTG), (1685, Y_RTG), (1770, Y1)],
    [(1420, Y1), (1505, Y_ISEC), (1685, Y_ISEC), (1770, Y1)],
    # benchmark metrics, Discordance, down the right edge, Report to the outputs
    [(1770, Y1), (2800, Y1, 60), (2800, Y7, 60), (2130, Y7)],
    # normalized VCFs to Known variants, on to the consistency checks of the Report
    [(1300, Y1, 40), (1370, Y1, 40), (1370, Y5), (2560, Y5), (2560, Y7)],
    # target BEDs: common region, callable bases, into the uniform filters
    [(705, Y3), (1140, Y3), (1140, Y1)],
    # known-variant sheet: normalized like the calls, matched in known_variants.py
    [(1495, YK), (1680, YK), (1680, Y5)],
]

# stations: (label, x, y, label side); "\n" breaks the label into lines
STATIONS = [
    # 1 Check
    ("nf-schema", 200, Y1, "above"),
    ("check_\nsheets.py", 320, Y1, "below"),
    ("check_\nreference.py", 470, Y1, "below"),
    # 2 Normalize
    ("bcftools\nview", 680, Y1, "above"),
    ("bcftools\nannotate", 790, Y1, "below"),
    ("bcftools\nnorm", 915, Y1, "below"),
    ("harmonize_\nvcf.py", 1010, Y1, "above"),
    ("bcftools\nview", 1140, Y1, "above"),
    ("variant_\naccounting.py", 1250, Y1, "below"),
    # 3 Regions
    ("bedtools\nmultiinter", 800, Y3, "below"),
    ("mosdepth", 930, Y3, "above"),
    ("bedtools\nintersect", 1050, Y3, "below"),
    # 4 Benchmark
    ("rtg\nvcfeval", 1595, Y_RTG, "below"),
    ("bcftools\nisec", 1595, Y_ISEC, "above"),
    ("benchmark_\nmetrics.py", 1840, Y1, "above"),
    ("aggregate_\nmetrics.py", 1930, Y1, "below"),
    # 5 Known variants
    ("bcftools\nnorm", 1560, YK, "above"),
    ("known_\nvariants.py", 1680, Y5, "below"),
    ("rtg\nvcfeval", 1815, Y5, "below"),
    ("known_\nsummary.py", 1915, Y5, "above"),
    # 6 Discordance
    ("read_\nevidence.py", 2060, Y1, "below"),
    ("classify_\ndiscordance.py", 2190, Y1, "above"),
    ("Ensembl\nVEP", 2320, Y1, "below"),
    ("msisensor-\npro", 2450, Y1, "below"),
    ("bedtools\nintersect", 2580, Y1, "below"),
    ("annotate_\nvariants.py", 2710, Y1, "above"),
    # 7 Report
    ("variantcomp_\nreport.py", 2690, Y7, "below"),
    ("consistency_\nchecks.py", 2560, Y7, "above"),
    ("report_\nhtml.py", 2430, Y7, "below"),
    ("datavzrd", 2310, Y7, "above"),
    ("MultiQC", 2200, Y7, "below"),
]

# line starts and ends: (icon label, caption, icon centre x, line y, terminus x)
STARTS = [
    ("vcf", "VCFs", 75, Y1, 125),
    ("bed", "target BEDs", 655, Y3, 705),
    ("csv", "known variants", 1445, YK, 1495),
    ("html", "reports", 2070, Y7, 2130),
]

# files read by one step (dashed): (icon label, caption, station x, station y, icon top y)
INPUTS = [
    ("fasta", "reference", 470, Y1, 110),
    ("fasta", "reference", 915, Y1, 110),
    ("bam", "alignments", 930, Y3, 805),
    ("fasta", "reference", 1595, Y_RTG, 92),
    ("fasta", "reference", 1815, Y5, 640),
    ("bam", "alignments", 2060, Y1, 110),
    ("vep", "VEP cache", 2320, Y1, 110),
    ("fasta", "reference", 2450, Y1, 110),
    ("bed", "strata", 2580, Y1, 110),
]

# short notes on the two benchmark branches: (text, x, y)
NOTES = [
    ("haplotype-aware", 1612, Y_RTG - 18),
    ("exact allele, position", 1612, Y_ISEC + 34),
]

# stage boxes: (number, word, hint for the legend, x0, y0, x1, y1, header corner, key step)
STAGES = [
    (1, "Check", "every input against the reference", 30, 90, 590, 500, "bl", False),
    (2, "Normalize", "one chain and filter set for all VCFs", 610, 90, 1310, 500, "bl", False),
    (3, "Regions", "common target region, callable bases", 610, 530, 1310, 970, "tl", False),
    (4, "Benchmark", "recall, precision vs. the truth technology", 1330, 90, 1970, 500, "br", True),
    (
        5,
        "Known variants",
        "required checklist, external truth sets",
        1330,
        530,
        1970,
        970,
        "tr",
        True,
    ),
    (6, "Discordance", "why calls differ; annotation", 1990, 90, 2900, 500, "bl", False),
    (7, "Report", "HTML report, table browser, MultiQC", 1990, 530, 2900, 970, "tl", False),
]


def unit(ax, ay, bx, by):
    """Unit vector from point a to point b."""
    length = math.hypot(bx - ax, by - ay)
    return (bx - ax) / length, (by - ay) / length


def line_path(points):
    """SVG path through the points with rounded bends."""
    xy = [(p[0], p[1]) for p in points]
    d = [f"M {xy[0][0]} {xy[0][1]}"]
    for i in range(1, len(xy) - 1):
        x, y = xy[i]
        u_in = unit(*xy[i - 1], x, y)
        u_out = unit(x, y, *xy[i + 1])
        radius = points[i][2] if len(points[i]) > 2 else RADIUS
        # turn angle and distance of the arc ends from the bend
        angle = math.acos(max(-1, min(1, u_in[0] * u_out[0] + u_in[1] * u_out[1])))
        t = radius * math.tan(angle / 2)
        # clockwise on screen if the cross product is positive
        sweep = 1 if u_in[0] * u_out[1] - u_in[1] * u_out[0] > 0 else 0
        d.append(f"L {x - t * u_in[0]:.1f} {y - t * u_in[1]:.1f}")
        d.append(f"A {radius} {radius} 0 0 {sweep} {x + t * u_out[0]:.1f} {y + t * u_out[1]:.1f}")
    d.append(f"L {xy[-1][0]} {xy[-1][1]}")
    return " ".join(d)


def on_line(x, y):
    """True if (x, y) lies on a horizontal or vertical segment of a line."""
    for points in LINES:
        for (ax, ay, *_), (bx, by, *_) in zip(points[:-1], points[1:], strict=True):
            if ay == by == y and min(ax, bx) <= x <= max(ax, bx):
                return True
            if ax == bx == x and min(ay, by) <= y <= max(ay, by):
                return True
    return False


def check_data():
    """Every station and every dashed input lies on a line; stops with an error otherwise."""
    for label, x, y, _side in STATIONS:
        if not on_line(x, y):
            sys.exit(f"ERROR: station {label!r} at ({x}, {y}) is not on a line")
    for icon, _caption, x, y, _top in INPUTS:
        if not any(sx == x and sy == y for _l, sx, sy, _s in STATIONS):
            sys.exit(f"ERROR: input {icon!r} at ({x}, {y}) does not reach a station")


def file_icon(cx, y, label):
    """nf-core style file icon centred on cx, top at y, with a black label box."""
    w, h, fold = 46, 58, 14
    x = cx - w / 2
    return (
        f'<g stroke="#000" stroke-width="4" fill="#fff" stroke-linejoin="round">'
        f'<path d="M {x} {y} L {x + w - fold} {y} L {x + w} {y + fold} L {x + w} {y + h} L {x} {y + h} Z"/>'
        f'<path d="M {x + w - fold} {y} L {x + w - fold} {y + fold} L {x + w} {y + fold}" fill="#000"/></g>'
        f'<rect x="{x - 6}" y="{y + 22}" width="{w + 12}" height="22" rx="4" fill="#000"/>'
        f'<text x="{cx}" y="{y + 38}" font-family="{FONT}" font-size="14" font-weight="bold" '
        f'text-anchor="middle" fill="#fff">{escape(label.upper())}</text>'
    )


def caption(cx, y, text):
    """Small grey caption centred on cx."""
    return (
        f'<text x="{cx}" y="{y}" font-family="{FONT}" font-size="16" fill="#444" '
        f'text-anchor="middle">{escape(text)}</text>'
    )


def terminus(x, y):
    """Black vertical bar at a line start or end."""
    return (
        f'<line x1="{x}" y1="{y - 14}" x2="{x}" y2="{y + 14}" stroke="#000" stroke-width="12" '
        'stroke-linecap="round"/>'
    )


def multiline(x, y, text, anchor="middle", size=21, step=24):
    """Bold monospace label; y is the baseline of the first line."""
    spans = "".join(
        f'<tspan x="{x}" dy="{0 if i == 0 else step}">{escape(part)}</tspan>'
        for i, part in enumerate(text.split("\n"))
    )
    return (
        f'<text x="{x}" y="{y}" font-family="{MONO}" font-size="{size}" font-weight="bold" '
        f'text-anchor="{anchor}">{spans}</text>'
    )


def station(label, x, y, side):
    """Station circle and its label above or below."""
    n_lines = label.count("\n") + 1
    if side == "above":
        baseline = y - 24 - 24 * (n_lines - 1)
    else:
        baseline = y + 40
    return (
        f'<circle cx="{x}" cy="{y}" r="12" fill="#fff" stroke="#000" stroke-width="5"/>'
        + multiline(x, baseline, label)
    )


def dashed_input(icon, text, x, y, top):
    """File icon above or below a station, its caption, and a dashed line to the station."""
    out = [file_icon(x, top, icon)]
    if top < y:
        # icon above: caption under the icon, dashed line down to the station
        out.append(caption(x, top + 78, text))
        y0, y1 = top + 88, y - 12
    else:
        # icon below: caption under the icon, dashed line up to the station
        out.append(caption(x, top + 78, text))
        y0, y1 = y + 12, top - 4
    out.insert(
        0,
        f'<line x1="{x}" y1="{y0}" x2="{x}" y2="{y1}" stroke="{DASH_COLOUR}" stroke-width="4" '
        'stroke-dasharray="9 7"/>',
    )
    return out


def stage_box(number, word, x0, y0, x1, y1, corner, key):
    """Grey (key step: green) box with its number and name in one corner (tl/tr/bl/br)."""
    fill = KEY_FILL if key else "#e9e9e9"
    out = [f'<rect x="{x0}" y="{y0}" width="{x1 - x0}" height="{y1 - y0}" rx="26" fill="{fill}"/>']
    # number and name as one text; its width is estimated to align it on the right
    width = 60 + 22 * len(word)
    x = x1 - 22 - width if corner in ("tr", "br") else x0 + 20
    y = y1 - 22 if corner in ("bl", "br") else y0 + 76
    out.append(
        f'<text x="{x}" y="{y}" font-family="{FONT}" font-weight="bold">'
        f'<tspan font-size="80" fill="#fff">{number}</tspan>'
        f'<tspan dx="14" font-size="32" fill="#333">{escape(word)}</tspan></text>'
    )
    return out


def legend():
    """Stage words and symbols, in the free space below stage 1."""
    x, y0 = 50, 572
    out = [
        f'<text x="{x}" y="{y0}" font-family="{FONT}" font-size="22" font-weight="bold">STAGE</text>'
    ]
    for i, (number, word, hint, *_box) in enumerate(STAGES):
        out.append(
            f'<text x="{x}" y="{y0 + 34 + 30 * i}" font-family="{FONT}" font-size="19">'
            f'<tspan font-weight="bold">{number} {escape(word)}</tspan>'
            f'<tspan fill="#555"> · {escape(hint)}</tspan></text>'
        )
    y = y0 + 280
    out.append(
        f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="22" font-weight="bold">SYMBOL</text>'
    )
    out.append(
        f'<circle cx="{x + 14}" cy="{y + 27}" r="11" fill="#fff" stroke="#000" stroke-width="5"/>'
    )
    out.append(
        f'<text x="{x + 40}" y="{y + 34}" font-family="{FONT}" font-size="19">step (main tool)</text>'
    )
    out.append(
        f'<line x1="{x + 230}" y1="{y + 27}" x2="{x + 290}" y2="{y + 27}" stroke="{DASH_COLOUR}" '
        'stroke-width="4" stroke-dasharray="9 7"/>'
    )
    out.append(
        f'<text x="{x + 302}" y="{y + 34}" font-family="{FONT}" font-size="19">file read by a step</text>'
    )
    out.append(f'<rect x="{x}" y="{y + 52}" width="28" height="24" rx="6" fill="{KEY_FILL}"/>')
    out.append(
        f'<text x="{x + 40}" y="{y + 71}" font-family="{FONT}" font-size="19">'
        "key step: the two comparisons</text>"
    )
    return out


def build_svg():
    check_data()
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">',
        f'<rect width="{WIDTH}" height="{HEIGHT}" fill="#fff"/>',
        f'<text x="30" y="62" font-family="{FONT}" font-size="40" font-weight="bold">'
        f'<tspan fill="{LINE_COLOUR}">core-unit-bioinformatics/</tspan>'
        "workflow-nxf-somatic-variantcomp</text>",
    ]
    for number, word, _hint, *box in STAGES:
        out.extend(stage_box(number, word, *box))
    # dashed inputs below the lines
    for args in INPUTS:
        out.extend(dashed_input(*args))
    for points in LINES:
        out.append(
            f'<path d="{line_path(points)}" fill="none" stroke="{LINE_COLOUR}" '
            f'stroke-width="{LINE_WIDTH}" stroke-linejoin="round"/>'
        )
    # line starts and ends: terminus bar, file icon next to it, caption below the icon
    for icon, text, cx, y, x_bar in STARTS:
        out.append(terminus(x_bar, y))
        out.append(file_icon(cx, y - 29, icon))
        out.append(caption(cx, y + 50, text))
    for text, x, y in NOTES:
        out.append(
            f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="17" font-style="italic" '
            f'fill="#333">{escape(text)}</text>'
        )
    # stations last, on top of the lines
    for args in STATIONS:
        out.append(station(*args))
    out.extend(legend())
    out.append(
        f'<text x="{WIDTH - 30}" y="{HEIGHT - 14}" font-family="{FONT}" font-size="14" '
        'fill="#777" text-anchor="end">Metro-map style: nf-core components (CC0), '
        "J. A. Fellows Yates, M. U. Garcia, L. Le Nezet; after the nf-core/rnaseq map. "
        "Every step: README.md</text>"
    )
    out.append("</svg>")
    return "\n".join(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--outdir", default="docs/images", type=Path, help="output directory")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    svg = args.outdir / "variantcomp_metro_map_overview.svg"
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
