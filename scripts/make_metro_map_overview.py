#!/usr/bin/env python3
"""Draw the overview metro map of the pipeline (nf-core/rnaseq style).

The overview complements the detailed map of `make_metro_map.py`: it shows only the
stages, the main tools at the stations, and how the inputs run through the pipeline.

Layout: two rows read like a snake. Row 1 runs from left to right (input checks, VCF
normalization, benchmark), the lines turn down at the right edge, and row 2 runs back
from right to left (discordance and annotation, reports). Two branches leave the main
rows: the regions/coverage branch (target BEDs and alignments) below row 1 and the
known-variant branch below row 2. Every stage is a grey box with its number.

Lines: every input has its own colour and runs in a bundle with the others; it leaves
the bundle where it is used and ends at the station of its last use. Optional inputs
are switched on by their params (legend). A station is a circle on one line or a pill
across the lines it uses. File icons mark inputs (line starts) and outputs (line ends).

The style follows the nf-core metro-map components (CC0; James A. Fellows Yates,
Maxime U. Garcia, Louis Le Nezet; https://nf-co.re/docs/community/brand/workflow-schematics)
and the nf-core/rnaseq metro map.

Usage:
  scripts/make_metro_map_overview.py --outdir docs/images   SVG (+ PNG with cairosvg)

Edit LINES, STATIONS and STAGES when the workflow changes.
"""

import argparse
import math
import sys
from pathlib import Path
from xml.sax.saxutils import escape

FONT = "Maven Pro, Helvetica, Arial, sans-serif"
MONO = "DejaVu Sans Mono, Menlo, Consolas, monospace"

WIDTH = 2500
HEIGHT = 1530
# distance between parallel lines, line width and default corner radius of a bend
SPACE = 12
LINE_WIDTH = 9
RADIUS = 50

# centerline coordinates: row 1, regions branch, row 2, known-variant branch, right turn
Y1, Y3, Y2, Y6 = 330, 550, 900, 1070
XR = 2380

# lines: key -> (colour, offset from the centerline in units of SPACE, legend text);
# a positive offset lies right of the direction of travel (below on row 1)
TRACKS = {
    "known": (
        "#ECDC86",
        -3,
        "Known variants, optional (--known_variants): required checklist or truth set",
    ),
    "annot": (
        "#E879E3",
        -2,
        "GIAB and stratification BEDs, repeats, optional (--giab_<name>, --stratification)",
    ),
    "vep": ("#FA6863", -1, "Ensembl VEP cache or GTF, optional (--vep)"),
    "vcf": ("#24B064", 0, "VCFs of all technologies (samplesheet): the main line"),
    "bed": ("#3d95fd", 1, "Target BEDs (samplesheet): the common comparison region"),
    "bam": (
        "#F18046",
        2,
        "Alignments, optional (samplesheet; --callable_regions, --read_evidence)",
    ),
}
REF_COLOUR = "#3F2B29"
REF_TEXT = "Reference FASTA (--fasta): every input is checked against it"

# centerline of every line; a point may carry its own corner radius as third value
TOP_INPUTS = [(340, 205), (340, Y1)]
LINES = {
    "vcf": [(130, Y1), (XR, Y1), (XR, Y2), (570, Y2)],
    "bed": [(130, Y1), (760, Y1), (970, Y3), (1400, Y3), (1610, Y1), (1680, Y1)],
    "bam": [
        (130, Y1),
        (760, Y1),
        (970, Y3),
        (1400, Y3),
        (1610, Y1),
        (XR, Y1),
        (XR, Y2),
        (2090, Y2),
    ],
    "vep": TOP_INPUTS + [(XR, Y1), (XR, Y2), (1950, Y2)],
    # the stratification line moves next to the VCF line once the VEP line has ended
    "annot": TOP_INPUTS
    + [(XR, Y1), (XR, Y2), (1900, Y2, 16), (1888, Y2 - SPACE, 16), (1700, Y2 - SPACE)],
    # the known-variant line leaves row 2 downwards and joins the report next to the VCFs
    "known": TOP_INPUTS
    + [(XR, Y1), (XR, Y6), (1400, Y6), (1206, Y2 - 2 * SPACE), (1120, Y2 - 2 * SPACE)],
}

# stations: (label, x, y of the centerline segment, lines, label side); "\n" breaks lines
STATIONS = [
    # 1 input checks
    ("nf-schema", 190, Y1, ["vcf", "bed", "bam"], "above"),
    ("check_\nsheets.py", 280, Y1, ["vcf", "bed", "bam"], "below"),
    ("check_\nreference.py\n(pysam)", 620, Y1, list(TRACKS), "below"),
    # 2 VCF normalization and uniform filters
    ("bcftools\nview", 870, Y1, ["vcf"], "above"),
    ("bcftools\nannotate", 1010, Y1, ["vcf"], "above"),
    ("bcftools\nreheader", 1150, Y1, ["vcf"], "above"),
    ("bcftools\nnorm", 1290, Y1, ["vcf"], "above"),
    ("harmonize_\nvcf.py", 1430, Y1, ["vcf"], "above"),
    ("bcftools\nview", 1680, Y1, ["vcf", "bed", "bam"], "above"),
    ("variant_\naccounting.py", 1790, Y1, ["vcf"], "below"),
    # 3 comparison regions and coverage
    ("samtools\nindex", 1020, Y3, ["bam"], "below"),
    ("bedtools\nmultiinter", 1110, Y3, ["bed"], "above"),
    ("mosdepth", 1210, Y3, ["bed", "bam"], "below"),
    ("bedtools\nintersect", 1310, Y3, ["bed", "bam"], "above"),
    # 4 benchmark
    ("rtg\nvcfeval", 1950, Y1, ["vcf"], "above"),
    ("bcftools\nisec", 2060, Y1, ["vcf"], "below"),
    ("benchmark_\nmetrics.py", 2160, Y1, ["vcf"], "above"),
    ("aggregate_\nmetrics.py", 2260, Y1, ["vcf"], "below"),
    # 5 discordance and annotation
    ("bcftools\nmpileup", 2230, Y2, ["bam"], "above"),
    ("classify_\ndiscordance.py", 2090, Y2, ["vcf", "bam"], "below"),
    ("Ensembl\nVEP", 1950, Y2, ["vcf", "vep"], "above"),
    ("msisensor-\npro", 1820, Y2, ["annot"], "below"),
    ("bedtools\nintersect", 1700, Y2, ["vcf", "annot"], "above"),
    ("annotate_\nvariants.py", 1580, Y2, ["vcf"], "below"),
    # 6 known variants
    ("bcftools norm", 2250, Y6, ["known"], "below"),
    ("known_\nvariants.py", 2080, Y6, ["known"], "below"),
    ("rtg vcfeval", 1900, Y6, ["known"], "below"),
    ("known_\nsummary.py", 1720, Y6, ["known"], "below"),
    # 7 reports
    ("variantcomp_\nreport.py", 1330, Y2, ["vcf"], "above"),
    ("consistency_\nchecks.py", 1120, Y2, ["vcf", "known"], "below"),
    ("datavzrd", 960, Y2, ["vcf"], "above"),
    ("report_\nhtml.py", 830, Y2, ["vcf"], "below"),
    ("MultiQC", 700, Y2, ["vcf"], "above"),
]

# stage boxes: (number, name, x0, y0, x1, y1, number position "tl"/"tr"/"bl")
STAGES = [
    (1, "Input checks: samplesheet, reference, every input", 30, 110, 720, 480, "tl"),
    (2, "VCF normalization and uniform filters", 740, 110, 1870, 420, "tl"),
    (3, "Comparison regions and coverage", 740, 440, 1640, 660, "bl"),
    (4, "Benchmark against the truth technology", 1890, 110, 2470, 420, "tr"),
    (5, "Discordance and annotation", 1480, 740, 2470, 1000, "tl"),
    (6, "Known variants (required / benchmark)", 1480, 1015, 2470, 1200, "bl"),
    (7, "Reports and end-of-run checks", 30, 740, 1460, 1000, "tl"),
]


def unit(ax, ay, bx, by):
    """Unit vector from point a to point b."""
    length = math.hypot(bx - ax, by - ay)
    return (bx - ax) / length, (by - ay) / length


def offset_vertices(points, offset):
    """Vertices of the line `offset` px right of the centerline (miter joins).

    Returns a list of (x, y, radius, (ux, uy) in, (ux, uy) out) per vertex.
    """
    xy = [(p[0], p[1]) for p in points]
    dirs = [unit(*xy[i], *xy[i + 1]) for i in range(len(xy) - 1)]
    out = []
    for i, (x, y) in enumerate(xy):
        u_in = dirs[i - 1] if i > 0 else dirs[0]
        u_out = dirs[i] if i < len(dirs) else dirs[-1]
        # right-hand normals (screen coordinates) of the incoming and outgoing segment
        n_in = (-u_in[1], u_in[0])
        n_out = (-u_out[1], u_out[0])
        # the miter vector has projection 1 on both normals
        dot = n_in[0] * n_out[0] + n_in[1] * n_out[1]
        mx, my = (n_in[0] + n_out[0]) / (1 + dot), (n_in[1] + n_out[1]) / (1 + dot)
        radius = points[i][2] if len(points[i]) > 2 else RADIUS
        # the inner line of a bend gets a smaller radius, so parallel lines stay concentric
        cross = u_in[0] * u_out[1] - u_in[1] * u_out[0]
        if cross > 0:
            radius -= offset
        elif cross < 0:
            radius += offset
        out.append((x + offset * mx, y + offset * my, max(radius, 4), u_in, u_out, cross))
    return out


def line_path(points, offset):
    """SVG path of one line with rounded bends."""
    vertices = offset_vertices(points, offset)
    d = [f"M {vertices[0][0]:.1f} {vertices[0][1]:.1f}"]
    for x, y, radius, u_in, u_out, cross in vertices[1:-1]:
        # turn angle and distance of the arc ends from the vertex
        angle = math.acos(max(-1, min(1, u_in[0] * u_out[0] + u_in[1] * u_out[1])))
        t = radius * math.tan(angle / 2)
        sweep = 1 if cross > 0 else 0
        d.append(f"L {x - t * u_in[0]:.1f} {y - t * u_in[1]:.1f}")
        d.append(
            f"A {radius:.1f} {radius:.1f} 0 0 {sweep} {x + t * u_out[0]:.1f} {y + t * u_out[1]:.1f}"
        )
    d.append(f"L {vertices[-1][0]:.1f} {vertices[-1][1]:.1f}")
    return " ".join(d)


def y_at(track, x, y_hint):
    """y of a line at x on its horizontal segment closest to y_hint (None if absent)."""
    vertices = offset_vertices(LINES[track], SPACE * TRACKS[track][1])
    best = None
    for (ax, ay, *_), (bx, by, *_) in zip(vertices[:-1], vertices[1:], strict=True):
        if abs(ay - by) < 0.5 and min(ax, bx) - 1 <= x <= max(ax, bx) + 1:
            if best is None or abs(ay - y_hint) < abs(best - y_hint):
                best = ay
    if best is not None and abs(best - y_hint) > 60:
        return None
    return best


def check_data():
    """Every station lies on all its lines; stops with an error otherwise."""
    for label, x, y, tracks, _side in STATIONS:
        for track in tracks:
            if y_at(track, x, y) is None:
                sys.exit(f"ERROR: station {label!r} at x={x}: line {track!r} does not pass")


def file_icon(cx, y, label, colour="#000", text_colour="#fff"):
    """nf-core style file icon centred on cx, top at y, with a label box."""
    w, h, fold = 46, 58, 14
    x = cx - w / 2
    return (
        f'<g stroke="#000" stroke-width="4" fill="#fff" stroke-linejoin="round">'
        f'<path d="M {x} {y} L {x + w - fold} {y} L {x + w} {y + fold} L {x + w} {y + h} L {x} {y + h} Z"/>'
        f'<path d="M {x + w - fold} {y} L {x + w - fold} {y + fold} L {x + w} {y + fold}" fill="#000"/></g>'
        f'<rect x="{x - 6}" y="{y + 22}" width="{w + 12}" height="22" rx="4" fill="{colour}" stroke="#000" stroke-width="2"/>'
        f'<text x="{cx}" y="{y + 38}" font-family="{FONT}" font-size="14" font-weight="bold" '
        f'text-anchor="middle" fill="{text_colour}">{escape(label.upper())}</text>'
    )


def caption(cx, y, text, size=15, anchor="middle"):
    """Small grey caption under an icon."""
    return (
        f'<text x="{cx}" y="{y}" font-family="{FONT}" font-size="{size}" fill="#444" '
        f'text-anchor="{anchor}">{escape(text)}</text>'
    )


def terminus(x0, y0, x1, y1):
    """Black rounded bar at a line start or end."""
    return (
        f'<line x1="{x0}" y1="{y0}" x2="{x1}" y2="{y1}" stroke="#000" stroke-width="12" '
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


def station(label, x, y, tracks, side):
    """Station (circle or pill) and its label above or below the whole bundle at x."""
    ys = [y_at(t, x, y) for t in tracks]
    top, bottom = min(ys), max(ys)
    out = []
    if top == bottom:
        out.append(
            f'<circle cx="{x}" cy="{top}" r="11" fill="#fff" stroke="#000" stroke-width="5"/>'
        )
    else:
        out.append(
            f'<rect x="{x - 11}" y="{top - 11}" width="22" height="{bottom - top + 22}" rx="11" '
            'fill="#fff" stroke="#000" stroke-width="5"/>'
        )
    # the label clears every line that passes the station, not only the ones it uses
    passing = [v for v in (y_at(t, x, y) for t in TRACKS) if v is not None]
    n_lines = label.count("\n") + 1
    if side == "above":
        baseline = min(passing + [top]) - 22 - 24 * (n_lines - 1)
    else:
        baseline = max(passing + [bottom]) + 38
    out.append(multiline(x, baseline, label))
    return "\n".join(out)


def build_svg():
    check_data()
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">',
        f'<rect width="{WIDTH}" height="{HEIGHT}" fill="#fff"/>',
        f'<text x="40" y="72" font-family="{FONT}" font-size="42" font-weight="bold">'
        f'<tspan fill="#24B064">core-unit-bioinformatics/</tspan>workflow-nxf-somatic-variantcomp</text>',
    ]
    # stage boxes with their numbers
    for number, _name, x0, y0, x1, y1, corner in STAGES:
        out.append(
            f'<rect x="{x0}" y="{y0}" width="{x1 - x0}" height="{y1 - y0}" rx="26" fill="#e9e9e9"/>'
        )
        nx, anchor = (x1 - 22, "end") if corner == "tr" else (x0 + 22, "start")
        ny = y1 - 18 if corner == "bl" else y0 + 82
        out.append(
            f'<text x="{nx}" y="{ny}" font-family="{FONT}" font-size="86" font-weight="bold" '
            f'fill="#fff" text-anchor="{anchor}">{number}</text>'
        )
    # reference: FASTA icon, its own short line with the indexing station, plugged into the bundle
    out.append(
        f'<line x1="540" y1="205" x2="540" y2="282" stroke="{REF_COLOUR}" stroke-width="{LINE_WIDTH}"/>'
    )
    out.append(terminus(526, 284, 554, 284))
    out.append(terminus(526, 205, 554, 205))
    out.append(file_icon(540, 132, "fasta"))
    out.append('<circle cx="540" cy="245" r="11" fill="#fff" stroke="#000" stroke-width="5"/>')
    out.append(multiline(560, 241, "samtools\nfaidx", anchor="start"))
    # lines, the outer ones first so that the main line is on top where lines cross
    for track in ["known", "annot", "vep", "bam", "bed", "vcf"]:
        colour, offset, _text = TRACKS[track]
        out.append(
            f'<path d="{line_path(LINES[track], SPACE * offset)}" fill="none" stroke="{colour}" '
            f'stroke-width="{LINE_WIDTH}" stroke-linejoin="round"/>'
        )
    # inputs: samplesheet (VCFs, BEDs, alignments) on the left, params from the top
    out.append(terminus(130, Y1 - 6, 130, Y1 + 2 * SPACE + 6))
    out.append(file_icon(72, Y1 - 18, "csv"))
    out.append(caption(72, Y1 + 62, "samplesheet"))
    out.append(caption(72, Y1 + 80, "VCF, BED, BAM"))
    out.append(terminus(346, 205, 382, 205))
    for cx, label, name in [(290, "gtf", "VEP"), (364, "bed", "strata"), (438, "csv", "known")]:
        out.append(file_icon(cx, 132, label))
        out.append(caption(cx, 126, name))
    # outputs: the reports at the end of the main line
    out.append(terminus(570, Y2 - 6, 570, Y2 + 6))
    for cx, name in [(500, "report"), (420, "tables"), (340, "MultiQC")]:
        out.append(file_icon(cx, Y2 - 29, "html"))
        out.append(caption(cx, Y2 + 52, name))
    # stations last, on top of the lines
    for args in STATIONS:
        out.append(station(*args))
    out.extend(legend())
    out.append("</svg>")
    return "\n".join(out)


def legend():
    """Stage list, line colours and the meaning of the symbols, below the map."""
    y0 = 1260
    out = [
        f'<text x="40" y="{y0}" font-family="{FONT}" font-size="24" font-weight="bold">STAGE</text>',
        f'<text x="760" y="{y0}" font-family="{FONT}" font-size="24" font-weight="bold">INPUT LINE</text>',
    ]
    for i, (number, name, *_box) in enumerate(STAGES):
        out.append(
            f'<text x="40" y="{y0 + 32 + 30 * i}" font-family="{FONT}" font-size="19">{number}. {escape(name)}</text>'
        )
    entries = [(c, text) for c, _o, text in TRACKS.values()][::-1]
    entries = [entries[0], entries[1]] + entries[2:] + [(REF_COLOUR, REF_TEXT)]
    for i, (colour, text) in enumerate(entries):
        y = y0 + 32 + 30 * i
        out.append(
            f'<line x1="760" y1="{y - 7}" x2="830" y2="{y - 7}" stroke="{colour}" stroke-width="{LINE_WIDTH}"/>'
        )
        out.append(
            f'<text x="846" y="{y}" font-family="{FONT}" font-size="19">{escape(text)}</text>'
        )
    x = 1830
    out.append(
        f'<text x="{x}" y="{y0}" font-family="{FONT}" font-size="24" font-weight="bold">SYMBOL</text>'
    )
    out.append(
        f'<circle cx="{x + 14}" cy="{y0 + 25}" r="11" fill="#fff" stroke="#000" stroke-width="5"/>'
    )
    out.append(
        f'<text x="{x + 46}" y="{y0 + 32}" font-family="{FONT}" font-size="19">step on one input</text>'
    )
    out.append(
        f'<rect x="{x + 3}" y="{y0 + 50}" width="22" height="46" rx="11" fill="#fff" stroke="#000" stroke-width="5"/>'
    )
    out.append(
        f'<text x="{x + 46}" y="{y0 + 80}" font-family="{FONT}" font-size="19">step that combines inputs</text>'
    )
    out.append(terminus(x + 14, y0 + 116, x + 14, y0 + 140))
    out.append(
        f'<text x="{x + 46}" y="{y0 + 135}" font-family="{FONT}" font-size="19">input or output file</text>'
    )
    out.append(
        f'<text x="{x}" y="{y0 + 178}" font-family="{FONT}" font-size="17" fill="#555">'
        "A line ends at its last use. All steps: README.md</text>"
    )
    out.append(
        f'<text x="40" y="{HEIGHT - 18}" font-family="{FONT}" font-size="15" fill="#666">'
        "Metro-map style: nf-core components (CC0), J. A. Fellows Yates, M. U. Garcia, L. Le Nezet; "
        "after the nf-core/rnaseq metro map. Detailed map: variantcomp_metro_map.svg</text>"
    )
    return out


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
