#!/usr/bin/env python3
"""Draw the nf-core-style metro map of the pipeline and its step table.

Layout: the pipeline runs from top to bottom; every subworkflow is a box and every row
in a box is one step with its number, what happens and the tools.

Lines: the green central line is the VCF of every technology, from the input check to
the report. Every other input has its own line: it starts at its file icon, runs through
the steps that prepare it and then merges into the step that uses it (a diagonal into
the station), so that each line can be followed as one path. A line can fork (the
alignments are used for the callable regions and for the read evidence). The reference
(FASTA) is used by many steps: after its own steps it feeds them with dashed lines.

The style follows the nf-core metro-map components (CC0; James A. Fellows Yates,
Maxime U. Garcia, Louis Le Nezet; https://nf-co.re/docs/community/brand/workflow-schematics):
nf-core line colours, white stations with black borders, file icons at the line starts.

Usage:
  scripts/make_metro_map.py --outdir docs/images     SVG (+ PNG with cairosvg)
  scripts/make_metro_map.py --markdown               step table for the README

Edit BANDS, LINES and FEEDS when the workflow changes; map and README table come from it.
"""

import argparse
import sys
from pathlib import Path
from xml.sax.saxutils import escape

FONT = "Maven Pro, Helvetica, Arial, sans-serif"
MONO = "DejaVu Sans Mono, Menlo, Consolas, monospace"

# tracks (vertical lines) from left to right: key -> (colour, icon label, legend text)
TRACKS = {
    "ref": ("#3F2B29", "fasta", "Reference genome (FASTA): feeds the steps that use it (dashed)"),
    "bed": ("#3d95fd", "bed", "Target BEDs (+ --regions_bed / --exclude_bed): the common region"),
    "core": ("#24B064", "vcf", "VCFs of all technologies: the central line to the results"),
    "bam": (
        "#F18046",
        "bam",
        "Alignments: callable regions (--callable_regions), read evidence (--read_evidence)",
    ),
    "annot": (
        "#E879E3",
        "bed",
        "GIAB (--giab_<name>) and other stratification BEDs, repeat regions",
    ),
    "vep": ("#FA6863", "gtf", "VEP cache or GTF (--vep): consequences, genes"),
    "known": (
        "#ECDC86",
        "csv",
        "Known variants (--known_variants): required checklist, truth sets",
    ),
}

# subworkflow boxes: (name, purpose, steps); step: (number, text, tools, track)
BANDS = [
    (
        "PIPELINE_INITIALISATION, CHECK_SAMPLESHEET",
        "read and check the samplesheet",
        [
            ("1.1", "Read the samplesheet, check every column", "nf-schema", "core"),
            (
                "1.2",
                "Check rows: truth technology, one BED per technology",
                "check_sheets.py",
                "core",
            ),
        ],
    ),
    (
        "PREPARE_REFERENCE",
        "index the reference genome",
        [
            ("2.1", "Unpack and index the reference FASTA", "gunzip, samtools faidx", "ref"),
            ("2.2", "Genome BED: one interval per compared contig", "gawk", "ref"),
        ],
    ),
    (
        "CHECK_INPUTS",
        "check every input against the reference, before anything else",
        [
            (
                "3.1",
                "VCFs: contigs, lengths, samples, REF alleles",
                "check_reference.py (pysam)",
                "core",
            ),
            ("3.2", "Alignment headers: contigs and lengths", "check_reference.py (pysam)", "bam"),
            ("3.3", "Target, regions and exclude BEDs", "check_reference.py", "bed"),
            (
                "3.4",
                "GIAB and stratification BEDs, --repeats_bed",
                "check_sheets.py, check_reference.py",
                "annot",
            ),
            (
                "3.5",
                "Known-variant sheet and lists: positions, REF alleles",
                "check_sheets.py, known_variants.py",
                "known",
            ),
        ],
    ),
    (
        "PREPARE_ALIGNMENTS",
        "index the alignments",
        [
            ("4.1", "Index alignments without an index", "samtools index", "bam"),
        ],
    ),
    (
        "PREPARE_REGIONS",
        "build the common comparison region",
        [
            (
                "5.1",
                "Sort and intersect all target BEDs",
                "bedtools sort, multiinter, gawk, merge",
                "bed",
            ),
            (
                "5.2",
                "Restrict to --regions_bed, remove --exclude_bed",
                "bedtools intersect, subtract",
                "bed",
            ),
            (
                "5.3",
                "Common region (genome if no BED) and its size",
                "bedtools merge, region_summary.py",
                "bed",
            ),
            ("5.4", "Size of the region after every step", "region_summary.py steps", "bed"),
        ],
    ),
    (
        "COVERAGE",
        "find the bases callable in all alignments of a sample",
        [
            (
                "6.1",
                "Common region in the alignment's contig names",
                "check_reference.py rename-bed",
                "bam",
            ),
            ("6.2", "Depth per base and per common-region interval", "mosdepth", "bam"),
            (
                "6.3",
                "Callable bases per alignment (tumor / normal threshold)",
                "mosdepth_callable.py",
                "bam",
            ),
            (
                "6.4",
                "Bases callable in all alignments of a sample",
                "bedtools multiinter, intersect",
                "bam",
            ),
            (
                "6.5",
                "Depth in the common and the target region",
                "mosdepth, coverage_summary.py",
                "bam",
            ),
        ],
    ),
    (
        "PREPROCESS_VCF, NORMALIZE_VCF",
        "normalize every VCF in the same way",
        [
            ("7.1", "Keep tumor (+ normal) columns, recompress with BGZF", "bcftools view", "core"),
            (
                "7.2",
                "Split PASS and non-PASS records (--filter_pass_only)",
                "bcftools view",
                "core",
            ),
            (
                "7.3",
                "Rename contigs; remove other contigs and gVCF blocks",
                "bcftools annotate, view",
                "core",
            ),
            ("7.4", "Rewrite the contig header from the reference", "bcftools reheader", "core"),
            (
                "7.5",
                "Split multiallelics, left-align, split MNVs, deduplicate",
                "bcftools norm",
                "core",
            ),
            ("7.6", "Same VAF, depth, alt reads for every caller", "harmonize_vcf.py", "core"),
            ("7.7", "All records together (PASS + non-PASS)", "bcftools concat", "core"),
            (
                "7.8",
                "Uniform filters; keep the common (or callable) region",
                "bcftools view",
                "core",
            ),
            (
                "7.9",
                "Records kept and removed per step and reason",
                "variant_accounting.py",
                "core",
            ),
        ],
    ),
    (
        "BENCHMARK",
        "compare every technology with the truth technology",
        [
            ("8.1", "Match calls haplotype-aware: TP, FN, FP", "rtg vcfeval", "core"),
            ("8.2", "Exact-allele and same-position matches", "bcftools isec", "core"),
            (
                "8.3",
                "Recall, precision, F1 per sample (overall, type, VAF)",
                "benchmark_metrics.py",
                "core",
            ),
            (
                "8.4",
                "Combine samples: summed (micro) and mean (macro)",
                "aggregate_metrics.py",
                "core",
            ),
        ],
    ),
    (
        "DISCORDANCE",
        "explain why variants differ",
        [
            (
                "9.1",
                "Read counts of every alignment at the compared positions",
                "bcftools query, mpileup + call, norm",
                "bam",
            ),
            ("9.2", "Reason for every FN/FP; detection power", "classify_discordance.py", "core"),
        ],
    ),
    (
        "ANNOTATE_VARIANTS",
        "annotate every compared variant",
        [
            ("10.1", "Compared variants as sites VCF", "annotate_variants.py to-vcf", "core"),
            (
                "10.2",
                "Consequence, gene, impact, variant class",
                "Ensembl VEP, bcftools +split-vep",
                "vep",
            ),
            ("10.3", "Repeat regions of the reference", "msisensor-pro scan", "annot"),
            (
                "10.4",
                "Repeat context and stratifications per variant",
                "bedtools map, intersect -c",
                "annot",
            ),
            (
                "10.5",
                "Merge all annotations into the variant table",
                "annotate_variants.py merge",
                "core",
            ),
        ],
    ),
    (
        "PREPARE_KNOWN_VARIANTS, KNOWN_VARIANTS",
        "check known variants",
        [
            ("11.1", "Normalize the lists like the calls", "bcftools norm", "known"),
            (
                "11.2",
                "Required: PASS/FAIL per sample and technology",
                "known_variants.py check",
                "known",
            ),
            ("11.3", "Required: summary and recall over all samples", "known_summary.py", "known"),
            (
                "11.4",
                "Benchmark: every technology vs the truth set",
                "rtg vcfeval, bcftools isec, benchmark_metrics.py",
                "known",
            ),
        ],
    ),
    (
        "REPORT",
        "summarize",
        [
            (
                "12.1",
                "Enrichment, gene and annotation tables, plots",
                "variantcomp_report.py",
                "core",
            ),
            (
                "12.2",
                "Consistency checks of all counts",
                "consistency_checks.py",
                "core",
            ),
            (
                "12.3",
                "Narrative report and table browser",
                "report_html.py, datavzrd",
                "core",
            ),
            ("12.4", "MultiQC report with every step summary", "MultiQC, mqc_tables.py", "core"),
            ("12.5", "Stop the run if a check failed", "consistency_checks.py gate", "core"),
        ],
    ),
]

# lines as paths through steps; a path that ends on another track merges into that station
LINES = [
    (
        "core",
        [
            "1.1",
            "1.2",
            "3.1",
            "7.1",
            "7.2",
            "7.3",
            "7.4",
            "7.5",
            "7.6",
            "7.7",
            "7.8",
            "7.9",
            "8.1",
            "8.2",
            "8.3",
            "8.4",
            "9.2",
            "10.1",
            "10.5",
            "12.1",
            "12.2",
            "12.3",
            "12.4",
            "12.5",
        ],
    ),
    ("ref", ["2.1", "2.2"]),
    ("bed", ["3.3", "5.1", "5.2", "5.3", "5.4", "6.1"]),
    ("bam", ["3.2", "4.1", "6.1", "6.2", "6.3", "6.4", "7.8"]),
    ("bam", ["6.4", "6.5", "12.3"]),
    ("bam", ["6.4", "9.1", "9.2"]),
    ("annot", ["3.4", "10.3", "10.4", "10.5"]),
    ("vep", ["10.2", "10.5"]),
    ("known", ["3.5", "11.1", "11.2", "11.3", "11.4", "12.2"]),
]
# the reference feeds these steps (dashed)
FEEDS = ["3.1", "5.1", "6.2", "7.5", "8.1", "10.2"]

# longer explanation of every step for the README table (the map shows the short text)
DETAILS = {
    "1.1": "The samplesheet is read as in the nf-core template: nf-schema checks every column (names, library types, file extensions). Relative paths are resolved against the samplesheet's directory; VCF and alignment indices are looked up next to the files.",
    "1.2": "Rules across rows: every sample has the truth technology and another technology, no duplicate rows, one target BED and library type per technology, alignments where the options need them. All problems are reported at once with line numbers.",
    "2.1": "A gzipped FASTA is unpacked and indexed; the `.fai` (contig names and lengths) is the reference every input is checked against.",
    "2.2": "One interval per compared contig (the mitochondrial contig is left out unless `--keep_mito`); the comparison region of WGS-only runs.",
    "3.1": "Contig names (safe renames chr/no-chr, MT/chrM), contig lengths (any difference is an error: another build), the tumor/normal columns and the REF allele of every record against the FASTA (`--max_ref_mismatch`). Switch off with `--check_vcf_contigs` / `--check_ref_alleles`.",
    "3.2": "Contig names and lengths of every alignment header (`--check_bam_contigs`).",
    "3.3": "BED format, contig names and interval ends of the target BEDs, `--regions_bed` and `--exclude_bed` (`--check_bed_contigs`).",
    "3.4": "The GIAB genome stratifications given as `--giab_<name>` files, the `--stratification` sheet (columns, unique names, files) and every stratification BED; a precomputed `--repeats_bed`.",
    "3.5": "The `--known_variants` sheet (columns, modes, samples, files) and every list: contigs, positions and REF alleles against the FASTA (`--check_known_variants`). Every later step starts only after all checks of 1.2 and 3.x passed.",
    "4.1": "Alignments without a `.bai`/`.csi`/`.crai` next to them are indexed once.",
    "5.1": "Every target BED is sorted, and only the bases targeted by every technology are kept (bedtools multiinter; WGS counts as genome-wide).",
    "5.2": "Optional: restriction to `--regions_bed` and removal of `--exclude_bed`.",
    "5.3": "The common comparison region (`regions/common_targets.bed`), merged, and its size per technology; it must have at least `--min_region_bp` bp.",
    "6.1": "The common region is translated to the alignment's contig names (e.g. `1` instead of `chr1`).",
    "6.2": "mosdepth classifies every base as below or above the threshold and computes the mean depth per common-region interval (duplicates kept for amplicon libraries).",
    "6.3": "Callable bases of each alignment: tumor depth ≥ `--min_coverage_tumor`, normal depth ≥ `--min_coverage_normal` (base mode; region mode uses the mean depth).",
    "5.4": "Size of the comparison region after every step (genome, intersection of the target BEDs, `--regions_bed`, `--exclude_bed`), in bp and percent of the genome (`regions/common_targets.steps.tsv`).",
    "6.5": "Depth of every alignment in the common region and in the target region of its technology (genome for WGS): mean and percent of bases ≥ each of `--coverage_thresholds`; also in the MultiQC General Statistics.",
    "6.4": "Per sample: bases callable in every alignment of every technology, inside the common region (at least `--min_region_bp`). The sample is compared only there, so missing coverage in one technology is not counted as a miss.",
    "7.1": "Only the tumor (and normal) columns are kept, in this order; written as bgzipped VCF with index (inputs may be plain gzip).",
    "7.2": "With `--filter_pass_only` (default): PASS records (FILTER `PASS`/`.`) are compared; the other records go through the same steps so that the discordance step can name their FILTER.",
    "7.3": "Contigs are renamed to the reference names; records on other contigs (and the mitochondrial contig) and gVCF reference blocks are removed.",
    "7.4": "Contig header lines are rewritten from the reference index (adds missing contig lines).",
    "7.5": "Multiallelic records are split, indels left-aligned, MNVs split into SNVs, records whose REF differs from the FASTA and duplicates removed: the same representation for every caller.",
    "7.6": "VAF, depth and alt reads are computed the same way for every caller (`--vaf_source`: FORMAT/AF, else AD) and stored as `VC_*` INFO tags; GT is set to tumor 0/1, normal 0/0.",
    "7.9": "For every input VCF: records kept, removed (with the reason: non-PASS, contig, gVCF block, REF mismatch, duplicate, symbolic, region, VAF/depth/alt-read filter) and added by splitting, per step; every removed record; checks that the steps add up.",
    "7.7": "PASS and non-PASS records together: the normalized VCF that the discordance step searches.",
    "7.8": "The same filters for every technology (`--min_vaf`, `--min_depth`, `--min_alt_reads`) and the comparison region: the common region, or with `--callable_regions` the sample's callable region.",
    "8.1": "Each technology is matched with the truth technology of the same sample, haplotype-aware, so different representations of one variant still match (TP, FN = missed, FP = only in the query).",
    "8.2": "Record-by-record comparison by CHROM/POS/REF/ALT and by position only.",
    "8.3": "Recall, precision and F1 per sample and technology: overall, SNV vs. indel, by truth VAF (detection limit) and query VAF (low-VAF false positives), exact and position matches.",
    "8.4": "Over all samples: from the summed counts (micro) and as mean ± sd of the per-sample values (macro). Main result: `benchmark/metrics.aggregated.tsv`.",
    "9.1": "With `--read_evidence` (default): the positions of all compared variants are looked up in every tumor and normal alignment (MQ ≥ `--evidence_min_mapq`, BQ ≥ `--evidence_min_baseq`), normalized like the calls.",
    "9.2": "Each FN/FP gets a reason from the other technology's view: called but non-PASS, called but filtered by the uniform filters, germline (alt reads in the normal), low coverage, alt reads below the caller's threshold, no alt reads, or absent. For truth variants, the binomial detection power gives the power-adjusted recall.",
    "10.1": "All compared variants of all comparisons as one sorted sites-only VCF.",
    "10.2": "Ensembl VEP (offline; `--vep_cache` or `--vep_gtf`, options as in cio-abcd/variantinterpretation): the most severe consequence and its gene, impact, variant class and biotype per variant.",
    "10.3": "msisensor-pro scans the reference for repeat (microsatellite) regions, or `--repeats_bed` is used (`--annotate_repeats`).",
    "10.4": "Every variant gets its repeat unit and length and a yes/no per GIAB and stratification BED.",
    "10.5": "Result: `discordance/variants.annotated.tsv` (consequence group synonymous / non-synonymous / other, gene, impact, variant class, repeat, strata).",
    "11.1": "The checked known-variant lists are normalized like the calls (split, left-align, atomize).",
    "11.2": "Required mode: each listed variant must be found in every technology of the sample: PASS, or FAIL as filtered, outside the comparison region or not called (with read evidence if available).",
    "11.3": "Required mode: one table over all lists, samples and technologies with recall per technology, a warning per failed variant, optionally a stop (`--fail_on_missing_required`).",
    "11.4": "Benchmark mode: every technology is compared with the external truth set inside its confident regions, with the metrics of the main benchmark.",
    "12.1": "Are discordant variants enriched for indels, low VAF, non-synonymous consequences, repeats, strata or genes (Fisher's exact test, Benjamini-Hochberg from SciPy)? Discordance per annotation value, gene table, plots.",
    "12.2": "At the end: do the counts of all steps add up (accounting per VCF, metrics vs. variant tables, aggregation over samples, discordance classes, regions, required variants)? `report/consistency_checks.tsv`.",
    "12.3": "`report/variantcomp_report.html`: Summary with the key plots, Introduction, Results (region steps, input variants and their filtering, metrics, truth variants, every discordance class explained, annotation, known variants) and Methods (tools, parameters, coverage, checks); `report/tables/`: datavzrd browser of all variant tables.",
    "12.4": "One report with metrics, discordance, annotation, required variants, every step summary (input checks, bcftools norm, harmonization, regions and their steps, variant accounting, coverage, consistency checks) and software versions.",
    "12.5": "A failed consistency check stops the run with an error, after all reports were written.",
}

# geometry
TRACK_X = {"ref": 80, "bed": 145, "core": 230, "bam": 315, "annot": 380, "vep": 445, "known": 510}
LABEL_X = 575
TOOLS_X = 1170
WIDTH = 1640
ROW = 38
BAND_HEAD = 44
BAND_GAP = 14
TOP = 340
ICON_Y = 215


def layout():
    """Positions: boxes [(y0, y1, band)], {step number: (y, step)}, total height."""
    boxes, steps = [], {}
    y = TOP
    for band in BANDS:
        y0 = y
        y += BAND_HEAD
        for step in band[2]:
            steps[step[0]] = (y + ROW / 2, step)
            y += ROW
        y += 10
        boxes.append((y0, y, band))
        y += BAND_GAP
    return boxes, steps, y + 40


def check_data(steps):
    """Steps of LINES/FEEDS exist; every step lies on a line; DETAILS for every step."""
    on_lines = set()
    for track, path in LINES:
        if track not in TRACKS:
            sys.exit(f"ERROR: unknown track {track}")
        for number in path:
            if number not in steps:
                sys.exit(f"ERROR: line {track}: unknown step {number}")
        on_lines.update(path)
    for number in FEEDS:
        if number not in steps:
            sys.exit(f"ERROR: FEEDS: unknown step {number}")
    missing = sorted(set(steps) - on_lines, key=lambda n: [int(x) for x in n.split(".")])
    if missing:
        sys.exit(f"ERROR: steps on no line: {', '.join(missing)}")
    no_details = [n for n in steps if n not in DETAILS]
    if no_details:
        sys.exit(f"ERROR: no DETAILS for steps {', '.join(no_details)}")


def path_points(track, path, steps):
    """Points of a line: vertical on its track, diagonal into a station on another track."""
    x_line = TRACK_X[track]
    points = [(x_line, ICON_Y + 50)] if path[0] == first_step(track) else []
    for number in path:
        y, step = steps[number]
        x_station = TRACK_X[step[3]]
        if x_station == x_line:
            points.append((x_line, y))
        else:
            # merge: down the own track to just above the station, then into it
            points.append((x_line, y - ROW * 0.9))
            points.append((x_station, y))
    return points


def first_step(track):
    """The first step of the first path of a track (the line starts at its icon there)."""
    return next(path[0] for t, path in LINES if t == track)


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
        f'text-anchor="middle" fill="{"#fff" if colour == "#3F2B29" else "#000"}">{escape(label)}</text>'
    )


def build_svg():
    boxes, steps, height = layout()
    check_data(steps)
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{height}" viewBox="0 0 {WIDTH} {height}">',
        f'<rect width="{WIDTH}" height="{height}" fill="#fff"/>',
        f'<text x="40" y="56" font-family="{FONT}" font-size="30" font-weight="bold">'
        "core-unit-bioinformatics/workflow-nxf-somatic-variantcomp</text>",
        f'<text x="40" y="88" font-family="{FONT}" font-size="17" fill="#444">'
        "Each box is a subworkflow, each row a step (top to bottom). The green line is the VCFs; every "
        "other input line merges into the step that uses it.</text>",
    ]
    # legend
    for i, (colour, _icon, text) in enumerate(TRACKS.values()):
        y = 128 + 28 * i
        out.append(
            f'<line x1="{LABEL_X}" y1="{y}" x2="{LABEL_X + 50}" y2="{y}" stroke="{colour}" stroke-width="10" stroke-linecap="round"/>'
        )
        out.append(
            f'<text x="{LABEL_X + 66}" y="{y + 6}" font-family="{FONT}" font-size="16">{escape(text)}</text>'
        )
    # subworkflow boxes
    for y0, y1, (name, purpose, _steps) in boxes:
        out.append(
            f'<rect x="30" y="{y0}" width="{WIDTH - 60}" height="{y1 - y0}" rx="12" fill="#f7f7f7" stroke="#9a9a9a" stroke-width="2" stroke-dasharray="7 5"/>'
        )
        out.append(
            f'<text x="{LABEL_X}" y="{y0 + 29}" font-family="{MONO}" font-size="17" font-weight="bold" fill="#333">{escape(name)}</text>'
        )
        out.append(
            f'<text x="{WIDTH - 50}" y="{y0 + 29}" font-family="{FONT}" font-size="16" font-style="italic" fill="#555" text-anchor="end">{escape(purpose)}</text>'
        )
    # reference feeds (dashed), below the lines
    ref_colour = TRACKS["ref"][0]
    y_ref_end = steps["2.2"][0]
    for number in FEEDS:
        y, step = steps[number]
        x = TRACK_X[step[3]]
        out.append(
            f'<path d="M {TRACK_X["ref"]} {y_ref_end} L {TRACK_X["ref"]} {y - ROW * 0.6} L {x} {y}" fill="none" '
            f'stroke="{ref_colour}" stroke-width="4" stroke-dasharray="8 6" stroke-linejoin="round"/>'
        )
    # lines
    for track, path in LINES:
        points = path_points(track, path, steps)
        d = " ".join(f"{'M' if i == 0 else 'L'} {x} {y}" for i, (x, y) in enumerate(points))
        out.append(
            f'<path d="{d}" fill="none" stroke="{TRACKS[track][0]}" stroke-width="10" stroke-linecap="round" stroke-linejoin="round"/>'
        )
    # file icons at the line starts
    for track, (colour, icon, _text) in TRACKS.items():
        out.append(file_icon(TRACK_X[track], ICON_Y, icon, colour))
    # stations and labels
    for number, (y, (_n, text, tools, track)) in steps.items():
        x = TRACK_X[track]
        out.append(f'<circle cx="{x}" cy="{y}" r="11" fill="#fff" stroke="#000" stroke-width="5"/>')
        out.append(
            f'<line x1="{x + 16}" y1="{y}" x2="{LABEL_X - 10}" y2="{y}" stroke="#d0d0d0" stroke-width="1.5" stroke-dasharray="2 4"/>'
        )
        out.append(
            f'<text x="{LABEL_X}" y="{y + 6}" font-family="{FONT}" font-size="17"><tspan font-weight="bold">{number}</tspan>  {escape(text)}</text>'
        )
        out.append(
            f'<text x="{TOOLS_X}" y="{y + 6}" font-family="{MONO}" font-size="14" fill="#555">{escape(tools)}</text>'
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
    _boxes, steps, _height = layout()
    check_data(steps)
    names = {key: text.split(":")[0] for key, (_c, _i, text) in TRACKS.items()}
    rows = ["| Step | What happens | Tools | Line |", "| --- | --- | --- | --- |"]
    for name, purpose, band_steps in BANDS:
        rows.append(f"| **{name}** | _{purpose}_ | | |")
        for number, text, tools, track in band_steps:
            rows.append(
                f"| {number} | **{text}.** {DETAILS[number]} | `{tools}` | {names[track]} |"
            )
    return "\n".join(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--outdir", default="docs/images", type=Path, help="output directory")
    parser.add_argument("--markdown", action="store_true", help="print the README step table")
    args = parser.parse_args()
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
