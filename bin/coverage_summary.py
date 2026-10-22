#!/usr/bin/env python3
"""Coverage of every alignment in the common comparison region and in its target region.

For each alignment two regions are summarized:
  common   the common comparison region (mosdepth --by <common region> --thresholds)
  target   the target BED of its technology (second mosdepth run, --by <target BED>
           --thresholds); for a technology without target BED (WGS) the whole genome
           (mosdepth global distribution and summary)

INPUT
  --alignments   TSV with header: alignment (id <entry>_<tumor|normal>), sample, technology,
                 role
  --dir          folder with the mosdepth outputs, named by the alignment id <ID>:
                   <ID>.regions.bed.gz, <ID>.thresholds.bed.gz    common region: mean depth
                       per interval (contig, start, end, [name,] mean) and bases >= each
                       threshold (header '#chrom start end region 10X ...')
                   <ID>.target.regions.bed.gz, <ID>.target.thresholds.bed.gz   the same for
                       the target BED (targeted technologies only)
                   <ID>.mosdepth.global.dist.txt, <ID>.mosdepth.summary.txt     genome
                       (contig, depth, fraction of bases >= depth; 'total' summary row),
                       used when there is no target BED
  --thresholds   the depths of mosdepth --thresholds (e.g. 10,20,30,50,100)

OUTPUT (--prefix)
  <prefix>.coverage.tsv          one row per alignment and region: alignment, sample,
                                 technology, role (tumor/normal), region (common, target,
                                 genome), bp (bases of the region), mean_depth (bases-weighted
                                 mean), pct_ge_<t>X (percent of the bases with depth >= t, one
                                 column per threshold)
  <prefix>_generalstats_mqc.tsv  MultiQC General Statistics columns per alignment: mean depth
                                 and % >= each threshold in the common region and in the
                                 target region (genome for WGS)

Checks: every alignment has a common-region summary; mosdepth threshold columns match
--thresholds; a region without bases is an error.

Why a script: mosdepth writes the per-interval values; the bases-weighted region means and
the General Statistics table are bookkeeping across files.
"""

import argparse
import gzip
import sys
from pathlib import Path


def read_alignments(path):
    """[(alignment, sample, technology, role)] of the alignment table."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header[:4] != ["alignment", "sample", "technology", "role"]:
            sys.exit(f"ERROR: {path}: expected the header alignment, sample, technology, role")
        return [tuple(line.rstrip("\n").split("\t")[:4]) for line in fh if line.strip()]


def open_text(path):
    """Text handle of a plain or gzipped file."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def region_coverage(regions_path, thresholds_path, thresholds):
    """(bp, mean depth, {threshold: percent of bases >= threshold}) of a region from the
    mosdepth regions and thresholds files of one run."""
    bp, depth_sum = 0, 0.0
    with open_text(regions_path) as fh:
        for line_no, line in enumerate(fh, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                sys.exit(
                    f"ERROR: {regions_path}: line {line_no}: expected contig, start, end, mean"
                )
            length = int(fields[2]) - int(fields[1])
            # the mean depth is the last column (a name column may come before it)
            bp += length
            depth_sum += float(fields[-1]) * length
    covered = dict.fromkeys(thresholds, 0)
    with open_text(thresholds_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        expected = [f"{t}X" for t in thresholds]
        if header[4:] != expected:
            sys.exit(
                f"ERROR: {thresholds_path}: threshold columns {header[4:]}, expected {expected} "
                "(mosdepth --thresholds must match --thresholds)"
            )
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            for t, value in zip(thresholds, fields[4:], strict=True):
                covered[t] += int(value)
    if bp == 0:
        sys.exit(f"ERROR: {regions_path}: the region has no bases")
    return bp, depth_sum / bp, {t: 100.0 * covered[t] / bp for t in thresholds}


def genome_coverage(dist_path, summary_path, thresholds):
    """(bp, mean depth, {threshold: percent}) of the whole genome from mosdepth's global
    distribution ('total' rows) and summary ('total' row)."""
    fraction = {}
    with open(dist_path) as fh:
        for line in fh:
            contig, depth, value = line.rstrip("\n").split("\t")
            if contig == "total":
                fraction[int(depth)] = float(value)
    bp, mean = None, None
    with open(summary_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "total":
                bp, mean = int(fields[1]), float(fields[3])
    if bp is None:
        sys.exit(f"ERROR: {summary_path}: no 'total' row")
    # the distribution lists every depth up to the maximum; beyond it the fraction is 0
    return bp, mean, {t: 100.0 * fraction.get(t, 0.0) for t in thresholds}


def main():
    """Parse the command line and summarize the coverage of every alignment."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--alignments", required=True, help="alignment table (TSV)")
    parser.add_argument("--dir", default=".", help="folder of the mosdepth outputs")
    parser.add_argument("--thresholds", required=True, help="e.g. 10,20,30,50,100")
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()
    thresholds = [int(t) for t in args.thresholds.split(",")]
    folder = Path(args.dir)
    alignments = read_alignments(args.alignments)

    rows = []
    for aid, sample, technology, role in sorted(alignments):
        common = (folder / f"{aid}.regions.bed.gz", folder / f"{aid}.thresholds.bed.gz")
        if not all(p.exists() for p in common):
            sys.exit(f"ERROR: alignment {aid}: no mosdepth output of the common region in {folder}")
        regions = [("common", region_coverage(*common, thresholds))]
        # the target region of a targeted technology, else the genome
        target = (
            folder / f"{aid}.target.regions.bed.gz",
            folder / f"{aid}.target.thresholds.bed.gz",
        )
        genome = (
            folder / f"{aid}.mosdepth.global.dist.txt",
            folder / f"{aid}.mosdepth.summary.txt",
        )
        if all(p.exists() for p in target):
            regions.append(("target", region_coverage(*target, thresholds)))
        elif all(p.exists() for p in genome):
            regions.append(("genome", genome_coverage(*genome, thresholds)))
        for region, (bp, mean, pct) in regions:
            rows.append(
                [
                    aid,
                    sample,
                    technology,
                    role,
                    region,
                    str(bp),
                    f"{mean:.2f}",
                    *(f"{pct[t]:.2f}" for t in thresholds),
                ]
            )

    columns = [
        "alignment",
        "sample",
        "technology",
        "role",
        "region",
        "bp",
        "mean_depth",
        *(f"pct_ge_{t}X" for t in thresholds),
    ]
    with open(f"{args.prefix}.coverage.tsv", "w") as fh:
        fh.write("\t".join(columns) + "\n")
        fh.writelines("\t".join(row) + "\n" for row in rows)

    # MultiQC General Statistics: one row per alignment, columns per region kind
    stats = {}
    for row in rows:
        kind = "common" if row[4] == "common" else "target"
        values = stats.setdefault(row[0], {})
        values[f"{kind}_mean"] = row[6]
        for t, value in zip(thresholds, row[7:], strict=True):
            values[f"{kind}_ge_{t}X"] = value
    keys = [
        f"{kind}_{name}"
        for kind in ("common", "target")
        for name in ["mean", *(f"ge_{t}X" for t in thresholds)]
    ]
    with open(f"{args.prefix}_generalstats_mqc.tsv", "w") as fh:
        fh.write("# id: 'variantcomp_coverage'\n# plot_type: 'generalstats'\n# pconfig:\n")
        for key in keys:
            kind, name = key.split("_", 1)
            region = "common region" if kind == "common" else "target region (genome for WGS)"
            if name == "mean":
                fh.write(
                    f"#     {key}:\n#         title: 'Depth {kind}'\n#         description: 'Mean depth in the {region}'\n#         format: '{{:,.1f}}'\n"
                )
            else:
                fh.write(
                    f"#     {key}:\n#         title: '{name.removeprefix('ge_')} {kind}'\n#         description: '% bases with depth >= {name.removeprefix('ge_')} in the {region}'\n#         max: 100\n#         suffix: '%'\n#         hidden: {str(name != f'ge_{thresholds[min(2, len(thresholds) - 1)]}X').lower()}\n"
                )
        fh.write("Sample\t" + "\t".join(keys) + "\n")
        for aid in sorted(stats):
            fh.write(aid + "\t" + "\t".join(stats[aid].get(k, "NA") for k in keys) + "\n")
    print(
        f"INFO: coverage of {len(alignments)} alignments summarized ({len(rows)} rows)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
