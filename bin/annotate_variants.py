#!/usr/bin/env python3
"""Annotate the compared variants: VEP consequences, repeats and stratifications.

The annotation itself is done by tools (Ensembl VEP, msisensor-pro, bedtools); this
script converts the variant table into their input formats and merges their results
back into the table.

Subcommands
  to-vcf          variant table -> sorted sites-only VCF (bgzipped + .tbi) for VEP
      INPUT   --table: TSV with at least chrom, pos (1-based), ref, alt;
              --fai: reference index (contig order and lengths of the VCF header)
      OUTPUT  --output: VCF, one record per distinct CHROM/POS/REF/ALT
  to-bed          variant table -> BED4 of the REF span of every row, for bedtools
      OUTPUT  chrom, start (pos - 1), end (pos - 1 + len(ref)), row number (1-based)
  repeats-to-bed  msisensor-pro scan list -> BED of the repeat (microsatellite) regions
      INPUT   tab-separated with header: chromosome, location (0-based start),
              repeat_unit_length, repeat_times, repeat_unit_bases, ... (msisensor-pro 1.3)
      OUTPUT  chrom, start, end (start + unit length x times), repeat unit, repeat times
  merge           variant table + tool results -> annotated table (all input columns
                  plus the new ones below)
      --vep       TSV of `bcftools +split-vep -H -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t
                  %Consequence\\t%SYMBOL\\t%IMPACT\\t%VARIANT_CLASS\\t%BIOTYPE\\n' -s worst`
                  -> columns consequence (most severe SO term of the variant, '&'-joined
                  if VEP gives several), consequence_group (synonymous, non_synonymous,
                  other; see CONSEQUENCE_GROUPS), gene (SYMBOL, '.' if none), impact
                  (HIGH, MODERATE, LOW, MODIFIER), variant_class (SO class: SNV,
                  insertion, deletion, ...), biotype
                  Optional extra fields after these nine (any subset, in this order;
                  requested by conf/modules/annotation.config only when VEP writes them):
                  %Feature -> transcript; %HGVSc, %HGVSp (VEP --hgvs) -> hgvsc, hgvsp
                  without the transcript/protein ID prefix (c.35G>T, p.Gly12Val);
                  %Existing_variation (VEP cache, --check_existing) -> dbsnp (rs IDs),
                  cosmic (COSV/COSM IDs), comma-separated; %CLIN_SIG -> clin_sig
                  (comma-separated). These six columns are always written ('.' if the
                  field was not given or is empty), so the table has the same columns
                  with a VEP cache and with a GTF.
      --repeats   output of `bedtools map -c 4,5 -o distinct,max` of the variant BED
                  with the repeat BED -> columns in_repeat (yes/no), repeat_unit
                  (repeat unit bases, '.' if none), repeat_times
      --stratum   NAME=FILE, output of `bedtools intersect -c` of the variant BED with
                  a stratification BED -> column NAME = yes (> 0 intervals) / no

Every variant (row number, resp. CHROM/POS/REF/ALT) must appear in every tool result,
otherwise the script stops: a missing row would silently lose an annotation.
"""

import argparse
import csv
import re
import sys

REQUIRED = ["chrom", "pos", "ref", "alt"]
NAME_RE = re.compile(r"^[A-Za-z0-9_]+$")
VEP_COLUMNS = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "Consequence",
    "SYMBOL",
    "IMPACT",
    "VARIANT_CLASS",
    "BIOTYPE",
]
# optional split-vep fields after VEP_COLUMNS (in this order)
OPTIONAL_VEP_COLUMNS = ["Feature", "HGVSc", "HGVSp", "Existing_variation", "CLIN_SIG"]
# columns written from the optional fields
OPTIONAL_COLUMNS = ["transcript", "hgvsc", "hgvsp", "dbsnp", "cosmic", "clin_sig"]

# Sequence Ontology consequence terms (Ensembl VEP) -> group of the report
CONSEQUENCE_GROUPS = {
    "synonymous": {"synonymous_variant", "stop_retained_variant", "start_retained_variant"},
    "non_synonymous": {
        "missense_variant",
        "stop_gained",
        "stop_lost",
        "start_lost",
        "frameshift_variant",
        "inframe_insertion",
        "inframe_deletion",
        "protein_altering_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "transcript_ablation",
        "incomplete_terminal_codon_variant",
    },
}


def read_table(path):
    """Header and rows of the variant table; checks the required columns."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        missing = [c for c in REQUIRED if c not in header]
        if missing:
            sys.exit(
                f"ERROR: {path}: missing column(s) {', '.join(missing)} "
                f"(header: {', '.join(header)})"
            )
        rows = list(reader)
    # positions are positive integers, alleles are not empty
    for i, row in enumerate(rows, start=2):
        if not row["pos"].isdigit():
            sys.exit(f"ERROR: {path}: line {i}: pos '{row['pos']}' is not a positive integer")
        if not row["ref"] or not row["alt"]:
            sys.exit(f"ERROR: {path}: line {i}: empty ref or alt")
    return header, rows


def read_fai(path):
    """{contig: length} of a FASTA index, in reference order."""
    lengths = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            lengths[fields[0]] = int(fields[1])
    return lengths


def consequence_group(consequence):
    """Group of a VEP consequence ('&'-joined terms): the most protein-relevant group."""
    terms = set(consequence.split("&")) if consequence not in ("", ".") else set()
    if terms & CONSEQUENCE_GROUPS["non_synonymous"]:
        return "non_synonymous"
    if terms & CONSEQUENCE_GROUPS["synonymous"]:
        return "synonymous"
    return "other"


# ------------------------------------------------------------------ subcommands


def to_vcf(args):
    """Sorted sites-only VCF of the distinct variants of the table (VEP input)."""
    import pysam

    _header, rows = read_table(args.table)
    reference = read_fai(args.fai)
    order = {contig: i for i, contig in enumerate(reference)}
    # distinct variants; every contig must be a reference contig
    variants = sorted(
        {(r["chrom"], int(r["pos"]), r["ref"], r["alt"]) for r in rows},
        key=lambda v: (order.get(v[0], len(order)), v[1], v[2], v[3]),
    )
    unknown = sorted({v[0] for v in variants if v[0] not in reference})
    if unknown:
        sys.exit(f"ERROR: {args.table}: contig(s) {', '.join(unknown)} not in {args.fai}")
    header = pysam.VariantHeader()
    for name, length in reference.items():
        header.contigs.add(name, length=length)
    with pysam.VariantFile(args.output, "wz", header=header) as out:
        for chrom, pos, ref, alt in variants:
            out.write(out.new_record(contig=chrom, start=pos - 1, alleles=(ref, alt)))
    pysam.tabix_index(args.output, preset="vcf", force=True)
    print(f"INFO: {args.table}: {len(variants)} distinct variants written", file=sys.stderr)


def to_bed(args):
    """BED4 of the REF span of every table row; the name is the row number."""
    _header, rows = read_table(args.table)
    with open(args.output, "w") as fh:
        for i, row in enumerate(rows, start=1):
            # positions are 1-based in the table, BED is 0-based half-open
            start = int(row["pos"]) - 1
            fh.write(f"{row['chrom']}\t{start}\t{start + len(row['ref'])}\t{i}\n")
    print(f"INFO: {args.table}: {len(rows)} variants written as BED", file=sys.stderr)


def repeats_to_bed(args):
    """BED (chrom, start, end, unit, times) of the msisensor-pro scan list."""
    needed = ["chromosome", "location", "repeat_unit_length", "repeat_times", "repeat_unit_bases"]
    n = 0
    with open(args.list, newline="") as fin, open(args.output, "w") as out:
        reader = csv.DictReader(fin, delimiter="\t")
        missing = [c for c in needed if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(
                f"ERROR: {args.list}: missing column(s) {', '.join(missing)}; expected the "
                "output of msisensor-pro scan (version 1.3)"
            )
        for line_no, row in enumerate(reader, start=2):
            try:
                start = int(row["location"])
                length = int(row["repeat_unit_length"]) * int(row["repeat_times"])
            except ValueError:
                sys.exit(f"ERROR: {args.list}: line {line_no}: location/length are not integers")
            # the repeat covers unit length x times bases from the 0-based location
            out.write(
                f"{row['chromosome']}\t{start}\t{start + length}\t"
                f"{row['repeat_unit_bases']}\t{row['repeat_times']}\n"
            )
            n += 1
    print(f"INFO: {args.list}: {n} repeat regions written", file=sys.stderr)


def read_bedtools_result(path, n_rows, what, n_values=1):
    """{row number: [values of columns 5, 6, ...]} of a bedtools map/intersect result."""
    values = {}
    with open(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 + n_values:
                sys.exit(
                    f"ERROR: {path} ({what}): line {line_no}: expected {4 + n_values} columns "
                    f"(BED4 + result), found {len(fields)}"
                )
            row = int(fields[3])
            if row in values:
                sys.exit(f"ERROR: {path} ({what}): row {row} occurs twice")
            values[row] = fields[4 : 4 + n_values]
    # every variant must be annotated exactly once
    missing = [r for r in range(1, n_rows + 1) if r not in values]
    if missing:
        sys.exit(
            f"ERROR: {path} ({what}): {len(missing)} variant(s) missing (rows "
            f"{', '.join(map(str, missing[:5]))}{', ...' if len(missing) > 5 else ''})"
        )
    return values


def read_vep(path):
    """{(chrom, pos, ref, alt): VEP fields} of the split-vep TSV; checks the columns."""
    with open(path, newline="") as fh:
        # split-vep -H writes a header like '# [1]CHROM [2]POS ...'
        first = fh.readline().lstrip("#").strip().split("\t")
        names = [re.sub(r"^\[\d+\]", "", name) for name in first]
        # the nine fixed columns first, then optional fields in the order of
        # OPTIONAL_VEP_COLUMNS (a subsequence of it)
        extra = names[len(VEP_COLUMNS) :]
        in_order = [c for c in OPTIONAL_VEP_COLUMNS if c in extra]
        if names[: len(VEP_COLUMNS)] != VEP_COLUMNS or extra != in_order:
            sys.exit(
                f"ERROR: {path}: columns {', '.join(names)}; expected {', '.join(VEP_COLUMNS)} "
                f"followed by any of {', '.join(OPTIONAL_VEP_COLUMNS)} in this order "
                "(bcftools +split-vep -H with the format of conf/modules/annotation.config)"
            )
        annotations = {}
        for line in fh:
            fields = dict(zip(names, line.rstrip("\n").split("\t"), strict=True))
            key = (fields["CHROM"], int(fields["POS"]), fields["REF"], fields["ALT"])
            annotations[key] = fields
    return annotations


def strip_id(hgvs):
    """HGVS notation without the transcript/protein ID ('ENST1.2:c.35G>T' -> 'c.35G>T');
    '.' for an empty value."""
    if not hgvs or hgvs == ".":
        return "."
    return hgvs.split(":", 1)[1] if ":" in hgvs else hgvs


def optional_vep_values(fields):
    """The OPTIONAL_COLUMNS of one split-vep line ('.' for fields not given or empty)."""
    # Existing_variation: '&'-joined IDs of known variants (rs..., COSV..., COSM..., CM...)
    known = [v for v in fields.get("Existing_variation", "").split("&") if v not in ("", ".")]
    dbsnp = [v for v in known if v.startswith("rs")]
    cosmic = [v for v in known if v.startswith(("COSV", "COSM"))]
    # CLIN_SIG: '&'-joined ClinVar significances
    clin_sig = [v for v in fields.get("CLIN_SIG", "").split("&") if v not in ("", ".")]
    return {
        "transcript": fields.get("Feature") or ".",
        "hgvsc": strip_id(fields.get("HGVSc")),
        "hgvsp": strip_id(fields.get("HGVSp")),
        "dbsnp": ",".join(dbsnp) or ".",
        "cosmic": ",".join(cosmic) or ".",
        "clin_sig": ",".join(clin_sig) or ".",
    }


def merge(args):
    """Annotated table: input columns + VEP, repeat and stratification columns."""
    header, rows = read_table(args.table)
    n = len(rows)
    new_columns = []
    # VEP: one annotation per distinct variant (most severe consequence)
    if args.vep:
        vep = read_vep(args.vep)
        for i, row in enumerate(rows, start=2):
            key = (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
            if key not in vep:
                sys.exit(f"ERROR: {args.vep}: no VEP annotation for line {i} of the table: {key}")
            fields = vep[key]
            row["consequence"] = fields["Consequence"] or "."
            row["consequence_group"] = consequence_group(fields["Consequence"])
            row["gene"] = fields["SYMBOL"] or "."
            row["impact"] = fields["IMPACT"] or "."
            row["variant_class"] = fields["VARIANT_CLASS"] or "."
            row["biotype"] = fields["BIOTYPE"] or "."
            row.update(optional_vep_values(fields))
        new_columns += [
            "consequence",
            "consequence_group",
            "gene",
            "impact",
            "variant_class",
            "biotype",
            *OPTIONAL_COLUMNS,
        ]
        n_hgvsp = sum(r["hgvsp"] != "." for r in rows)
        n_dbsnp = sum(r["dbsnp"] != "." for r in rows)
        print(
            f"INFO: VEP: {n_hgvsp} of {n} variants with HGVSp, {n_dbsnp} with a dbSNP ID",
            file=sys.stderr,
        )
        groups = {
            g: sum(r["consequence_group"] == g for r in rows)
            for g in ("synonymous", "non_synonymous", "other")
        }
        print(f"INFO: VEP consequence groups: {groups}", file=sys.stderr)
    # repeats: bedtools map writes '.' if no repeat overlaps
    if args.repeats:
        repeats = read_bedtools_result(args.repeats, n, "repeats", n_values=2)
        for i, row in enumerate(rows, start=1):
            unit, times = repeats[i]
            row["in_repeat"] = "no" if unit == "." else "yes"
            row["repeat_unit"] = unit
            row["repeat_times"] = times
        new_columns += ["in_repeat", "repeat_unit", "repeat_times"]
        n_in = sum(r["in_repeat"] == "yes" for r in rows)
        print(f"INFO: repeats: {n_in} of {n} variants in a repeat region", file=sys.stderr)
    # stratifications: bedtools intersect -c writes the number of overlapping intervals
    for item in args.stratum:
        if "=" not in item:
            sys.exit(f"ERROR: --stratum '{item}' must be NAME=FILE")
        name, path = item.split("=", 1)
        if not NAME_RE.match(name) or name in header or name in new_columns:
            sys.exit(f"ERROR: stratum name '{name}' is invalid or already a column")
        counts = read_bedtools_result(path, n, f"stratum {name}")
        for i, row in enumerate(rows, start=1):
            try:
                row[name] = "yes" if int(counts[i][0]) > 0 else "no"
            except ValueError:
                sys.exit(
                    f"ERROR: {path} (stratum {name}): row {i}: count '{counts[i][0]}' is not an "
                    "integer; was the file made with bedtools intersect -c?"
                )
        new_columns.append(name)
        n_yes = sum(1 for row in rows if row[name] == "yes")
        print(f"INFO: stratum {name}: {n_yes} of {n} variants inside", file=sys.stderr)
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=header + new_columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    """Parse the command line and run the (sub)command."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    v = sub.add_parser("to-vcf", help="variant table -> sites-only VCF (VEP input)")
    v.add_argument("--table", required=True)
    v.add_argument("--fai", required=True, help="reference index (contig order)")
    v.add_argument("--output", required=True, help="output .vcf.gz")
    b = sub.add_parser("to-bed", help="variant table -> BED4 (row number as name)")
    b.add_argument("--table", required=True)
    b.add_argument("--output", required=True)
    r = sub.add_parser("repeats-to-bed", help="msisensor-pro scan list -> BED")
    r.add_argument("--list", required=True, help="msisensor-pro scan output")
    r.add_argument("--output", required=True)
    m = sub.add_parser("merge", help="merge the tool results into the variant table")
    m.add_argument("--table", required=True)
    m.add_argument("--vep", help="bcftools +split-vep TSV")
    m.add_argument("--repeats", help="bedtools map result with the repeat BED")
    m.add_argument("--stratum", nargs="*", default=[], help="NAME=bedtools intersect -c result")
    m.add_argument("--output", required=True)
    args = parser.parse_args()
    {"to-vcf": to_vcf, "to-bed": to_bed, "repeats-to-bed": repeats_to_bed, "merge": merge}[
        args.command
    ](args)


if __name__ == "__main__":
    main()
