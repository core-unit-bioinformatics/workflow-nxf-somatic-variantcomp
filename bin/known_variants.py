#!/usr/bin/env python3
"""Known variants: preparation and the 'required' check.

prepare  Read known variants from one or more VCFs or TSVs (columns chrom, pos,
         ref, alt and optionally gene, label, vaf; 1-based positions), map contig
         names to the reference (chr prefix, MT/chrM), check positions and alleles
         and write a sorted, sites-only VCF (INFO/KV_GENE, INFO/KV_LABEL
         and, if a VAF is known, INFO/VC_T_VAF, so recall can be stratified by the
         truth VAF). VCF VAF tags: --vaf-tag, else the first of TVAF, VAF, AF.
         Normalize the output with bcftools norm afterwards (left-align, split; the
         REF alleles are checked there against the FASTA: --check-ref x, and the
         pipeline stops if any REF allele differs, see check_reference.py ref-alleles).

check    For the 'required' mode: every known variant must be found in every
         technology of a sample. For each technology the variant is looked up in
         the normalized VCF (all calls), the compared VCF (calls passing the
         uniform filters in the comparison region), the comparison region and,
         if available, the tumor read evidence. Status per variant:

           PASS                 in the compared calls of the technology
           FAIL_filtered        called, but failed FILTER or the uniform thresholds
           FAIL_outside_region  not in the comparison (common/callable) region
           FAIL_not_called      in the region, but not called
                                (read evidence tells whether reads support it)

         Matching is by CHROM/POS/REF/ALT after normalization; a call at the same
         position with another allele is reported in `same_position_call`.

Outputs of check (prefix):
  <prefix>.known_variants.tsv   one row per technology and known variant
  <prefix>.known_summary.tsv    per technology: required, passed, all_passed
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_reference import map_contig, read_fai  # noqa: E402

TSV_REQUIRED = ["chrom", "pos", "ref", "alt"]


# a 'required' list is a short clinical checklist (user requirement: fewer than 10)
MAX_REQUIRED = 10


def open_text(path):
    """Text handle of a plain or gzipped file."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def first_value(value):
    """First element of a tuple-like INFO value (Number=A/R/.), or the value itself."""
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def read_known(path, vaf_tag=None):
    """Yield (chrom, pos, ref, alt, gene, label, vaf) from a VCF or TSV."""
    if str(path).endswith((".vcf", ".vcf.gz", ".bcf")):
        with pysam.VariantFile(str(path)) as vcf:
            # VAF tag: --vaf-tag, else the first of TVAF, VAF, AF defined in the header
            tags = [vaf_tag] if vaf_tag else ["TVAF", "VAF", "AF"]
            tag = next((t for t in tags if t in vcf.header.info), None)
            if vaf_tag and tag is None:
                sys.exit(f"ERROR: {path}: --vaf-tag {vaf_tag} is not defined in the VCF header")
            print(
                (
                    f"INFO: {path}: VAF from INFO/{tag}"
                    if tag
                    else f"INFO: {path}: no VAF tag "
                    "(TVAF, VAF, AF) in the header; VAF-based strata will be 'unknown'"
                ),
                file=sys.stderr,
            )
            # gene and label: INFO/GENE or KV_GENE; the ID column, else INFO/KV_LABEL
            gene_tag = next((t for t in ("GENE", "KV_GENE") if t in vcf.header.info), None)
            label_tag = "KV_LABEL" if "KV_LABEL" in vcf.header.info else None
            for rec in vcf:
                # one known variant per ALT allele
                gene = (rec.info.get(gene_tag) if gene_tag else None) or "."
                label = rec.id or (rec.info.get(label_tag) if label_tag else None) or "."
                vaf = first_value(rec.info.get(tag)) if tag else None
                for alt in rec.alts or ():
                    yield rec.chrom, rec.pos, rec.ref, alt, str(gene), str(label), vaf
        return
    with open_text(path) as fh:
        # CSV or TSV list, by file extension
        delimiter = "\t" if str(path).endswith((".tsv", ".tsv.gz")) else ","
        reader = csv.DictReader(
            (line for line in fh if not line.startswith("##")), delimiter=delimiter
        )
        header = [h.lower().lstrip("#") for h in reader.fieldnames or []]
        missing = [c for c in TSV_REQUIRED if c not in header]
        if missing:
            sys.exit(f"ERROR: {path}: missing column(s) {', '.join(missing)}")
        # column names are case-insensitive (CHROM, Chrom, chrom)
        reader.fieldnames = header
        # line numbers count the header line (line 1 = header)
        for line_no, row in enumerate(reader, start=2):
            where = f"{path}: line {line_no}"
            # 1-based position: a positive integer
            if not (row["pos"] or "").strip().isdigit():
                sys.exit(f"ERROR: {where}: pos '{row['pos']}' is not a positive integer")
            # optional VAF: a number within [0, 1]; NA, '.' or empty = unknown
            vaf = (row.get("vaf") or "").strip()
            if vaf in ("", "NA", "."):
                vaf = None
            else:
                try:
                    vaf = float(vaf)
                except ValueError:
                    sys.exit(f"ERROR: {where}: vaf '{vaf}' is not a number")
                if not 0 <= vaf <= 1:
                    sys.exit(f"ERROR: {where}: vaf {vaf} is outside [0, 1]")
            yield (
                row["chrom"].strip(),
                int(row["pos"]),
                row["ref"].strip().upper(),
                row["alt"].strip().upper(),
                row.get("gene") or ".",
                row.get("label") or ".",
                vaf,
            )


def prepare(args):
    # reference contigs and their order (for sorting the output)
    reference = read_fai(args.fai)
    order = {c: i for i, c in enumerate(reference)}
    records, errors = [], []
    # all variants of all input files of the list
    known = [v for path in args.input for v in read_known(path, args.vaf_tag)]
    if not known:
        sys.exit(f"ERROR: {', '.join(args.input)}: the known-variant list is empty")
    for chrom, pos, ref, alt, gene, label, vaf in known:
        # contig name on the reference (chr prefix, MT/chrM); unknown contigs are an error
        ref_chrom = map_contig(chrom, reference)
        if ref_chrom is None:
            errors.append(
                f"{chrom}:{pos} {ref}>{alt} ({label}): contig '{chrom}' is not in the reference"
            )
            continue
        # position must lie on the contig (the REF alleles are checked by bcftools norm afterwards)
        if pos < 1 or pos - 1 + len(ref) > reference[ref_chrom]:
            errors.append(
                f"{chrom}:{pos} {ref}>{alt} ({label}): position beyond the end of {ref_chrom} "
                f"(length {reference[ref_chrom]}): another genome build?"
            )
            continue
        # alleles must be plain bases (symbolic alleles can't be checked against a call)
        if not set(ref + alt) <= set("ACGTN") or ref == alt:
            errors.append(
                f"{chrom}:{pos} {ref}>{alt} ({label}): REF/ALT must be different A/C/G/T/N alleles"
            )
            continue
        records.append((ref_chrom, pos, ref, alt, gene, label, vaf))
    if errors:
        sys.exit("ERROR: known variants not usable:\n  " + "\n  ".join(errors))
    # reference order, then position and alleles (needed for tabix and bcftools)
    records.sort(key=lambda r: (order[r[0]], r[1], r[2], r[3]))

    # sites-only VCF header: reference contigs and the KV_* / VC_T_VAF INFO tags
    header = pysam.VariantHeader()
    for name, length in reference.items():
        header.contigs.add(name, length=length)
    header.info.add("KV_GENE", 1, "String", "Gene of the known variant")
    header.info.add("KV_LABEL", 1, "String", "Label of the known variant")
    header.info.add("VC_T_VAF", 1, "Float", "Tumor VAF of the known variant (if given)")
    with pysam.VariantFile(args.output, "wz", header=header) as out:
        seen = set()
        for chrom, pos, ref, alt, gene, label, vaf in records:
            # the same variant listed twice (e.g. in two input files) is written once
            if (chrom, pos, ref, alt) in seen:
                continue
            seen.add((chrom, pos, ref, alt))
            info = {"KV_GENE": gene.replace(" ", "_"), "KV_LABEL": label.replace(" ", "_")}
            if vaf is not None:
                info["VC_T_VAF"] = float(vaf)
            rec = out.new_record(
                contig=chrom, start=pos - 1, alleles=(ref, alt), id=label, info=info
            )
            out.write(rec)
    pysam.tabix_index(args.output, preset="vcf", force=True)
    n_vaf = sum(1 for r in records if r[6] is not None)
    print(
        f"INFO: {', '.join(args.input)}: {len(seen)} known variants written "
        f"({len(records) - len(seen)} duplicates removed, {n_vaf} with a VAF)",
        file=sys.stderr,
    )


def load_calls(path):
    """{(chrom, pos, ref, alt): FILTER, VAF, depth, alt reads} of a harmonized VCF."""
    calls = {}
    with pysam.VariantFile(path) as vcf:
        # the VC_* tags come from harmonize_vcf.py
        if "VC_T_VAF" not in vcf.header.info:
            sys.exit(
                f"ERROR: {path}: INFO/VC_T_VAF is not defined; expected a VCF harmonized by "
                "harmonize_vcf.py"
            )
        for rec in vcf:
            for alt in rec.alts or ():
                calls[(rec.chrom, rec.pos, rec.ref, alt)] = {
                    "filter": ";".join(rec.filter.keys()) or ".",
                    "vaf": rec.info.get("VC_T_VAF"),
                    "dp": rec.info.get("VC_T_DP"),
                    "alt": rec.info.get("VC_T_ALT"),
                }
    return calls


def load_regions(path):
    """{chrom: [(start, end), ...]} of a BED; header and empty lines are skipped."""
    regions = {}
    with open_text(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or not fields[1].isdigit() or not fields[2].isdigit():
                sys.exit(f"ERROR: {path}: line {line_no}: expected chrom, start, end (integers)")
            regions.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    return regions


def in_regions(regions, chrom, pos, ref):
    """True if the REF span of a variant overlaps an interval (0-based half-open)."""
    start, end = pos - 1, pos - 1 + len(ref)
    return any(s < end and start < e for s, e in regions.get(chrom, ()))


def load_evidence(path):
    """Pileup evidence of one alignment: ({(chrom, pos, ref, alt): alt reads}, {(chrom, pos): depth})."""
    if not path:
        return None
    alts, depth = {}, {}
    with pysam.VariantFile(path) as vcf:
        # one alignment = one sample with FORMAT/AD (bcftools mpileup -a AD,DP)
        if len(vcf.header.samples) != 1 or "AD" not in vcf.header.formats:
            sys.exit(
                f"ERROR: {path}: expected one sample with FORMAT/AD (bcftools mpileup -a AD,DP)"
            )
        sample = next(iter(vcf.header.samples))
        for rec in vcf:
            fmt = rec.samples[sample]
            ad = fmt.get("AD") if "AD" in rec.format else None
            dp = fmt.get("DP") if "DP" in rec.format else None
            # several records at one position: keep the largest depth
            if dp is not None:
                depth[(rec.chrom, rec.pos)] = max(depth.get((rec.chrom, rec.pos), 0), dp)
            for i, alt in enumerate(rec.alts or (), start=1):
                if ad is not None and i < len(ad) and ad[i] is not None:
                    alts[(rec.chrom, rec.pos, rec.ref, alt)] = ad[i]
    return alts, depth


def fmt(value):
    """Table value: NA for missing, 4 digits for floats."""
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def check(args):
    with pysam.VariantFile(args.known) as vcf:
        # the known-variant VCF is the output of 'prepare' (+ bcftools norm)
        missing = [t for t in ("KV_GENE", "KV_LABEL") if t not in vcf.header.info]
        if missing:
            sys.exit(
                f"ERROR: {args.known}: INFO/{', INFO/'.join(missing)} not defined; expected "
                "the output of 'known_variants.py prepare'"
            )
        known = [
            (
                rec.chrom,
                rec.pos,
                rec.ref,
                alt,
                rec.info.get("KV_GENE", "."),
                rec.info.get("KV_LABEL", "."),
            )
            for rec in vcf
            for alt in rec.alts or ()
        ]
    if not known:
        sys.exit(f"ERROR: {args.known}: no known variants")
    # a required list is a short checklist; a long list is probably a benchmark set
    if len(known) >= MAX_REQUIRED:
        print(
            f"WARNING: {args.known}: {len(known)} required variants; 'required' is meant for a "
            f"checklist of fewer than {MAX_REQUIRED} variants (use mode 'benchmark' for truth sets)",
            file=sys.stderr,
        )
    # one value of each option per technology, in the same order
    lengths = {len(args.technology), len(args.normalized), len(args.compare), len(args.regions)}
    if len(lengths) != 1:
        sys.exit(
            "ERROR: give --technology, --normalized, --compare and --regions once per "
            "technology, in the same order"
        )
    if len(set(args.technology)) != len(args.technology):
        sys.exit(f"ERROR: --technology {' '.join(args.technology)}: a technology is given twice")
    # optional tumor evidence: technology=VCF, technology must be one of --technology
    evidence_paths = {}
    for item in args.evidence:
        tech, sep, path = item.partition("=")
        if not sep or tech not in args.technology:
            sys.exit(
                f"ERROR: --evidence '{item}' must be TECHNOLOGY=VCF with a technology of "
                f"--technology ({', '.join(args.technology)})"
            )
        evidence_paths[tech] = path

    rows, summary = [], []
    # technologies in name order: the inputs arrive in any order, the tables must not depend on it
    inputs = sorted(zip(args.technology, args.normalized, args.compare, args.regions, strict=True))
    for tech, norm_path, compare_path, region_path in inputs:
        normalized = load_calls(norm_path)
        compared = load_calls(compare_path)
        regions = load_regions(region_path)
        evidence = load_evidence(evidence_paths.get(tech))
        # all normalized alleles per position, to report calls with another allele
        by_position = {}
        for key in normalized:
            by_position.setdefault(key[:2], []).append(f"{key[2]}>{key[3]}")
        passed = 0
        for chrom, pos, ref, alt, gene, label in known:
            key = (chrom, pos, ref, alt)
            region = in_regions(regions, chrom, pos, ref)
            call = normalized.get(key)
            # status: compared call > outside the region > filtered call > not called
            if key in compared:
                status = "PASS"
                passed += 1
            elif not region:
                status = "FAIL_outside_region"
            elif call:
                status = "FAIL_filtered"
            else:
                status = "FAIL_not_called"
            # read evidence in the tumor alignment (NA without alignment)
            ev_dp = ev_alt = None
            if evidence is not None:
                ev_dp = evidence[1].get((chrom, pos), 0)
                ev_alt = evidence[0].get(key, 0)
            others = [a for a in by_position.get((chrom, pos), []) if a != f"{ref}>{alt}"]
            rows.append(
                {
                    "sample": args.sample,
                    "technology": tech,
                    "label": label,
                    "gene": gene,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "status": status,
                    "in_region": "yes" if region else "no",
                    "called": "yes" if call else "no",
                    "filter": call["filter"] if call else "NA",
                    "vaf": fmt(call["vaf"]) if call else "NA",
                    "depth": fmt(call["dp"]) if call else "NA",
                    "alt_reads": fmt(call["alt"]) if call else "NA",
                    "evidence_depth": fmt(ev_dp),
                    "evidence_alt_reads": fmt(ev_alt),
                    "same_position_call": ",".join(others) or ".",
                }
            )
        print(
            f"INFO: {args.sample} {tech}: {passed} of {len(known)} required variants PASS",
            file=sys.stderr,
        )
        summary.append(
            {
                "sample": args.sample,
                "technology": tech,
                "required": len(known),
                "passed": passed,
                "failed": len(known) - passed,
                "all_passed": "yes" if passed == len(known) else "no",
            }
        )

    # both tables have at least one row (known variants and technologies are not empty)
    for name, data in (("known_variants", rows), ("known_summary", summary)):
        with open(f"{args.prefix}.{name}.tsv", "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=list(data[0]), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare", help="known variants -> sites-only VCF on the reference")
    p.add_argument(
        "--input",
        nargs="+",
        required=True,
        help="VCF(s) or TSV(s) (chrom, pos, ref, alt[, gene, label, vaf])",
    )
    p.add_argument("--vaf-tag", help="INFO tag with the VAF in VCF input (default: auto)")
    p.add_argument("--fai", required=True)
    p.add_argument("--output", required=True, help="output .vcf.gz")
    c = sub.add_parser("check", help="'required' check of one sample")
    c.add_argument("--known", required=True, help="normalized known-variant VCF")
    c.add_argument("--sample", required=True)
    c.add_argument("--technology", nargs="+", required=True)
    c.add_argument("--normalized", nargs="+", required=True, help="normalized VCF per technology")
    c.add_argument("--compare", nargs="+", required=True, help="compared VCF per technology")
    c.add_argument("--regions", nargs="+", required=True, help="comparison region per technology")
    c.add_argument(
        "--evidence", nargs="*", default=[], help="technology=tumor evidence VCF (optional)"
    )
    c.add_argument("--prefix", required=True)
    args = parser.parse_args()
    prepare(args) if args.command == "prepare" else check(args)


if __name__ == "__main__":
    main()
