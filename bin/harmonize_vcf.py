#!/usr/bin/env python3
"""Harmonize a normalized somatic tumor/normal VCF across callers.

INPUT (supported)
  VCF 4.1-4.3 (VCF/BCF, read with pysam/HTSlib) after `bcftools norm -m-any`: one ALT
  per record, left-aligned, the tumor sample first and the optional normal sample
  second (order fixed by `bcftools view -s` in the pipeline). The sample columns need
  FORMAT/AF and/or FORMAT/AD (Number=R); FORMAT/DP is used if present.
  Tested: GATK MuTect2 3.8 (SEQC2 VCFs: AD, AF Number=1, DP) and synthetic DRAGEN-like
  VCFs (AD, AF Number=A, DP). Expected to work, untested: DRAGEN 4.x, GATK Mutect2 4.x,
  other callers writing FORMAT/AF or FORMAT/AD. Not supported: callers without both
  (e.g. Strelka2 tier counts AU/CU/GU/TU, TAR/TIR).

Output INFO tags (computed the same way for every technology)
  VC_T_VAF, VC_T_DP, VC_T_ALT   tumor VAF, depth and alt read count
  VC_N_VAF, VC_N_DP, VC_N_ALT   the same for the normal (only with --normal)
  VC_T_SRC, VC_N_SRC            source of the VAF: AF (FORMAT/AF) or AD (allelic depths)
  VC_TYPE                       SNV, MNV, INS, DEL or COMPLEX
  VC_LEN                        len(ALT) - len(REF) (0 for SNV/MNV)

VAF source per sample (--vaf-source; pipeline param --vaf_source)
  auto (default)  FORMAT/AF if present (not missing, not NaN), else AD[1] / sum(AD)
  af              FORMAT/AF only
  ad              AD[1] / sum(AD) only
  In every mode:
    depth     = FORMAT/DP if present, else sum(AD)
    alt reads = AD[1] if AD is present, else round(AF * DP)
    AD 0,0 (no reads at all, e.g. filtered GATK MuTect2 3.8 calls with AF=1) is source
    `zero_depth`: the VAF is undefined and left missing, also if the caller wrote an AF.
  Records without VAF keep the tags missing ('.').

Checks (an error stops the script with exit status 1)
  - sample names and order are exactly --tumor [--normal];
  - neither FORMAT/AF nor FORMAT/AD in the header;
  - a record with more than one ALT (the input must be split with bcftools norm);
  - AD with other than 2 values after splitting, negative AD or DP, AF outside [0, 1]
    (NaN, written by some callers when AF is undefined, counts as missing).
  If more than --max-missing-vaf of the tumor records have no VAF (zero_depth records are
  only counted), a WARNING is printed and the header line `##variantcomp_vaf=unavailable`
  is written: the VAF-based results (VAF strata, detection power) of this VCF are then
  reported as NA by the following steps.

Records whose ALT is a symbolic allele (<NON_REF>, <*>, <DEL>, ...) are dropped and
counted. The tumor GT is set to 0/1 and the normal GT to 0/0 (every record is a somatic
call of the tumor; callers fill GT differently). FILTER, ID, QUAL, INFO and all other
FORMAT fields are copied unchanged.

OUTPUT
  --output    the harmonized VCF (bgzipped; --index writes a .tbi)
  --summary   two columns with header `key value`: records_read, records_written,
              records_dropped_symbolic, vaf_available (yes/no), and per sample
              (prefix T_ / N_) source_AF, source_AD, source_zero_depth, source_none,
              depth_DP, depth_AD_sum, depth_none, af_nan

Why a script: `bcftools +fill-tags` computes VAF from AD, but not the AF/AD choice with
logging, the per-record source tag, the GT normalization for split somatic records or
the strict checks.
"""

import argparse
import math
import sys
from collections import Counter

import pysam

INFO_LINES = [
    ("VC_T_VAF", "1", "Float", "Tumor variant allele fraction (FORMAT/AF or AD[1]/sum(AD))"),
    ("VC_T_DP", "1", "Integer", "Tumor depth (FORMAT/DP, else sum(AD); variantcomp)"),
    ("VC_T_ALT", "1", "Integer", "Tumor alt read count (AD[1], else round(AF*DP); variantcomp)"),
    ("VC_T_SRC", "1", "String", "Source of the tumor VAF: AF or AD (variantcomp)"),
    ("VC_N_VAF", "1", "Float", "Normal variant allele fraction (FORMAT/AF or AD[1]/sum(AD))"),
    ("VC_N_DP", "1", "Integer", "Normal depth (FORMAT/DP, else sum(AD); variantcomp)"),
    ("VC_N_ALT", "1", "Integer", "Normal alt read count (AD[1], else round(AF*DP); variantcomp)"),
    ("VC_N_SRC", "1", "String", "Source of the normal VAF: AF or AD (variantcomp)"),
    ("VC_TYPE", "1", "String", "Variant type: SNV, MNV, INS, DEL or COMPLEX (variantcomp)"),
    ("VC_LEN", "1", "Integer", "len(ALT) - len(REF) (variantcomp)"),
]
# per-sample summary keys (prefix T_ / N_)
SAMPLE_KEYS = [
    "source_AF",
    "source_AD",
    "source_zero_depth",
    "source_none",
    "depth_DP",
    "depth_AD_sum",
    "depth_none",
    "af_nan",
]


class HarmonizeError(Exception):
    """A record violates an input requirement; the message names the problem."""


def variant_type(ref, alt):
    """SNV / MNV (same length), INS / DEL (anchor base shared), else COMPLEX."""
    if len(ref) == len(alt):
        return "SNV" if len(ref) == 1 else "MNV"
    # insertion: REF is the anchor base, ALT starts with it
    if len(ref) == 1 and alt.startswith(ref):
        return "INS"
    # deletion: ALT is the anchor base, REF starts with it
    if len(alt) == 1 and ref.startswith(alt):
        return "DEL"
    return "COMPLEX"


def first(value):
    """First element of a tuple-like FORMAT value (Number=A/R/.), or the value itself."""
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def sample_fields(sample, fmt_keys):
    """(AD, DP, AF, af_was_nan) of one sample; checks the values (HarmonizeError)."""
    # FORMAT values that the record does not have are None
    ad = sample.get("AD") if "AD" in fmt_keys else None
    dp = first(sample.get("DP")) if "DP" in fmt_keys else None
    af = first(sample.get("AF")) if "AF" in fmt_keys else None
    # NaN = the caller could not compute AF: missing, not invalid
    af_nan = af is not None and math.isnan(af)
    if af_nan:
        af = None
    # AD with all values missing ('.,.') counts as missing
    if ad is not None and all(x is None for x in ad):
        ad = None
    # a split record has exactly one ref and one alt depth
    if ad is not None and len(ad) != 2:
        raise HarmonizeError(
            f"FORMAT/AD has {len(ad)} values {ad}; expected 2 (ref, alt) in a split record; "
            "was the VCF split with bcftools norm -m-any?"
        )
    if ad is not None and any(x is not None and x < 0 for x in ad):
        raise HarmonizeError(f"FORMAT/AD {ad} has a negative value")
    if dp is not None and dp < 0:
        raise HarmonizeError(f"FORMAT/DP {dp} is negative")
    if af is not None and not 0 <= af <= 1:
        raise HarmonizeError(f"FORMAT/AF {af} is outside [0, 1]")
    return ad, dp, af, af_nan


def read_counts(sample, fmt_keys, vaf_source="auto"):
    """(vaf, depth, alt, source, depth_source, af_nan) of one sample (see the module doc)."""
    ad, dp, af, af_nan = sample_fields(sample, fmt_keys)
    ad_complete = ad is not None and None not in ad
    total = sum(ad) if ad_complete else None
    # depth: FORMAT/DP, else the sum of the allelic depths
    if dp is not None:
        depth, depth_source = dp, "DP"
    elif total is not None:
        depth, depth_source = total, "AD_sum"
    else:
        depth, depth_source = None, "none"
    # alt reads: AD[1], else estimated from AF and DP
    if ad_complete:
        alt = ad[1]
    elif af is not None and dp is not None:
        alt = round(af * dp)
    else:
        alt = None
    # no reads at all: the VAF is undefined, whatever the caller wrote as AF
    if total == 0:
        return None, depth, alt, "zero_depth", depth_source, af_nan
    # the VAF from the allowed sources, in the order of the mode
    use_af = vaf_source in ("auto", "af") and af is not None
    use_ad = vaf_source in ("auto", "ad") and ad_complete
    if use_af:
        return af, depth, alt, "AF", depth_source, af_nan
    if use_ad:
        return ad[1] / total, depth, alt, "AD", depth_source, af_nan
    return None, depth, alt, "none", depth_source, af_nan


def open_input(path, tumor, normal):
    """The input VCF (pysam); checks sample names/order and the FORMAT fields."""
    vcf = pysam.VariantFile(path)
    samples = list(vcf.header.samples)
    expected = [tumor] + ([normal] if normal else [])
    if samples != expected:
        sys.exit(
            f"ERROR: {path}: expected the samples {expected} in this order, found {samples}; "
            "check tumor_id/normal_id in the samplesheet"
        )
    if "AD" not in vcf.header.formats and "AF" not in vcf.header.formats:
        sys.exit(
            f"ERROR: {path}: neither FORMAT/AD nor FORMAT/AF is defined; VAFs can't be "
            "computed. Supported callers write FORMAT/AF or FORMAT/AD."
        )
    # AD should be declared per allele (Number=R); anything else is reported
    if "AD" in vcf.header.formats and vcf.header.formats["AD"].number != "R":
        print(
            f"WARNING: {path}: FORMAT/AD is declared with Number="
            f"{vcf.header.formats['AD'].number}, expected R",
            file=sys.stderr,
        )
    return vcf


def where(path, rec):
    """Record position for messages: file: chrom:pos REF>ALT."""
    return f"{path}: {rec.chrom}:{rec.pos} {rec.ref}>{','.join(rec.alts or ('.',))}"


def check_record(path, rec):
    """True if the record is written; False for symbolic ALT; error for unsplit records."""
    # the input must be split: exactly one ALT allele per record
    if len(rec.alts or ()) != 1:
        sys.exit(
            f"ERROR: {where(path, rec)}: {len(rec.alts or ())} ALT alleles; split multiallelic "
            "records first (bcftools norm -m-any)"
        )
    # symbolic alleles (<NON_REF>, <*>, <DEL>, ...) are not sequence variants
    return not rec.alts[0].startswith("<")


def count_tumor_vaf(args):
    """First pass: (records written, tumor records without VAF) for the availability flag."""
    written = missing = 0
    with open_input(args.input, args.tumor, args.normal) as vcf:
        for rec in vcf:
            if not check_record(args.input, rec):
                continue
            written += 1
            try:
                source = read_counts(rec.samples[args.tumor], rec.format.keys(), args.vaf_source)[3]
            except HarmonizeError as err:
                sys.exit(f"ERROR: {where(args.input, rec)}: sample {args.tumor}: {err}")
            missing += source == "none"
    return written, missing


def make_header(header, vaf_available):
    """Output header: VC_* INFO lines (replacing old ones), GT line, VAF availability."""
    header = header.copy()
    for tag, number, kind, desc in INFO_LINES:
        # a previous harmonization run must not leave stale definitions
        if tag in header.info:
            header.info.remove_header(tag)
        header.info.add(tag, number, kind, desc)
    if "GT" not in header.formats:
        header.formats.add("GT", "1", "String", "Genotype")
    if not vaf_available:
        header.add_line("##variantcomp_vaf=unavailable")
    return header


def copy_record(vcf_out, rec, samples, skipped_format):
    """New output record with all fields of rec except GT (set by the caller)."""
    new = vcf_out.new_record(
        contig=rec.chrom,
        start=rec.start,
        stop=rec.stop,
        alleles=rec.alleles,
        id=rec.id,
        qual=rec.qual,
        filter=list(rec.filter.keys()),
        info=dict(rec.info),
    )
    for name in samples:
        for key in rec.format.keys():
            # GT is rewritten for every sample
            if key == "GT":
                continue
            try:
                new.samples[name][key] = rec.samples[name][key]
            except (TypeError, ValueError):
                # a value that does not fit its header definition (e.g. wrong Number)
                # can't be copied; it is counted and reported, not silently lost
                skipped_format[key] += 1
    return new


def harmonize_record(args, rec, new, roles, counts):
    """Set the VC_* tags and GT of one output record; count the sources per sample."""
    values = {}
    for role, name in roles:
        try:
            vaf, depth, alt, source, depth_source, af_nan = read_counts(
                rec.samples[name], rec.format.keys(), args.vaf_source
            )
        except HarmonizeError as err:
            sys.exit(f"ERROR: {where(args.input, rec)}: sample {name}: {err}")
        # count which source was used for this sample
        counts[f"{role}_source_{source}"] += 1
        counts[f"{role}_depth_{depth_source}"] += 1
        counts[f"{role}_af_nan"] += af_nan
        values.update({f"VC_{role}_VAF": vaf, f"VC_{role}_DP": depth, f"VC_{role}_ALT": alt})
        if source in ("AF", "AD"):
            values[f"VC_{role}_SRC"] = source
    # somatic call: tumor heterozygous, normal reference
    new.samples[args.tumor]["GT"] = (0, 1)
    if args.normal:
        new.samples[args.normal]["GT"] = (0, 0)
    values["VC_TYPE"] = variant_type(rec.ref, rec.alts[0])
    values["VC_LEN"] = len(rec.alts[0]) - len(rec.ref)
    # missing values stay missing ('.') in the INFO column
    for key, value in values.items():
        if value is not None:
            new.info[key] = round(value, 6) if isinstance(value, float) else value


def write_summary(path, counts, roles, vaf_available):
    """The key<TAB>value summary (see the module doc)."""
    with open(path, "w") as fh:
        fh.write("key\tvalue\n")
        for key in ("records_read", "records_written", "records_dropped_symbolic"):
            fh.write(f"{key}\t{counts[key]}\n")
        fh.write(f"vaf_available\t{'yes' if vaf_available else 'no'}\n")
        for role, _name in roles:
            fh.writelines(f"{role}_{key}\t{counts[f'{role}_{key}']}\n" for key in SAMPLE_KEYS)


def log_summary(args, counts, roles, skipped_format):
    """INFO/WARNING lines about the sources used and problems found."""
    print(
        f"INFO: {args.input}: {counts['records_read']} records read, "
        f"{counts['records_written']} written, {counts['records_dropped_symbolic']} dropped "
        f"(symbolic ALT); VAF source mode '{args.vaf_source}'",
        file=sys.stderr,
    )
    for role, name in roles:
        print(
            f"INFO: {name}: VAF from AF in {counts[f'{role}_source_AF']}, from AD in "
            f"{counts[f'{role}_source_AD']}, no reads in {counts[f'{role}_source_zero_depth']}, "
            f"missing in {counts[f'{role}_source_none']} records; depth from DP in "
            f"{counts[f'{role}_depth_DP']}, from sum(AD) in {counts[f'{role}_depth_AD_sum']}; "
            f"AF NaN in {counts[f'{role}_af_nan']} records",
            file=sys.stderr,
        )
    if counts["T_source_zero_depth"]:
        print(
            f"WARNING: {args.tumor}: {counts['T_source_zero_depth']} records with no tumor reads "
            "(AD 0,0); their VAF is undefined (missing)",
            file=sys.stderr,
        )
    for key, n in sorted(skipped_format.items()):
        print(
            f"WARNING: FORMAT/{key} could not be copied in {n} sample values (value does not "
            "match its header definition)",
            file=sys.stderr,
        )


def parse_args():
    """Command-line options; checks their ranges."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--input", required=True, help="normalized VCF (tumor first)")
    parser.add_argument("--output", required=True, help="output VCF (.vcf.gz)")
    parser.add_argument("--tumor", required=True, help="tumor sample name")
    parser.add_argument("--normal", help="normal sample name (omit for tumor-only)")
    parser.add_argument("--index", action="store_true", help="write a tabix index")
    parser.add_argument("--summary", help="write key<TAB>value counts to this file")
    parser.add_argument(
        "--vaf-source",
        choices=["auto", "af", "ad"],
        default="auto",
        help="VAF from FORMAT/AF (af), the allelic depths (ad) or AF, else AD (auto)",
    )
    parser.add_argument(
        "--max-missing-vaf",
        type=float,
        default=0.05,
        help="fraction of tumor records without VAF above which the VAF is unavailable",
    )
    args = parser.parse_args()
    if not 0 <= args.max_missing_vaf <= 1:
        parser.error("--max-missing-vaf must be within [0, 1]")
    return args


def main():
    """Parse the command line and run the (sub)command."""
    args = parse_args()
    roles = [("T", args.tumor)] + ([("N", args.normal)] if args.normal else [])

    # first pass: is the tumor VAF available for (nearly) every record?
    n_records, n_missing = count_tumor_vaf(args)
    vaf_available = not n_records or n_missing / n_records <= args.max_missing_vaf
    if not vaf_available:
        print(
            f"WARNING: {args.input}: {n_missing} of {n_records} tumor records have no VAF "
            f"(--vaf_source {args.vaf_source}; allowed fraction --max_missing_vaf "
            f"{args.max_missing_vaf}); the VAF-based results of this VCF are NA",
            file=sys.stderr,
        )

    # second pass: write the harmonized records
    counts, skipped_format = Counter(), Counter()
    with open_input(args.input, args.tumor, args.normal) as vcf_in:
        header = make_header(vcf_in.header, vaf_available)
        samples = list(vcf_in.header.samples)
        with pysam.VariantFile(args.output, "wz", header=header) as vcf_out:
            for rec in vcf_in:
                counts["records_read"] += 1
                if not check_record(args.input, rec):
                    counts["records_dropped_symbolic"] += 1
                    continue
                new = copy_record(vcf_out, rec, samples, skipped_format)
                harmonize_record(args, rec, new, roles, counts)
                vcf_out.write(new)
                counts["records_written"] += 1
    if args.index:
        pysam.tabix_index(args.output, preset="vcf", force=True)

    log_summary(args, counts, roles, skipped_format)
    if args.summary:
        write_summary(args.summary, counts, roles, vaf_available)


if __name__ == "__main__":
    main()
