#!/usr/bin/env python3
"""Harmonize a normalized somatic tumor/normal VCF across callers.

Input: a VCF after `bcftools norm -m-any` (one ALT per record, left-aligned), with
the tumor sample first and the optional normal sample second (order fixed by
`bcftools view -s` in the pipeline). For every record the script adds
caller-independent INFO tags, computed the same way for every technology:

  VC_T_VAF, VC_T_DP, VC_T_ALT   tumor VAF, depth and alt read count
  VC_N_VAF, VC_N_DP, VC_N_ALT   the same for the normal (only with --normal)
  VC_T_SRC, VC_N_SRC            which FORMAT fields the values came from (AD | AF_DP)
  VC_TYPE                       SNV, MNV, INS, DEL or COMPLEX
  VC_LEN                        len(ALT) - len(REF) (0 for SNV/MNV)

Field selection per sample (in this order; the first rule that applies is used):

  1. AD     FORMAT/AD has a ref and an alt value (both not missing):
              alt   = AD[1]
              VAF   = AD[1] / sum(AD)        (sum over all AD values of the record)
              depth = FORMAT/DP if present and not missing, else sum(AD)
            DP is preferred as depth because it is what the caller reports as
            coverage; callers filter reads differently for AD and DP (e.g. DRAGEN
            and Mutect2 count only informative reads in AD), so DP may be larger
            than sum(AD). VAF always uses sum(AD) so that numerator and denominator
            are counted on the same reads.
  2. AF_DP  no usable AD, but FORMAT/AF and FORMAT/DP are present:
              VAF = AF, depth = DP, alt = round(AF * DP)
            (an approximation; it is logged and counted).
  3. none   neither: the three tags of the sample are left missing.
  AD of 0,0 (no reads, e.g. filtered GATK MuTect2 3.8 calls) is source `zero_depth`:
  VAF is undefined and left missing, alt = 0, depth from DP or 0.

Checks (an error stops the script with exit status 1):
  - sample names and order are exactly --tumor [--normal];
  - every record has exactly one ALT (the input must be split with bcftools norm);
  - AD has 2 values (ref, alt) after splitting; AD values are >= 0;
  - AF is within [0, 1] (NaN, written by some callers when AF is undefined, is
    treated as missing and counted); DP is >= 0;
  - the fraction of tumor records without VAF (rule 3; zero_depth records are only
    counted, their VAF is undefined, not missing) is at most
    --max-missing-vaf (default 0: every tumor record needs a VAF). VAF strata,
    detection power and discordance classes depend on it. Callers without
    FORMAT/AD or FORMAT/AF (e.g. Strelka2 tier counts) are not supported.
A missing normal VAF is only counted (the normal can have no reads at a site).

Records whose ALT is a symbolic allele (<NON_REF>, <*>, <DEL>, ...) are dropped and
counted: they are not sequence variants.

Every record of a somatic VCF is a call of the tumor sample, so the tumor GT is
set to 0/1 and the normal GT to 0/0 (callers fill GT differently, and splitting
multiallelic records can leave alleles with GT 0/0). FILTER, ID, QUAL, INFO and all
other FORMAT fields are copied unchanged.

Outputs:
  --output         the harmonized VCF (bgzipped; --index writes a .tbi)
  --summary        key<TAB>value counts: records read/written/dropped and, per
                   sample, how often each rule was used and how often depth came
                   from DP or from sum(AD)

Why a script: `bcftools +fill-tags` computes VAF from AD but neither the fallback to
AF/DP with logging, the per-record source tag, the GT normalization for split
somatic records nor the strict checks.
"""

import argparse
import math
import sys
from collections import Counter

import pysam

INFO_LINES = [
    (
        "VC_T_VAF",
        "1",
        "Float",
        "Tumor variant allele fraction (AD[1]/sum(AD), else FORMAT/AF; variantcomp)",
    ),
    ("VC_T_DP", "1", "Integer", "Tumor depth (FORMAT/DP, else sum(AD); variantcomp)"),
    ("VC_T_ALT", "1", "Integer", "Tumor alt read count (AD[1], else round(AF*DP); variantcomp)"),
    ("VC_T_SRC", "1", "String", "Source of the tumor values: AD or AF_DP (variantcomp)"),
    (
        "VC_N_VAF",
        "1",
        "Float",
        "Normal variant allele fraction (AD[1]/sum(AD), else FORMAT/AF; variantcomp)",
    ),
    ("VC_N_DP", "1", "Integer", "Normal depth (FORMAT/DP, else sum(AD); variantcomp)"),
    ("VC_N_ALT", "1", "Integer", "Normal alt read count (AD[1], else round(AF*DP); variantcomp)"),
    ("VC_N_SRC", "1", "String", "Source of the normal values: AD or AF_DP (variantcomp)"),
    ("VC_TYPE", "1", "String", "Variant type: SNV, MNV, INS, DEL or COMPLEX (variantcomp)"),
    ("VC_LEN", "1", "Integer", "len(ALT) - len(REF) (variantcomp)"),
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


def read_counts(sample, fmt_keys, counts=None, role=""):
    """(vaf, depth, alt, source, depth_source) of one sample, following rules 1-3.

    Raises HarmonizeError for values that violate the input requirements. A FORMAT/AF
    of NaN (written e.g. by GATK MuTect2 3.8 for a normal without alt reads) is treated
    as missing and counted in counts[f"{role}_af_nan"].
    """
    # FORMAT values that the record does not have are None
    ad = sample.get("AD") if "AD" in fmt_keys else None
    dp = first(sample.get("DP")) if "DP" in fmt_keys else None
    af = first(sample.get("AF")) if "AF" in fmt_keys else None
    # NaN = the caller could not compute AF: missing, not invalid
    if af is not None and math.isnan(af):
        af = None
        if counts is not None:
            counts[f"{role}_af_nan"] += 1
    # a split record has exactly one ref and one alt depth
    if ad is not None and any(x is not None for x in ad) and len(ad) != 2:
        raise HarmonizeError(
            f"FORMAT/AD has {len(ad)} values {ad}; expected 2 (ref, alt) "
            "in a split record; was the VCF split with bcftools norm -m-any?"
        )
    if ad is not None and any(x is not None and x < 0 for x in ad):
        raise HarmonizeError(f"FORMAT/AD {ad} has a negative value")
    if dp is not None and dp < 0:
        raise HarmonizeError(f"FORMAT/DP {dp} is negative")
    if af is not None and not 0 <= af <= 1:
        raise HarmonizeError(f"FORMAT/AF {af} is outside [0, 1]")

    # rule 1: allelic depths
    if ad is not None and None not in ad:
        total = sum(ad)
        alt = ad[1]
        depth, depth_source = (dp, "DP") if dp is not None else (total, "AD_sum")
        # no reads at all (AD 0,0): the VAF is undefined; FORMAT/AF is not used instead
        if total == 0:
            return None, depth, alt, "zero_depth", depth_source
        return alt / total, depth, alt, "AD", depth_source
    # rule 2: allele fraction and depth
    if af is not None and dp is not None:
        return af, dp, round(af * dp), "AF_DP", "DP"
    # rule 3: nothing usable; keep a depth if the caller reported one
    return None, dp, None, "none", "DP" if dp is not None else "none"


def add_header_lines(header):
    """Header with the VC_* INFO lines (replacing old ones) and a GT FORMAT line."""
    header = header.copy()
    for tag, number, kind, desc in INFO_LINES:
        # a previous harmonization run must not leave stale definitions
        if tag in header.info:
            header.info.remove_header(tag)
        header.info.add(tag, number, kind, desc)
    if "GT" not in header.formats:
        header.formats.add("GT", "1", "String", "Genotype")
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


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--input", required=True, help="normalized VCF (tumor first)")
    parser.add_argument("--output", required=True, help="output VCF (.vcf.gz)")
    parser.add_argument("--tumor", required=True, help="tumor sample name")
    parser.add_argument("--normal", help="normal sample name (omit for tumor-only)")
    parser.add_argument("--index", action="store_true", help="write a tabix index")
    parser.add_argument("--summary", help="write key<TAB>value counts to this file")
    parser.add_argument(
        "--max-missing-vaf",
        type=float,
        default=0.0,
        help="maximum fraction of tumor records without VAF (default 0)",
    )
    args = parser.parse_args()
    if not 0 <= args.max_missing_vaf <= 1:
        parser.error("--max-missing-vaf must be within [0, 1]")

    vcf_in = pysam.VariantFile(args.input)
    # the sample columns must be exactly tumor (+ normal), in this order
    samples = list(vcf_in.header.samples)
    expected = [args.tumor] + ([args.normal] if args.normal else [])
    if samples != expected:
        sys.exit(
            f"ERROR: {args.input}: expected the samples {expected} in this order, found "
            f"{samples}; check tumor_id/normal_id in the samplesheet"
        )
    # AD should be declared per allele (Number=R); anything else is reported
    if "AD" in vcf_in.header.formats and vcf_in.header.formats["AD"].number != "R":
        print(
            f"WARNING: {args.input}: FORMAT/AD is declared with Number="
            f"{vcf_in.header.formats['AD'].number}, expected R",
            file=sys.stderr,
        )
    if "AD" not in vcf_in.header.formats and "AF" not in vcf_in.header.formats:
        sys.exit(
            f"ERROR: {args.input}: neither FORMAT/AD nor FORMAT/AF is defined; VAFs can't "
            "be computed. Supported callers write FORMAT/AD (or FORMAT/AF + FORMAT/DP)."
        )

    header = add_header_lines(vcf_in.header)
    roles = [("T", args.tumor)] + ([("N", args.normal)] if args.normal else [])
    counts = Counter()
    skipped_format = Counter()
    with pysam.VariantFile(args.output, "wz", header=header) as vcf_out:
        for rec in vcf_in:
            counts["records_read"] += 1
            where = f"{args.input}: {rec.chrom}:{rec.pos} {rec.ref}>{','.join(rec.alts or ('.',))}"
            # the input must be split: exactly one ALT allele per record
            if len(rec.alts or ()) != 1:
                sys.exit(
                    f"ERROR: {where}: {len(rec.alts or ())} ALT alleles; split multiallelic "
                    "records first (bcftools norm -m-any)"
                )
            # symbolic alleles (<NON_REF>, <*>, <DEL>, ...) are not sequence variants
            if rec.alts[0].startswith("<"):
                counts["records_dropped_symbolic"] += 1
                continue
            new = copy_record(vcf_out, rec, samples, skipped_format)
            values = {}
            for role, name in roles:
                try:
                    vaf, depth, alt, source, depth_source = read_counts(
                        rec.samples[name], rec.format.keys(), counts, role
                    )
                except HarmonizeError as err:
                    sys.exit(f"ERROR: {where}: sample {name}: {err}")
                # count which rule and which depth source was used for this sample
                counts[f"{role}_source_{source}"] += 1
                counts[f"{role}_depth_{depth_source}"] += 1
                values.update(
                    {f"VC_{role}_VAF": vaf, f"VC_{role}_DP": depth, f"VC_{role}_ALT": alt}
                )
                if source not in ("none", "zero_depth"):
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
            vcf_out.write(new)
            counts["records_written"] += 1
    if args.index:
        pysam.tabix_index(args.output, preset="vcf", force=True)

    # report what was done
    written = counts["records_written"]
    print(
        f"INFO: {args.input}: {counts['records_read']} records read, {written} written, "
        f"{counts['records_dropped_symbolic']} dropped (symbolic ALT)",
        file=sys.stderr,
    )
    for role, name in roles:
        print(
            f"INFO: {name}: values from AD in {counts[f'{role}_source_AD']}, from AF/DP in "
            f"{counts[f'{role}_source_AF_DP']}, no reads in {counts[f'{role}_source_zero_depth']}, "
            f"missing in {counts[f'{role}_source_none']} records; "
            f"depth from DP in {counts[f'{role}_depth_DP']}, from sum(AD) in "
            f"{counts[f'{role}_depth_AD_sum']} records",
            file=sys.stderr,
        )
    if counts["T_source_zero_depth"]:
        print(
            f"WARNING: {args.tumor}: {counts['T_source_zero_depth']} records with no tumor reads "
            "(AD 0,0); their VAF is undefined (missing)",
            file=sys.stderr,
        )
    for role, name in roles:
        if counts[f"{role}_af_nan"]:
            print(
                f"INFO: {name}: FORMAT/AF is NaN in {counts[f'{role}_af_nan']} records (treated as "
                "missing; AD is used where present)",
                file=sys.stderr,
            )
    if counts["T_source_AF_DP"]:
        print(
            f"WARNING: {args.tumor}: {counts['T_source_AF_DP']} records without FORMAT/AD; VAF taken "
            "from FORMAT/AF, alt reads approximated as round(AF * DP)",
            file=sys.stderr,
        )
    for key, n in sorted(skipped_format.items()):
        print(
            f"WARNING: FORMAT/{key} could not be copied in {n} sample values (value does not "
            "match its header definition)",
            file=sys.stderr,
        )
    if args.summary:
        with open(args.summary, "w") as fh:
            fh.write("key\tvalue\n")
            keys = ["records_read", "records_written", "records_dropped_symbolic"]
            for role, _name in roles:
                keys += [f"{role}_source_{s}" for s in ("AD", "AF_DP", "zero_depth", "none")]
                keys += [f"{role}_depth_{s}" for s in ("DP", "AD_sum", "none")]
                keys.append(f"{role}_af_nan")
            fh.writelines(f"{key}\t{counts[key]}\n" for key in keys)
    # strict: VAF-based results need a tumor VAF for (almost) every record
    missing = counts["T_source_none"]
    if written and missing / written > args.max_missing_vaf:
        sys.exit(
            f"ERROR: {args.input}: {missing} of {written} tumor records have no VAF (no usable "
            f"FORMAT/AD, no FORMAT/AF + FORMAT/DP); allowed fraction --max-missing-vaf "
            f"{args.max_missing_vaf}. Check that the caller writes FORMAT/AD."
        )


if __name__ == "__main__":
    main()
