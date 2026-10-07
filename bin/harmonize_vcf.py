#!/usr/bin/env python3
"""Harmonize a normalized somatic tumor/normal VCF across callers.

Expects a biallelic (split), left-aligned VCF with the tumor sample first and the
optional normal sample second. For every record it adds caller-independent INFO
tags, computed the same way for every technology:

  VC_T_VAF, VC_T_DP, VC_T_ALT   tumor VAF, depth and alt read count
  VC_N_VAF, VC_N_DP, VC_N_ALT   the same for the normal (if present)
  VC_TYPE                       SNV, MNV, INS, DEL or COMPLEX
  VC_LEN                        len(ALT) - len(REF) (0 for SNV/MNV)

Read counts come from FORMAT/AD (ref, alt). Depth is FORMAT/DP if present,
otherwise the sum of AD. VAF is alt / sum(AD); if AD is missing, FORMAT/AF and
FORMAT/DP are used. Tags that can't be computed are left missing.

Records with a symbolic ALT allele (<NON_REF>, <*>, ...) are dropped.

Every record of a somatic VCF is a call of the tumor sample, so the tumor GT is
set to 0/1 and the normal GT to 0/0 (callers differ in how they fill GT, and
splitting multiallelic records can leave alleles with GT 0/0). FILTER is kept.
"""

import argparse
import sys

import pysam

INFO_LINES = [
    ("VC_T_VAF", "1", "Float", "Tumor variant allele fraction (alt / sum(AD); variantcomp)"),
    ("VC_T_DP", "1", "Integer", "Tumor depth (FORMAT/DP or sum(AD); variantcomp)"),
    ("VC_T_ALT", "1", "Integer", "Tumor alt read count (AD[1]; variantcomp)"),
    ("VC_N_VAF", "1", "Float", "Normal variant allele fraction (alt / sum(AD); variantcomp)"),
    ("VC_N_DP", "1", "Integer", "Normal depth (FORMAT/DP or sum(AD); variantcomp)"),
    ("VC_N_ALT", "1", "Integer", "Normal alt read count (AD[1]; variantcomp)"),
    ("VC_TYPE", "1", "String", "Variant type: SNV, MNV, INS, DEL or COMPLEX (variantcomp)"),
    ("VC_LEN", "1", "Integer", "len(ALT) - len(REF) (variantcomp)"),
]


def variant_type(ref, alt):
    if len(ref) == len(alt):
        return "SNV" if len(ref) == 1 else "MNV"
    if len(ref) == 1 and alt.startswith(ref):
        return "INS"
    if len(alt) == 1 and ref.startswith(alt):
        return "DEL"
    return "COMPLEX"


def first(value):
    """Return the first element of a tuple-like FORMAT value, or the value itself."""
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def read_counts(sample, fmt_keys):
    """Return (vaf, depth, alt) for one sample; elements may be None."""
    ad = sample.get("AD") if "AD" in fmt_keys else None
    dp = first(sample.get("DP")) if "DP" in fmt_keys else None
    af = first(sample.get("AF")) if "AF" in fmt_keys else None
    if ad is not None and len(ad) >= 2 and None not in ad[:2]:
        total = sum(x for x in ad if x is not None)
        alt = ad[1]
        vaf = alt / total if total > 0 else None
        return vaf, dp if dp is not None else total, alt
    if af is not None and dp is not None:
        return af, dp, round(af * dp)
    return af, dp, None


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--input", required=True, help="normalized VCF (tumor first)")
    parser.add_argument("--output", required=True, help="output VCF (.vcf.gz)")
    parser.add_argument("--tumor", required=True, help="tumor sample name")
    parser.add_argument("--normal", help="normal sample name (omit for tumor-only)")
    parser.add_argument("--index", action="store_true", help="write a tabix index")
    args = parser.parse_args()

    vcf_in = pysam.VariantFile(args.input)
    samples = list(vcf_in.header.samples)
    expected = [args.tumor] + ([args.normal] if args.normal else [])
    if samples != expected:
        sys.exit(f"ERROR: expected samples {expected} in this order, found {samples}")

    header = vcf_in.header.copy()
    for tag, number, kind, desc in INFO_LINES:
        if tag in header.info:
            header.info.remove_header(tag)
        header.info.add(tag, number, kind, desc)
    if "GT" not in header.formats:
        header.formats.add("GT", "1", "String", "Genotype")

    n_records = n_multi = n_symbolic = 0
    with pysam.VariantFile(args.output, "wz", header=header) as vcf_out:
        for rec in vcf_in:
            n_records += 1
            if len(rec.alts or ()) != 1:
                n_multi += 1
                continue
            if rec.alts[0].startswith("<"):
                n_symbolic += 1
                continue
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
            fmt_keys = list(rec.format.keys())
            for name in samples:
                for key in fmt_keys:
                    if key == "GT":
                        continue
                    try:
                        new.samples[name][key] = rec.samples[name][key]
                    except (TypeError, ValueError):
                        pass
            t_vaf, t_dp, t_alt = read_counts(rec.samples[args.tumor], fmt_keys)
            new.samples[args.tumor]["GT"] = (0, 1)
            values = {"VC_T_VAF": t_vaf, "VC_T_DP": t_dp, "VC_T_ALT": t_alt}
            if args.normal:
                n_vaf, n_dp, n_alt = read_counts(rec.samples[args.normal], fmt_keys)
                new.samples[args.normal]["GT"] = (0, 0)
                values.update({"VC_N_VAF": n_vaf, "VC_N_DP": n_dp, "VC_N_ALT": n_alt})
            ref, alt = rec.ref, rec.alts[0]
            values["VC_TYPE"] = variant_type(ref, alt)
            values["VC_LEN"] = len(alt) - len(ref)
            for key, value in values.items():
                if value is not None:
                    new.info[key] = round(value, 6) if isinstance(value, float) else value
            vcf_out.write(new)
    if args.index:
        pysam.tabix_index(args.output, preset="vcf", force=True)
    if n_symbolic:
        print(f"dropped {n_symbolic} records with a symbolic ALT allele", file=sys.stderr)
    if n_multi:
        print(
            f"WARNING: skipped {n_multi} multiallelic records (split them first)", file=sys.stderr
        )
    print(f"records: {n_records - n_multi - n_symbolic}", file=sys.stderr)


if __name__ == "__main__":
    main()
