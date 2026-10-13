#!/usr/bin/env python3
"""Generate the synthetic edge-case test dataset (test tier T0).

The dataset consists of a small random reference genome, DRAGEN-like somatic
tumor/normal VCFs for two samples and three technologies, target BED files and
a samplesheet. Variants are designed to exercise every normalization and
comparison step of the pipeline. The intended outcome of each variant is
written to ``design.tsv`` so that tests can check the pipeline's results.

Technologies: ``wgs`` (truth, no BED), ``wes`` (capture BED), ``panel``
(amplicon BED, a subset of the WES BED). The comparison region therefore
equals the panel BED.

Edge cases (see ``design.tsv`` for the full list):

- MNV reported as one record (wes) vs. two SNV records (wgs, panel)
- deletion in a homopolymer, left-aligned (wgs, wes) vs. right-shifted (panel)
- insertion in a dinucleotide repeat, left-aligned vs. right-shifted
- multiallelic record (wes) with one true and one extra allele
- ALT allele list containing ``<NON_REF>`` and a gVCF reference block (wes)
- duplicated record (panel)
- FILTER != PASS for a true variant (wes)
- variants outside the common target region
- known-variant lists for S1 (``known_variants.csv``): a 'required' list (one
  variant found everywhere, one filtered in wes, one outside the common region)
  and the designed true variants as a 'benchmark' truth set
- gene annotation for VEP (``annotation/genes.gtf.gz`` + .tbi: one protein-coding
  transcript per gene, exon = CDS) and stratification BEDs for the homopolymer and
  the dinucleotide repeat (``stratification.csv``)
- sample S2: contigs without ``chr`` prefix and no contig header lines
  (panel), wrong chrM length in the header (wes), normal column before tumor
  column (wgs)

Usage:
    make_synthetic_data.py --outdir <testdata-repo>/synthetic
"""

import argparse
import csv
import random
from pathlib import Path

import pysam

SEED = 20261007
CONTIGS = {"chr1": 30000, "chr2": 20000, "chrM": 16569}
NOCHR = {"chr1": "1", "chr2": "2", "chrM": "MT"}
TECHS = {
    "wgs": {"library_type": "wgs", "bed": None},
    "wes": {"library_type": "capture", "bed": "wes.bed"},
    "panel": {"library_type": "amplicon", "bed": "panel.bed"},
}
# BED intervals are 0-based, half-open
BEDS = {
    "wes.bed": [("chr1", 1000, 20000), ("chr2", 1000, 15000)],
    "panel.bed": [("chr1", 3990, 6000), ("chr1", 8000, 9000), ("chr2", 2000, 6000)],
}
# gene annotation (BED4) and stratification BEDs for the discordance annotation
GENES = [
    ("chr1", 4000, 4500, "GENE_A"),
    ("chr1", 4500, 6000, "GENE_B"),
    ("chr1", 8000, 9000, "GENE_C"),
    ("chr2", 2000, 6000, "GENE_D"),
]
HOMOPOLYMER = ("chr1", 5000, "A" * 12)  # 1-based start
DINUC_REPEAT = ("chr2", 3000, "CA" * 8)  # 1-based start


def make_reference(rng):
    """Return a dict contig -> sequence with embedded repeat contexts."""
    seqs = {c: [rng.choice("ACGT") for _ in range(n)] for c, n in CONTIGS.items()}
    for chrom, start, motif in (HOMOPOLYMER, DINUC_REPEAT):
        for i, base in enumerate(motif):
            seqs[chrom][start - 1 + i] = base
    # make sure the repeats are not extended by flanking bases
    seqs["chr1"][HOMOPOLYMER[1] - 2] = "G"
    seqs["chr1"][HOMOPOLYMER[1] - 1 + len(HOMOPOLYMER[2])] = "T"
    seqs["chr2"][DINUC_REPEAT[1] - 2] = "G"
    seqs["chr2"][DINUC_REPEAT[1] - 1 + len(DINUC_REPEAT[2])] = "T"
    return {c: "".join(s) for c, s in seqs.items()}


def other_base(ref, offset=1):
    bases = "ACGT"
    return bases[(bases.index(ref) + offset) % 4]


class Designer:
    """Collect VCF records per (sample, tech) and the design table."""

    def __init__(self, ref):
        self.ref = ref
        self.records = {}
        self.design = []

    def base(self, chrom, pos):
        return self.ref[chrom][pos - 1]

    def seq(self, chrom, pos, length):
        return self.ref[chrom][pos - 1 : pos - 1 + length]

    def add(self, sample, tech, chrom, pos, ref, alts, vaf=0.3, filt="PASS", dp=200, info=""):
        assert self.seq(chrom, pos, len(ref)) == ref, (chrom, pos, ref)
        self.records.setdefault((sample, tech), []).append(
            dict(chrom=chrom, pos=pos, ref=ref, alts=alts, vaf=vaf, filt=filt, dp=dp, info=info)
        )

    def note(self, sample, event, chrom, pos, description, expectation):
        self.design.append(
            dict(
                sample=sample,
                event=event,
                chrom=chrom,
                pos=pos,
                description=description,
                expectation=expectation,
            )
        )


def design_variants(d):
    # ---------------------------------------------------------------- S1 ---
    s = "S1"
    # E01: SNV in all technologies
    c, p = "chr1", 4100
    alt = other_base(d.base(c, p))
    for t, v in (("wgs", 0.30), ("wes", 0.32), ("panel", 0.29)):
        d.add(s, t, c, p, d.base(c, p), [alt], vaf=v)
    d.add(s, "panel", c, p, d.base(c, p), [alt], vaf=0.29)  # E09 duplicate record
    d.note(s, "E01", c, p, "SNV in all technologies", "TP wes; TP panel")
    d.note(s, "E09", c, p, "panel record duplicated", "deduplicated to one record")

    # E02: SNV missing in panel; wes ALT list carries <NON_REF>
    c, p = "chr1", 4200
    alt = other_base(d.base(c, p))
    d.add(s, "wgs", c, p, d.base(c, p), [alt])
    d.add(s, "wes", c, p, d.base(c, p), [alt, "<NON_REF>"])
    d.note(
        s, "E02", c, p, "SNV absent in panel; wes ALT has <NON_REF>", "TP wes; FN panel (absent)"
    )

    # E03: low-VAF SNV only in the deep panel
    c, p = "chr1", 4300
    d.add(s, "panel", c, p, d.base(c, p), [other_base(d.base(c, p))], vaf=0.03, dp=2000)
    d.note(s, "E03", c, p, "low-VAF SNV only in panel", "FP panel")

    # E04: MNV (two adjacent SNVs)
    c, p = "chr1", 4400
    r1, r2 = d.base(c, p), d.base(c, p + 1)
    a1, a2 = other_base(r1), other_base(r2, 2)
    for t in ("wgs", "panel"):
        d.add(s, t, c, p, r1, [a1])
        d.add(s, t, c, p + 1, r2, [a2])
    d.add(s, "wes", c, p, r1 + r2, [a1 + a2])
    d.note(
        s,
        "E04",
        c,
        p,
        "MNV: one record in wes, two SNVs in wgs/panel",
        "TP wes (2 atoms); TP panel",
    )

    # E05: deletion of one A in the homopolymer
    c, start, motif = HOMOPOLYMER
    left = start - 1
    for t in ("wgs", "wes"):
        d.add(s, t, c, left, d.seq(c, left, 2), [d.base(c, left)])
    right = start + len(motif) - 2
    d.add(s, "panel", c, right, "AA", ["A"])
    d.note(
        s,
        "E05",
        c,
        left,
        "homopolymer deletion; panel right-shifted",
        "TP wes; TP panel after norm",
    )

    # E06: multiallelic site in wes (second allele is extra)
    c, p = "chr1", 4600
    r = d.base(c, p)
    a_true, a_extra = other_base(r, 1), other_base(r, 2)
    d.add(s, "wgs", c, p, r, [a_true])
    d.add(s, "panel", c, p, r, [a_true])
    d.add(s, "wes", c, p, r, [a_true, a_extra])
    d.note(s, "E06", c, p, "multiallelic in wes", "TP wes + FP wes (extra allele); TP panel")

    # E07: SNV outside the common region (outside the panel BED)
    c, p = "chr1", 15000
    alt = other_base(d.base(c, p))
    d.add(s, "wgs", c, p, d.base(c, p), [alt])
    d.add(s, "wes", c, p, d.base(c, p), [alt])
    d.note(s, "E07", c, p, "outside common target region", "excluded")

    # E08: true variant filtered in wes
    c, p = "chr1", 4700
    alt = other_base(d.base(c, p))
    d.add(s, "wgs", c, p, d.base(c, p), [alt])
    d.add(s, "wes", c, p, d.base(c, p), [alt], filt="weak_evidence")
    d.add(s, "panel", c, p, d.base(c, p), [alt])
    d.note(s, "E08", c, p, "FILTER=weak_evidence in wes", "FN wes (present_but_filtered); TP panel")

    # E10: gVCF reference block in wes
    c, p = "chr1", 4800
    d.add(s, "wes", c, p, d.base(c, p), ["<NON_REF>"], vaf=0.0, info="END=4810")
    d.note(s, "E10", c, p, "gVCF reference block in wes", "removed")

    # E11: 3-bp deletion in all technologies
    c, p = "chr2", 2500
    for t in TECHS:
        d.add(s, t, c, p, d.seq(c, p, 4), [d.base(c, p)])
    d.note(s, "E11", c, p, "3-bp deletion in all technologies", "TP wes; TP panel")

    # E12: CA insertion in dinucleotide repeat; wes right-shifted; panel absent
    c, start, motif = DINUC_REPEAT
    left = start - 1
    d.add(s, "wgs", c, left, d.base(c, left), [d.base(c, left) + "CA"])
    right = start + len(motif) - 1
    d.add(s, "wes", c, right, d.base(c, right), [d.base(c, right) + "CA"])
    d.note(
        s,
        "E12",
        c,
        left,
        "repeat insertion; wes right-shifted; panel absent",
        "TP wes after norm; FN panel",
    )

    # E13: SNV in the second panel interval
    c, p = "chr1", 8500
    alt = other_base(d.base(c, p))
    for t in TECHS:
        d.add(s, t, c, p, d.base(c, p), [alt])
    d.note(s, "E13", c, p, "SNV in second panel interval", "TP wes; TP panel")

    # E14: SNV between panel intervals (inside WES BED)
    c, p = "chr1", 7000
    alt = other_base(d.base(c, p))
    d.add(s, "wgs", c, p, d.base(c, p), [alt])
    d.add(s, "wes", c, p, d.base(c, p), [alt])
    d.note(s, "E14", c, p, "between panel intervals", "excluded")

    # ---------------------------------------------------------------- S2 ---
    s = "S2"
    # F01: SNV in all technologies
    c, p = "chr1", 4150
    alt = other_base(d.base(c, p))
    for t in TECHS:
        d.add(s, t, c, p, d.base(c, p), [alt])
    d.note(s, "F01", c, p, "SNV in all technologies", "TP wes; TP panel")

    # F02: SNV missing in wes
    c, p = "chr1", 5500
    alt = other_base(d.base(c, p))
    d.add(s, "wgs", c, p, d.base(c, p), [alt])
    d.add(s, "panel", c, p, d.base(c, p), [alt])
    d.note(s, "F02", c, p, "SNV absent in wes", "FN wes (absent); TP panel")

    # F03: SNV in all, high VAF in panel
    c, p = "chr2", 4000
    alt = other_base(d.base(c, p))
    for t, v in (("wgs", 0.4), ("wes", 0.45), ("panel", 0.5)):
        d.add(s, t, c, p, d.base(c, p), [alt], vaf=v)
    d.note(s, "F03", c, p, "SNV in all technologies", "TP wes; TP panel")

    # F04: SNV in wes and panel, absent in truth
    c, p = "chr2", 4500
    alt = other_base(d.base(c, p))
    d.add(s, "wes", c, p, d.base(c, p), [alt], vaf=0.08)
    d.add(s, "panel", c, p, d.base(c, p), [alt], vaf=0.06)
    d.note(s, "F04", c, p, "SNV absent in truth", "FP wes; FP panel")

    # F05: mitochondrial SNV (outside all BEDs; tests MT/chrM naming)
    c, p = "chrM", 100
    alt = other_base(d.base(c, p))
    for t in TECHS:
        d.add(s, t, c, p, d.base(c, p), [alt], vaf=0.9)
    d.note(s, "F05", c, p, "chrM SNV; panel uses MT", "excluded")


# per (sample, tech) VCF peculiarities
QUIRKS = {
    ("S2", "panel"): {"nochr": True, "contig_header": False},
    ("S2", "wes"): {"contig_lengths": {"chrM": 16570}},
    ("S2", "wgs"): {"normal_first": True},
}

VCF_HEADER = """##fileformat=VCFv4.2
##source=synthetic-DRAGEN-like (workflow-nxf-somatic-variantcomp tests)
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=weak_evidence,Description="Somatic variant score does not meet threshold">
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele not already represented at this location by REF and ALT">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=SQ,Number=1,Type=Float,Description="Somatic quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
"""


def format_sample(rec, tumor):
    n_alt = len(rec["alts"])
    dp = rec["dp"]
    if tumor:
        alt_counts = [round(dp * rec["vaf"]) if a != "<NON_REF>" else 0 for a in rec["alts"]]
        if n_alt > 1 and rec["alts"][1] != "<NON_REF>":
            alt_counts = [round(dp * rec["vaf"] / n_alt)] * n_alt
        gt, sq = ("0/0", "30.00") if rec["alts"] == ["<NON_REF>"] else ("0/1", "40.00")
    else:
        alt_counts = [0] * n_alt
        gt, sq = "0/0", "0.10"
    ref_count = dp - sum(alt_counts)
    af = ",".join(f"{a / dp:.3f}" for a in alt_counts)
    ad = ",".join(str(x) for x in [ref_count, *alt_counts])
    return f"{gt}:{sq}:{ad}:{af}:{dp}"


def write_vcf(path, sample, tech, records, quirks):
    nochr = quirks.get("nochr", False)
    tumor_id, normal_id = f"{sample}_{tech}_T", f"{sample}_{tech}_N"
    lengths = {**CONTIGS, **quirks.get("contig_lengths", {})}
    order = list(CONTIGS)
    header = VCF_HEADER
    if quirks.get("contig_header", True):
        for c in order:
            header += f"##contig=<ID={NOCHR[c] if nochr else c},length={lengths[c]}>\n"
    cols = [normal_id, tumor_id] if quirks.get("normal_first") else [tumor_id, normal_id]
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(cols) + "\n"
    lines = []
    for rec in sorted(records, key=lambda r: (order.index(r["chrom"]), r["pos"])):
        info = f"DP={rec['dp'] * 2}" + (f";{rec['info']}" if rec["info"] else "")
        samples = {tumor_id: format_sample(rec, True), normal_id: format_sample(rec, False)}
        lines.append(
            "\t".join(
                [
                    NOCHR[rec["chrom"]] if nochr else rec["chrom"],
                    str(rec["pos"]),
                    ".",
                    rec["ref"],
                    ",".join(rec["alts"]),
                    ".",
                    rec["filt"],
                    info,
                    "GT:SQ:AD:AF:DP",
                    samples[cols[0]],
                    samples[cols[1]],
                ]
            )
        )
    plain = path.with_suffix("")
    plain.write_text(header + "\n".join(lines) + "\n")
    pysam.tabix_compress(str(plain), str(path), force=True)
    plain.unlink()
    pysam.tabix_index(str(path), preset="vcf", force=True)
    return tumor_id, normal_id


def write_known_variants(out, d):
    """Known-variant lists for S1: a 'required' list and a 'benchmark' truth set."""
    known = out / "known"
    known.mkdir(exist_ok=True)
    wgs = {(r["chrom"], r["pos"]): r for r in d.records[("S1", "wgs")]}

    def line(pos, chrom="chr1", gene=".", label="."):
        r = wgs[(chrom, pos)]
        return f"{chrom},{pos},{r['ref']},{r['alts'][0]},{gene},{label},{r['vaf']}\n"

    header = "chrom,pos,ref,alt,gene,label,vaf\n"
    # E01 found everywhere, E08 filtered in wes, E07 outside the common region
    (known / "required_S1.csv").write_text(
        header
        + line(4100, gene="GENE_A", label="E01")
        + line(4700, gene="GENE_B", label="E08")
        + line(15000, label="E07")
    )
    # designed true variants of S1 in the common region (as an external truth set)
    truth = [
        r for (c, p), r in sorted(wgs.items()) if (c, p) not in {("chr1", 15000), ("chr1", 7000)}
    ]
    (known / "truth_S1.csv").write_text(header + "".join(line(r["pos"], r["chrom"]) for r in truth))
    (out / "known_variants.csv").write_text(
        "name,sample,mode,variants,regions\n"
        "required_demo,S1,required,known/required_S1.csv,n/a\n"
        "designed_truth,S1,benchmark,known/truth_S1.csv,bed/panel.bed\n"
    )


def write_gtf(path, genes):
    """Bgzipped, tabix-indexed GTF: one protein-coding transcript per gene (exon = CDS).

    GTF is 1-based and inclusive; the CDS is shortened to a multiple of 3 bases, so that
    VEP can translate it (synonymous / missense consequences of the designed variants).
    """
    lines = []
    for chrom, start, end, gene in sorted(genes):
        cds_end = start + (end - start) // 3 * 3
        attributes = (
            f'gene_id "{gene}"; transcript_id "{gene}.1"; gene_name "{gene}"; '
            'gene_type "protein_coding"; gene_biotype "protein_coding"; '
            'transcript_type "protein_coding"; transcript_biotype "protein_coding";'
        )
        for feature, f_end in (("gene", end), ("transcript", end), ("exon", end), ("CDS", cds_end)):
            frame = "0" if feature == "CDS" else "."
            lines.append(
                f"{chrom}\tsynthetic\t{feature}\t{start + 1}\t{f_end}\t.\t+\t{frame}\t{attributes}\n"
            )
    path.write_text("".join(lines))
    # bgzip + tabix (gff preset: columns 1, 4, 5)
    pysam.tabix_index(str(path), preset="gff", force=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--outdir", required=True, type=Path, help="output directory")
    args = parser.parse_args()

    out = args.outdir
    for sub in ("reference", "vcf", "bed", "annotation"):
        (out / sub).mkdir(parents=True, exist_ok=True)

    ref = make_reference(random.Random(SEED))
    fasta = out / "reference" / "synthetic.fa"
    with fasta.open("w") as fh:
        for c, s in ref.items():
            fh.write(f">{c}\n")
            fh.writelines(s[i : i + 60] + "\n" for i in range(0, len(s), 60))
    pysam.faidx(str(fasta))

    for name, intervals in BEDS.items():
        with (out / "bed" / name).open("w") as fh:
            fh.writelines(f"{c}\t{s}\t{e}\n" for c, s, e in intervals)
    write_gtf(out / "annotation" / "genes.gtf", GENES)
    # stratifications: the designed repeat contexts (1-based starts -> BED; one base of flank)
    strata = {
        "homopolymer": HOMOPOLYMER,
        "tandem_repeat": DINUC_REPEAT,
    }
    rows = []
    for name, (chrom, start, motif) in strata.items():
        bed = out / "bed" / f"stratification_{name}.bed"
        bed.write_text(f"{chrom}\t{start - 2}\t{start - 1 + len(motif) + 1}\n")
        rows.append(f"{name},bed/{bed.name}\n")
    (out / "stratification.csv").write_text("name,bed\n" + "".join(rows))

    d = Designer(ref)
    design_variants(d)

    rows = []
    for (sample, tech), records in sorted(d.records.items()):
        vcf = out / "vcf" / f"{sample}_{tech}.vcf.gz"
        tumor_id, normal_id = write_vcf(vcf, sample, tech, records, QUIRKS.get((sample, tech), {}))
        bed = TECHS[tech]["bed"]
        rows.append(
            dict(
                sample=sample,
                technology=tech,
                library_type=TECHS[tech]["library_type"],
                tumor_id=tumor_id,
                normal_id=normal_id,
                caller="synthetic-dragen-like",
                vcf=f"vcf/{vcf.name}",
                tumor_bam="n/a",
                normal_bam="n/a",
                target_bed=f"bed/{bed}" if bed else "n/a",
            )
        )
    with (out / "samplesheet.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    with (out / "design.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=list(d.design[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(d.design)

    write_known_variants(out, d)


if __name__ == "__main__":
    main()
