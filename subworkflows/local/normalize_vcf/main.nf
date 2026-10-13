//
// Steps 3-7 of the normalization chain for one set of records (PREPROCESS_VCF calls it
// for the PASS records and, separately, for the non-PASS records):
//
//  3. rename contigs to the reference (bcftools annotate --rename-chrs)
//  4. keep reference contigs; drop records whose only ALT is symbolic (bcftools view)
//  5. rewrite contig header lines from the .fai (bcftools reheader --fai)
//  6. split multiallelics, left-align, atomize MNVs, remove REF mismatches, deduplicate
//     (bcftools norm; its log with the counts goes to the report)
//  7. harmonize VAF/DP/alt counts, variant type and GT; drop remaining symbolic alleles
//     (harmonize_vcf.py)
//

include { BCFTOOLS_ANNOTATE                    } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CLEAN } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_REHEADER                    } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_NORM                        } from '../../../modules/nf-core/bcftools/norm/main'
include { HARMONIZE_VCF                        } from '../../../modules/local/harmonize_vcf/main'

workflow NORMALIZE_VCF {

    take:
    ch_vcf        // channel: [ meta, vcf, tbi ] tumor (+ normal) columns selected, BGZF
    ch_vcf_rename // channel: [ meta, rename map ] of the reference check (CHECK_INPUTS)
    ch_reference  // value channel: [ meta, fasta, fai, gzi ]
    genome_bed    // value channel: BED of all compared reference contigs

    main:
    // 3. contig renaming to the reference names
    BCFTOOLS_ANNOTATE(
        ch_vcf
            .join(ch_vcf_rename)
            .map { meta, vcf, tbi, rename -> [meta, vcf, tbi, [], [], [], [], rename] }
    )

    // 4. compared reference contigs only; drop pure gVCF reference blocks
    BCFTOOLS_VIEW_CLEAN(
        BCFTOOLS_ANNOTATE.out.vcf.join(BCFTOOLS_ANNOTATE.out.index),
        [],
        genome_bed,
        []
    )

    // 5. contig header lines from the reference .fai
    BCFTOOLS_REHEADER(
        BCFTOOLS_VIEW_CLEAN.out.vcf.map { meta, vcf -> [meta, vcf, [], []] },
        ch_reference.map { meta, _fasta, fai, _gzi -> [meta, fai] }
    )

    // 6. split, left-align, atomize, remove REF mismatches, deduplicate
    BCFTOOLS_NORM(
        BCFTOOLS_REHEADER.out.vcf.join(BCFTOOLS_REHEADER.out.index),
        ch_reference.map { meta, fasta, _fai, _gzi -> [meta, fasta] }
    )

    // 7. caller-independent tags (VC_T_VAF, VC_T_DP, ...)
    HARMONIZE_VCF(BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.index))

    emit:
    vcf       = HARMONIZE_VCF.out.vcf     // channel: [ meta, vcf, tbi ]
    norm_logs = BCFTOOLS_NORM.out.log     // channel: [ meta, log ] counts of bcftools norm
    summaries = HARMONIZE_VCF.out.summary // channel: [ meta, tsv ] harmonization counts
}
