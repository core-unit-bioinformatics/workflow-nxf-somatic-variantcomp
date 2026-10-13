//
// Reference check, normalization, harmonization and uniform filtering of every VCF
//
// Order (do not change without discussion, see CLAUDE.md):
//  1. reference check of contigs, contig lengths, sample names and REF alleles: done
//     before for all inputs (subworkflow CHECK_INPUTS); its rename maps come in here
//  2. select tumor (+ normal) in this order, recompress to BGZF (bcftools view)
//  2b. --filter_pass_only (default): split into PASS records (FILTER PASS or '.') and
//     non-PASS records (bcftools view -f / -e). Only the PASS records are compared; the
//     non-PASS records are normalized the same way, so that DISCORDANCE can tell that a
//     missed variant was called but not PASS (class filtered_non_pass)
//  3.-7. contig renaming, cleaning, contig header, bcftools norm, harmonization
//     (subworkflow NORMALIZE_VCF, once per branch)
//  ->  all normalized records: PASS and non-PASS branch concatenated (bcftools concat)
//  8. uniform filters + restriction to the comparison region (bcftools view)
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SELECT   } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_PASS     } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_NON_PASS } from '../../../modules/nf-core/bcftools/view/main'
include { NORMALIZE_VCF as NORMALIZE_PASS         } from '../normalize_vcf/main'
include { NORMALIZE_VCF as NORMALIZE_NON_PASS     } from '../normalize_vcf/main'
include { NORMALIZE_VCF as NORMALIZE_ALL          } from '../normalize_vcf/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_ALL  } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_COMPARE  } from '../../../modules/nf-core/bcftools/view/main'

workflow PREPROCESS_VCF {

    take:
    ch_samplesheet   // channel: samplesheet entries [ meta, vcf, vcf_index, tumor_bam, tumor_index, normal_bam, normal_index, target_bed ]
    ch_vcf_rename    // channel: [ meta, rename map ] of the reference check (CHECK_INPUTS)
    ch_reference     // value channel: [ meta, fasta, fai, gzi ]
    genome_bed       // value channel: BED of all compared reference contigs
    ch_regions       // channel: [ meta, bed ] comparison region per entry (meta.id must match)
    filter_pass_only // boolean: compare only PASS records (--filter_pass_only)

    main:
    // [ meta, vcf, index or [] ]
    def ch_vcf = ch_samplesheet.map { entry -> [entry[0], entry[1], entry[2]] }

    // 2. tumor (+ normal) in a fixed order, BGZF output (input may be plain gzip)
    BCFTOOLS_VIEW_SELECT(ch_vcf, [], [], [])
    def ch_selected = BCFTOOLS_VIEW_SELECT.out.vcf.join(BCFTOOLS_VIEW_SELECT.out.index)

    // [ meta, vcf, tbi ] all normalized records, and the records that are compared
    def ch_normalized = channel.empty()
    def ch_to_compare = channel.empty()
    def ch_norm_logs = channel.empty()
    def ch_summaries = channel.empty()
    if (filter_pass_only) {
        // 2b. PASS records (compared) and non-PASS records (only looked up by DISCORDANCE)
        BCFTOOLS_VIEW_PASS(ch_selected, [], [], [])
        BCFTOOLS_VIEW_NON_PASS(ch_selected, [], [], [])

        // 3.-7. the same normalization for both branches
        NORMALIZE_PASS(
            BCFTOOLS_VIEW_PASS.out.vcf.join(BCFTOOLS_VIEW_PASS.out.index),
            ch_vcf_rename,
            ch_reference,
            genome_bed
        )
        NORMALIZE_NON_PASS(
            BCFTOOLS_VIEW_NON_PASS.out.vcf.join(BCFTOOLS_VIEW_NON_PASS.out.index),
            ch_vcf_rename,
            ch_reference,
            genome_bed
        )

        // all normalized records: both branches concatenated (sorted, --allow-overlaps)
        BCFTOOLS_CONCAT_ALL(
            NORMALIZE_PASS.out.vcf
                .join(NORMALIZE_NON_PASS.out.vcf)
                .map { meta, pass_vcf, pass_tbi, non_pass_vcf, non_pass_tbi -> [meta, [pass_vcf, non_pass_vcf], [pass_tbi, non_pass_tbi]] }
        )
        ch_normalized = BCFTOOLS_CONCAT_ALL.out.vcf.join(BCFTOOLS_CONCAT_ALL.out.index)
        ch_to_compare = NORMALIZE_PASS.out.vcf
        ch_norm_logs = NORMALIZE_PASS.out.norm_logs.mix(NORMALIZE_NON_PASS.out.norm_logs)
        ch_summaries = NORMALIZE_PASS.out.summaries.mix(NORMALIZE_NON_PASS.out.summaries)
    } else {
        // 3.-7. all records together
        NORMALIZE_ALL(ch_selected, ch_vcf_rename, ch_reference, genome_bed)
        ch_normalized = NORMALIZE_ALL.out.vcf
        ch_to_compare = NORMALIZE_ALL.out.vcf
        ch_norm_logs = NORMALIZE_ALL.out.norm_logs
        ch_summaries = NORMALIZE_ALL.out.summaries
    }

    // 8. uniform filters and restriction to the comparison region of the entry
    def ch_compare_input = ch_to_compare
        .map { meta, vcf, tbi -> [meta.id, meta, vcf, tbi] }
        .join(ch_regions.map { meta, bed -> [meta.id, bed] })
        .multiMap { _id, meta, vcf, tbi, bed ->
            vcf: [meta, vcf, tbi]
            bed: bed
        }
    BCFTOOLS_VIEW_COMPARE(ch_compare_input.vcf, [], ch_compare_input.bed, [])

    emit:
    normalized = ch_normalized                                                       // channel: [ meta, vcf, tbi ] all normalized records (PASS and non-PASS)
    compare    = BCFTOOLS_VIEW_COMPARE.out.vcf.join(BCFTOOLS_VIEW_COMPARE.out.index) // channel: [ meta, vcf, tbi ] filtered, in region
    norm_logs  = ch_norm_logs                                                        // channel: [ meta, log ] counts of bcftools norm (per branch)
    summaries  = ch_summaries                                                        // channel: [ meta, tsv ] harmonization counts (per branch)
}
