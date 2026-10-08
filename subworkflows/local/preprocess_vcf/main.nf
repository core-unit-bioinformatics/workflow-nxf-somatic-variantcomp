//
// Reference check, normalization, harmonization and uniform filtering of every VCF
//
// Order (do not change without discussion, see CLAUDE.md):
//  1. reference check of contigs, contig lengths and sample names
//     (facts: bcftools view -h, bcftools query; decision: check_reference.py)
//  2. select tumor (+ normal) in this order, recompress to BGZF (bcftools view)
//  3. rename contigs to the reference (bcftools annotate --rename-chrs)
//  4. keep reference contigs; drop records whose only ALT is symbolic (bcftools view)
//  5. rewrite contig header lines from the .fai (bcftools reheader --fai)
//  6. split multiallelics, left-align, atomize MNVs, remove REF mismatches, deduplicate
//     (bcftools norm), then the REF-mismatch rate check (check_reference.py ref-alleles)
//  7. harmonize VAF/DP/alt counts, variant type and GT; drop remaining symbolic alleles
//     (harmonize_vcf.py)
//  8. uniform filters + restriction to the comparison region (bcftools view)
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_HEADER  } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_CONTIGS } from '../../../modules/nf-core/bcftools/query/main'
include { CHECK_REFERENCE as CHECK_REFERENCE_VCF } from '../../../modules/local/check_reference/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SELECT  } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_ANNOTATE                      } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_REHEADER                      } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CLEAN   } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_NORM                          } from '../../../modules/nf-core/bcftools/norm/main'
include { CHECK_REF_ALLELES                      } from '../../../modules/local/check_ref_alleles/main'
include { HARMONIZE_VCF                          } from '../../../modules/local/harmonize_vcf/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_COMPARE } from '../../../modules/nf-core/bcftools/view/main'

workflow PREPROCESS_VCF {

    take:
    ch_samplesheet // channel: samplesheet entries [ meta, vcf, vcf_index, tumor_bam, tumor_index, normal_bam, normal_index, target_bed ]
    ch_reference   // value channel: [ meta, fasta, fai, gzi ]
    genome_bed     // value channel: BED of all reference contigs
    ch_regions     // channel: [ meta, bed ] comparison region per entry (meta.id must match)

    main:
    // [ meta, vcf, index or [] ]
    def ch_vcf = ch_samplesheet.map { entry -> [entry[0], entry[1], entry[2]] }

    // 1. facts for the reference check: the header (contigs, samples) and the contig of every record
    BCFTOOLS_VIEW_HEADER(ch_vcf, [], [], [])
    BCFTOOLS_QUERY_CONTIGS(ch_vcf, [], [], [])
    CHECK_REFERENCE_VCF(
        BCFTOOLS_VIEW_HEADER.out.vcf.join(BCFTOOLS_QUERY_CONTIGS.out.output).map { meta, header, contigs -> [meta, [header, contigs]] },
        'vcf',
        ch_reference.map { meta, _fasta, fai, _gzi -> [meta, fai] },
        []
    )

    // 2. tumor (+ normal) in a fixed order, BGZF output (input may be plain gzip);
    //    starts only after the reference check passed
    BCFTOOLS_VIEW_SELECT(
        ch_vcf.join(CHECK_REFERENCE_VCF.out.report).map { meta, vcf, index, _report -> [meta, vcf, index] },
        [],
        [],
        []
    )

    // 3. contig renaming to the reference names
    BCFTOOLS_ANNOTATE(
        BCFTOOLS_VIEW_SELECT.out.vcf
            .join(BCFTOOLS_VIEW_SELECT.out.index)
            .join(CHECK_REFERENCE_VCF.out.rename)
            .map { meta, vcf, tbi, rename -> [meta, vcf, tbi, [], [], [], [], rename] }
    )

    // 4. reference contigs only; drop pure gVCF reference blocks
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

    // 6. split, left-align, atomize, remove REF mismatches, deduplicate; then the REF check
    BCFTOOLS_NORM(
        BCFTOOLS_REHEADER.out.vcf.join(BCFTOOLS_REHEADER.out.index),
        ch_reference.map { meta, fasta, _fai, _gzi -> [meta, fasta] }
    )
    CHECK_REF_ALLELES(BCFTOOLS_NORM.out.log)

    // 7. caller-independent tags (VC_T_VAF, VC_T_DP, ...); starts only after the REF check passed
    HARMONIZE_VCF(
        BCFTOOLS_NORM.out.vcf
            .join(BCFTOOLS_NORM.out.index)
            .join(CHECK_REF_ALLELES.out.report)
            .map { meta, vcf, tbi, _report -> [meta, vcf, tbi] }
    )

    // 8. uniform filters and restriction to the comparison region of the entry
    def ch_compare_input = HARMONIZE_VCF.out.vcf
        .map { meta, vcf, tbi -> [meta.id, meta, vcf, tbi] }
        .join(ch_regions.map { meta, bed -> [meta.id, bed] })
        .multiMap { _id, meta, vcf, tbi, bed ->
            vcf: [meta, vcf, tbi]
            bed: bed
        }
    BCFTOOLS_VIEW_COMPARE(ch_compare_input.vcf, [], ch_compare_input.bed, [])

    emit:
    normalized = HARMONIZE_VCF.out.vcf                                               // channel: [ meta, vcf, tbi ] all records
    compare    = BCFTOOLS_VIEW_COMPARE.out.vcf.join(BCFTOOLS_VIEW_COMPARE.out.index) // channel: [ meta, vcf, tbi ] filtered, in region
    reports    = CHECK_REFERENCE_VCF.out.report.mix(CHECK_REF_ALLELES.out.report)    // channel: [ meta, tsv ]
}
