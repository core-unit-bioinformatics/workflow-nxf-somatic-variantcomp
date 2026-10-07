//
// Explain and annotate discordant variants
//
// 1. read evidence: pileup of every alignment at all compared variant positions of its sample
// 2. classification of every FN/FP (present_but_filtered, germline_in_normal, low_coverage,
//    subthreshold_evidence, no_evidence, absent) and power-adjusted recall
// 3. annotation of all compared variants with genes and stratifications (repeats, ...)
//

include { CHECK_REFERENCE as CHECK_REFERENCE_ANNOTATION } from '../../../modules/local/check_reference/main'
include { READ_EVIDENCE                         } from '../../../modules/local/read_evidence/main'
include { CLASSIFY_DISCORDANCE                  } from '../../../modules/local/classify_discordance/main'
include { AGGREGATE_METRICS as AGGREGATE_POWER  } from '../../../modules/local/aggregate_metrics/main'
include { ANNOTATE_VARIANTS                     } from '../../../modules/local/annotate_variants/main'

workflow DISCORDANCE {

    take:
    ch_samplesheet  // channel: samplesheet entries [ meta, vcf, vcf_index, tumor_bam, tumor_bai, normal_bam, normal_bai, target_bed ]
    ch_normalized   // channel: [ meta, vcf, tbi ] all normalized records per entry
    ch_compare      // channel: [ meta, vcf, tbi ] compared records per entry
    ch_variants     // channel: [ meta(pair), tsv ] per-variant table of each comparison
    ch_reference    // value channel: [ meta, fasta, fai, gzi ]
    gene_bed        // path: optional gene BED (BED4, name = gene)
    ch_strata       // channel: [ meta(id: name), bed ] stratification BEDs

    main:
    // annotation BEDs: reference check, contig renaming (extra columns are kept)
    def ch_annotation_beds = ch_strata.map { meta, bed -> [meta + [role: 'stratum'], bed] }
    if (gene_bed) {
        ch_annotation_beds = ch_annotation_beds.mix(channel.of([[id: 'genes', role: 'genes'], gene_bed]))
    }
    CHECK_REFERENCE_ANNOTATION(ch_annotation_beds, 'bed', ch_reference, [])
    def ch_gene_bed = CHECK_REFERENCE_ANNOTATION.out.bed
        .filter { meta, _bed -> meta.role == 'genes' }
        .map { _meta, bed -> bed }
        .collect()
        .ifEmpty([])
    def ch_strata_beds = CHECK_REFERENCE_ANNOTATION.out.bed
        .filter { meta, _bed -> meta.role == 'stratum' }
        .map { _meta, bed -> bed }
        .collect()
        .ifEmpty([])

    // compared positions per sample
    def ch_sites = ch_compare
        .map { meta, vcf, _tbi -> [meta.sample, vcf] }
        .groupTuple()

    // alignments (tumor and normal) per entry; amplicon data keeps duplicates (see modules.config)
    def ch_alignments = ch_samplesheet.flatMap { entry ->
        def meta = entry[0]
        def alignments = []
        if (entry[3]) {
            alignments << [meta + [id: "${meta.id}_tumor".toString(), entry_id: meta.id, role: 'tumor'], entry[3], entry[4]]
        }
        if (entry[5]) {
            alignments << [meta + [id: "${meta.id}_normal".toString(), entry_id: meta.id, role: 'normal'], entry[5], entry[6]]
        }
        alignments
    }

    READ_EVIDENCE(
        ch_alignments
            .filter { _meta, _bam, bai -> bai }
            .map { meta, bam, bai -> [meta.sample, meta, bam, bai] }
            .combine(ch_sites, by: 0)
            .map { _sample, meta, bam, bai, sites -> [meta, bam, bai, sites] },
        ch_reference.map { meta, fasta, fai, _gzi -> [meta, fasta, fai] }
    )

    // per entry: normalized VCF, tumor evidence (or []), normal evidence (or [])
    def ch_entry_ids = ch_samplesheet.map { entry -> [entry[0].id] }
    def ch_tumor_evidence = READ_EVIDENCE.out.vcf
        .filter { meta, _vcf, _tbi -> meta.role == 'tumor' }
        .map { meta, vcf, _tbi -> [meta.entry_id, vcf] }
    def ch_normal_evidence = READ_EVIDENCE.out.vcf
        .filter { meta, _vcf, _tbi -> meta.role == 'normal' }
        .map { meta, vcf, _tbi -> [meta.entry_id, vcf] }
    def ch_entry_data = ch_entry_ids
        .join(ch_normalized.map { meta, vcf, _tbi -> [meta.id, vcf] })
        .join(ch_tumor_evidence, remainder: true)
        .join(ch_normal_evidence, remainder: true)
        .map { id, norm, tumor, normal -> [id, norm, tumor ?: [], normal ?: []] }

    // assemble the inputs of each comparison: truth entry data + query entry data
    def ch_classify_input = ch_variants
        .map { meta, tsv -> ["${meta.sample}_${meta.truth}".toString(), meta, tsv] }
        .combine(ch_entry_data, by: 0)
        .map { _truth_id, meta, tsv, t_norm, t_tumor, t_normal ->
            ["${meta.sample}_${meta.query}".toString(), meta, tsv, t_norm, t_tumor, t_normal]
        }
        .combine(ch_entry_data, by: 0)
        .map { _query_id, meta, tsv, t_norm, t_tumor, t_normal, q_norm, q_tumor, q_normal ->
            [meta, tsv, t_norm, q_norm, t_tumor, q_tumor, [t_normal, q_normal].flatten()]
        }

    CLASSIFY_DISCORDANCE(ch_classify_input)

    AGGREGATE_POWER(
        CLASSIFY_DISCORDANCE.out.power
            .map { _meta, tsv -> tsv }
            .collect()
            .map { files -> [[id: 'power'], files] }
    )

    // all compared variants of all comparisons in one table, annotated
    def ch_all_variants = CLASSIFY_DISCORDANCE.out.discordance
        .map { _meta, tsv -> tsv }
        .collectFile(name: 'variants.tsv', keepHeader: true, skip: 1, sort: true)
        .map { tsv -> [[id: 'variants'], tsv] }

    ANNOTATE_VARIANTS(
        ch_all_variants,
        ch_gene_bed,
        ch_strata_beds,
        ch_reference.map { _meta, _fasta, fai, _gzi -> fai }
    )

    emit:
    discordance = CLASSIFY_DISCORDANCE.out.discordance // channel: [ meta(pair), tsv ]
    power       = AGGREGATE_POWER.out.aggregated       // channel: [ meta, tsv ]
    annotated   = ANNOTATE_VARIANTS.out.tsv            // channel: [ meta, tsv ] all compared variants
}
