//
// Explain discordant variants
//
//  1. read evidence: the alleles of all compared variants of a sample (bcftools query),
//     reads supporting REF, ALT or another allele of each of them in every alignment of
//     the sample (read_evidence.py; indels by haplotype sequence, independent of the
//     placement of the indel in the alignment)
//  2. classification of every FN/FP and power-adjusted recall (classify_discordance.py)
//  3. one table with every compared variant of every comparison (annotated next by
//     the subworkflow ANNOTATE_VARIANTS)
//

include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_SITES        } from '../../../modules/nf-core/bcftools/query/main'
include { READ_EVIDENCE                                 } from '../../../modules/local/read_evidence/main'
include { CLASSIFY_DISCORDANCE                          } from '../../../modules/local/classify_discordance/main'
include { AGGREGATE_METRICS as AGGREGATE_POWER          } from '../../../modules/local/aggregate_metrics/main'

workflow DISCORDANCE {

    take:
    ch_samplesheet  // channel: samplesheet entries [ meta, vcf, vcf_index, tumor_bam, tumor_index, normal_bam, normal_index, target_bed ] (alignments indexed)
    ch_normalized   // channel: [ meta, vcf, tbi ] all normalized records per entry
    ch_compare      // channel: [ meta, vcf, tbi ] compared records per entry
    ch_variants     // channel: [ meta(pair), tsv ] per-variant table of each comparison
    ch_reference    // value channel: [ meta, fasta, fai, gzi ]
    ch_extra_sites  // channel: [ meta(id, sample), vcf ] further positions for the read evidence (known variants)
    read_evidence   // boolean: count the reads of the alignments at the compared positions (--read_evidence)

    main:

    //
    // 1. read evidence
    //
    // (--read_evidence; without, the classification uses the calls only)
    // [ meta(alignment: id, entry_id, role, ...), vcf, tbi ]
    def ch_evidence = channel.empty()
    if (read_evidence) {
        // alleles (CHROM, POS, REF, ALT) of every compared variant and of the extra sites
        BCFTOOLS_QUERY_SITES(
            ch_compare.mix(ch_extra_sites.map { meta, vcf -> [meta, vcf, []] }),
            [],
            [],
            []
        )
        // one sites file per sample (concatenation; read_evidence.py removes duplicates)
        def ch_sites = BCFTOOLS_QUERY_SITES.out.output
            .collectFile { meta, sites -> ["${meta.sample}.sites.tsv", sites] }
            .map { sites -> [sites.name.replaceAll(/\.sites\.tsv$/, ''), sites] }

        // every alignment of the samplesheet: [ meta(id, entry_id, role, sample, ...), alignment, index ]
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

        // allele-specific read counts of each alignment at the variants of its sample
        READ_EVIDENCE(
            ch_alignments
                .map { meta, bam, index -> [meta.sample, meta, bam, index] }
                .combine(ch_sites, by: 0)
                .map { _sample, meta, bam, index, sites -> [meta, bam, index, sites] },
            ch_reference
        )
        ch_evidence = READ_EVIDENCE.out.vcf
    }

    //
    // 2. classification per comparison
    //
    // per entry: [ entry_id, normalized VCF, tumor evidence or [], normal evidence or [] ]
    def ch_tumor_evidence = ch_evidence
        .filter { meta, _vcf, _tbi -> meta.role == 'tumor' }
        .map { meta, vcf, _tbi -> [meta.entry_id, vcf] }
    def ch_normal_evidence = ch_evidence
        .filter { meta, _vcf, _tbi -> meta.role == 'normal' }
        .map { meta, vcf, _tbi -> [meta.entry_id, vcf] }
    def ch_entry_data = ch_samplesheet
        .map { entry -> [entry[0].id] }
        .join(ch_normalized.map { meta, vcf, _tbi -> [meta.id, vcf] })
        .join(ch_tumor_evidence, remainder: true)
        .join(ch_normal_evidence, remainder: true)
        .map { id, norm, tumor, normal -> [id, norm, tumor ?: [], normal ?: []] }

    // inputs of each comparison: truth entry data + query entry data
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

    //
    // 3. one table with every compared variant of every comparison
    //
    // one table with every variant of every comparison
    def ch_all_variants = CLASSIFY_DISCORDANCE.out.discordance
        .map { _meta, tsv -> tsv }
        // sorted by file name (<sample>_<truth>_vs_<query>): the row order must not depend on the work directories
        .collectFile(name: 'variants.tsv', keepHeader: true, skip: 1, sort: { tsv -> tsv.name })
        .map { tsv -> [[id: 'variants'], tsv] }

    emit:
    discordance = CLASSIFY_DISCORDANCE.out.discordance // channel: [ meta(pair), tsv ]
    power       = AGGREGATE_POWER.out.aggregated       // channel: [ meta, tsv ]
    variants    = ch_all_variants                      // channel: [ meta(id: 'variants'), tsv ] all compared variants
    evidence    = ch_evidence                          // channel: [ meta(alignment), vcf, tbi ]
}
