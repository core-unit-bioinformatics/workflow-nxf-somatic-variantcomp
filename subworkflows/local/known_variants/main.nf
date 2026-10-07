//
// Known variants: 'required' check and 'benchmark' against an external truth set
//

include { KNOWN_CHECK                                } from '../../../modules/local/known_check/main'
include { RTGTOOLS_VCFEVAL as KNOWN_VCFEVAL          } from '../../../modules/nf-core/rtgtools/vcfeval/main'
include { BENCHMARK_METRICS as KNOWN_METRICS         } from '../../../modules/local/benchmark_metrics/main'
include { AGGREGATE_METRICS as KNOWN_AGGREGATE       } from '../../../modules/local/aggregate_metrics/main'

workflow KNOWN_VARIANTS {

    take:
    ch_known       // channel: [ meta(id, name, sample, mode), vcf, tbi, regions or [] ] per list and sample
    ch_normalized  // channel: [ meta, vcf, tbi ] per samplesheet entry
    ch_compare     // channel: [ meta, vcf, tbi ] per samplesheet entry
    ch_regions     // channel: [ meta, bed ] comparison region per samplesheet entry
    ch_evidence    // channel: [ meta(alignment), vcf, tbi ] read evidence
    ch_sdf         // channel: [ meta, sdf ] reference for rtg

    main:
    // per entry: [ sample, meta, normalized, compare, compare_tbi, region, tumor evidence or [] ]
    def ch_entries = ch_normalized
        .map { meta, vcf, _tbi -> [meta.id, meta, vcf] }
        .join(ch_compare.map { meta, vcf, tbi -> [meta.id, vcf, tbi] })
        .join(ch_regions.map { meta, bed -> [meta.id, bed] })
        .join(
            ch_evidence.filter { meta, _vcf, _tbi -> meta.role == 'tumor' }.map { meta, vcf, _tbi -> [meta.entry_id, vcf] },
            remainder: true
        )
        .filter { _id, meta, _norm, _comp, _ctbi, _region, _evidence -> meta != null }
        .map { _id, meta, norm, comp, ctbi, region, evidence -> [meta.sample, meta, norm, comp, ctbi, region, evidence ?: []] }

    //
    // required: every variant must be called (and pass) in every technology
    //
    def ch_by_sample = ch_entries
        .map { sample, meta, norm, comp, _ctbi, region, evidence -> [sample, meta.technology, norm, comp, region, evidence] }
        .groupTuple()
    KNOWN_CHECK(
        ch_known
            .filter { meta, _vcf, _tbi, _bed -> meta.mode == 'required' }
            .map { meta, vcf, _tbi, _bed -> [meta.sample, meta, vcf] }
            .combine(ch_by_sample, by: 0)
            .map { _sample, meta, vcf, techs, norms, comps, regions, evidence ->
                [meta, vcf, techs, norms, comps, regions, evidence.flatten()]
            }
    )

    //
    // benchmark: external truth set vs. every technology (truth technology included)
    //
    def ch_bench = ch_known
        .filter { meta, _vcf, _tbi, _bed -> meta.mode == 'benchmark' }
        .map { meta, vcf, tbi, bed -> [meta.sample, meta, vcf, tbi, bed] }
        .combine(ch_entries, by: 0)
        .map { sample, kmeta, kvcf, ktbi, kbed, emeta, _norm, comp, ctbi, region, _evidence ->
            def meta = [
                id         : "${sample}_${kmeta.name}_vs_${emeta.technology}".toString(),
                sample     : sample,
                truth      : kmeta.name,
                query      : emeta.technology,
                truth_tumor: 'ALT',
                query_tumor: emeta.tumor_id,
                pair_type  : 'known',
            ]
            [meta, comp, ctbi, kvcf, ktbi, kbed, region]
        }
    KNOWN_VCFEVAL(ch_bench, ch_sdf)
    KNOWN_METRICS(
        KNOWN_VCFEVAL.out.baseline_vcf
            .join(KNOWN_VCFEVAL.out.fn_vcf)
            .join(KNOWN_VCFEVAL.out.tp_vcf)
            .join(KNOWN_VCFEVAL.out.fp_vcf)
            .join(ch_bench.map { meta, comp, _ctbi, kvcf, _ktbi, kbed, region -> [meta, kvcf, comp, [kbed, region].flatten()] })
    )
    KNOWN_AGGREGATE(
        KNOWN_METRICS.out.metrics
            .map { _meta, tsv -> tsv }
            .collect()
            .map { files -> [[id: 'known_benchmark'], files] }
    )

    emit:
    required_variants  = KNOWN_CHECK.out.variants         // channel: [ meta, tsv ]
    required_summary   = KNOWN_CHECK.out.summary          // channel: [ meta, tsv ]
    benchmark_metrics  = KNOWN_METRICS.out.metrics        // channel: [ meta, tsv ]
    benchmark_per_pair = KNOWN_AGGREGATE.out.per_pair     // channel: [ meta, tsv ]
    benchmark          = KNOWN_AGGREGATE.out.aggregated   // channel: [ meta, tsv ]
}
