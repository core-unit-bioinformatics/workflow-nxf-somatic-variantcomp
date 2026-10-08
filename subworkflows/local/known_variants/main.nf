//
// Known variants: 'required' check and 'benchmark' against an external truth set
//

include { KNOWN_CHECK                                } from '../../../modules/local/known_check/main'
include { RTGTOOLS_VCFEVAL as KNOWN_VCFEVAL          } from '../../../modules/nf-core/rtgtools/vcfeval/main'
include { BCFTOOLS_ISEC as KNOWN_ISEC_EXACT           } from '../../../modules/nf-core/bcftools/isec/main'
include { BCFTOOLS_ISEC as KNOWN_ISEC_POSITION        } from '../../../modules/nf-core/bcftools/isec/main'
include { BENCHMARK_METRICS as KNOWN_METRICS         } from '../../../modules/local/benchmark_metrics/main'
include { KNOWN_SUMMARY                              } from '../../../modules/local/known_summary/main'
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
    // record matching (truth list first, then the technology's calls)
    def ch_isec_input = ch_bench.map { meta, comp, ctbi, kvcf, ktbi, _kbed, _region -> [meta, [kvcf, comp], [ktbi, ctbi], [], [], []] }
    KNOWN_ISEC_EXACT(ch_isec_input)
    KNOWN_ISEC_POSITION(ch_isec_input)
    KNOWN_METRICS(
        KNOWN_VCFEVAL.out.baseline_vcf
            .join(KNOWN_VCFEVAL.out.fn_vcf)
            .join(KNOWN_VCFEVAL.out.tp_vcf)
            .join(KNOWN_VCFEVAL.out.fp_vcf)
            .join(KNOWN_ISEC_EXACT.out.results)
            .join(KNOWN_ISEC_POSITION.out.results)
            // the isec counts are restricted to the list regions and the comparison region
            .join(ch_bench.map { meta, _comp, _ctbi, _kvcf, _ktbi, kbed, region -> [meta, [kbed, region].flatten()] })
    )
    KNOWN_AGGREGATE(
        KNOWN_METRICS.out.metrics
            .map { _meta, tsv -> tsv }
            .collect()
            .map { files -> [[id: 'known_benchmark'], files] }
    )

    // all required checks in one summary; warnings for every missing variant, error with
    // --fail_on_missing_required (known_summary.py)
    KNOWN_SUMMARY(
        KNOWN_CHECK.out.variants
            .join(KNOWN_CHECK.out.summary)
            .map { meta, variants, summary -> [meta.name, variants, summary] }
            .toSortedList { a, b -> a[1].name <=> b[1].name }
            .filter { items -> items }
            .map { items ->
                [[id: 'known'], items.collect { item -> item[0] }, items.collect { item -> item[1] }, items.collect { item -> item[2] }]
            }
    )

    emit:
    required_variants  = KNOWN_SUMMARY.out.variants       // channel: [ meta, tsv ]
    required_summary   = KNOWN_SUMMARY.out.summary        // channel: [ meta, tsv ]
    multiqc            = KNOWN_SUMMARY.out.multiqc        // channel: path
    benchmark_metrics  = KNOWN_METRICS.out.metrics        // channel: [ meta, tsv ]
    benchmark_per_pair = KNOWN_AGGREGATE.out.per_pair     // channel: [ meta, tsv ]
    benchmark          = KNOWN_AGGREGATE.out.aggregated   // channel: [ meta, tsv ]
}
