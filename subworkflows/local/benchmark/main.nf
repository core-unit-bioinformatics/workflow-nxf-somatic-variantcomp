//
// Compare the calls of every technology with the truth technology of the same sample
//
// Matching: rtg vcfeval (haplotype-aware, --squash-ploidy, all records, decomposed);
// bcftools isec for exact (-c none) and position-only (-c all) record matching.
// With comparison_mode 'all_pairs', all non-truth technologies are additionally compared
// pairwise (the alphabetically first technology is the baseline).
//

include { RTGTOOLS_FORMAT   } from '../../../modules/nf-core/rtgtools/format/main'
include { RTGTOOLS_VCFEVAL  } from '../../../modules/nf-core/rtgtools/vcfeval/main'
include { BCFTOOLS_ISEC as BCFTOOLS_ISEC_EXACT    } from '../../../modules/nf-core/bcftools/isec/main'
include { BCFTOOLS_ISEC as BCFTOOLS_ISEC_POSITION } from '../../../modules/nf-core/bcftools/isec/main'
include { BENCHMARK_METRICS } from '../../../modules/local/benchmark_metrics/main'
include { AGGREGATE_METRICS } from '../../../modules/local/aggregate_metrics/main'

workflow BENCHMARK {

    take:
    ch_compare      // channel: [ meta, vcf, tbi ] per samplesheet entry (filtered, in the comparison region)
    ch_regions      // channel: [ meta, bed ] comparison region per samplesheet entry
    ch_reference    // value channel: [ meta, fasta, fai, gzi ]
    comparison_mode // string: 'truth' or 'all_pairs'

    main:
    RTGTOOLS_FORMAT(ch_reference.map { meta, fasta, _fai, _gzi -> [meta, fasta, [], []] })

    def ch_entries = ch_compare
        .map { meta, vcf, tbi -> [meta.id, meta, vcf, tbi] }
        .join(ch_regions.map { meta, bed -> [meta.id, bed] })
        .map { _id, meta, vcf, tbi, bed -> [meta.sample, meta, vcf, tbi, bed] }
    def ch_truth = ch_entries.filter { _sample, meta, _vcf, _tbi, _bed -> meta.is_truth }
    def ch_query = ch_entries.filter { _sample, meta, _vcf, _tbi, _bed -> !meta.is_truth }

    def ch_pairs = ch_truth.combine(ch_query, by: 0)
        .map { _sample, base, base_vcf, base_tbi, bed, call, call_vcf, call_tbi, _call_bed ->
            [pairMeta(base, call, 'truth'), call_vcf, call_tbi, base_vcf, base_tbi, bed]
        }
    if (comparison_mode == 'all_pairs') {
        ch_pairs = ch_pairs.mix(
            ch_query.combine(ch_query, by: 0)
                .filter { _sample, base, _bv, _bt, _bb, call, _cv, _ct, _cb -> base.technology < call.technology }
                .map { _sample, base, base_vcf, base_tbi, bed, call, call_vcf, call_tbi, _call_bed ->
                    [pairMeta(base, call, 'pairwise'), call_vcf, call_tbi, base_vcf, base_tbi, bed]
                }
        )
    }

    RTGTOOLS_VCFEVAL(
        ch_pairs.map { meta, call_vcf, call_tbi, base_vcf, base_tbi, bed -> [meta, call_vcf, call_tbi, base_vcf, base_tbi, [], bed] },
        RTGTOOLS_FORMAT.out.sdf
    )

    // record matching: truth first, then query (sites.txt flags: truth, query)
    def ch_isec_input = ch_pairs.map { meta, call_vcf, call_tbi, base_vcf, base_tbi, _bed -> [meta, [base_vcf, call_vcf], [base_tbi, call_tbi], [], [], []] }
    BCFTOOLS_ISEC_EXACT(ch_isec_input)
    BCFTOOLS_ISEC_POSITION(ch_isec_input)

    BENCHMARK_METRICS(
        RTGTOOLS_VCFEVAL.out.baseline_vcf
            .join(RTGTOOLS_VCFEVAL.out.fn_vcf)
            .join(RTGTOOLS_VCFEVAL.out.tp_vcf)
            .join(RTGTOOLS_VCFEVAL.out.fp_vcf)
            .join(BCFTOOLS_ISEC_EXACT.out.results)
            .join(BCFTOOLS_ISEC_POSITION.out.results)
            .map { meta, tpb, fn, tp, fp, exact, position -> [meta, tpb, fn, tp, fp, exact, position, []] }
    )

    AGGREGATE_METRICS(
        BENCHMARK_METRICS.out.metrics
            .map { _meta, tsv -> tsv }
            .collect()
            .map { files -> [[id: 'metrics'], files] }
    )

    emit:
    metrics    = BENCHMARK_METRICS.out.metrics    // channel: [ meta(pair), tsv ]
    variants   = BENCHMARK_METRICS.out.variants   // channel: [ meta(pair), tsv ]
    per_pair   = AGGREGATE_METRICS.out.per_pair   // channel: [ meta, tsv ]
    aggregated = AGGREGATE_METRICS.out.aggregated // channel: [ meta, tsv ]
    summary    = RTGTOOLS_VCFEVAL.out.summary     // channel: [ meta(pair), txt ]
    sdf        = RTGTOOLS_FORMAT.out.sdf          // channel: [ meta, sdf ]
}

//
// Meta map of a baseline/call pair
//
def pairMeta(base, call, pair_type) {
    return [
        id         : "${base.sample}_${base.technology}_vs_${call.technology}".toString(),
        sample     : base.sample,
        truth      : base.technology,
        query      : call.technology,
        truth_tumor: base.tumor_id,
        query_tumor: call.tumor_id,
        pair_type  : pair_type,
    ]
}
