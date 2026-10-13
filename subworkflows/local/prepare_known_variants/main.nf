//
// Known variant lists: normalization and one item per list and sample
//
// The lists were checked and converted to sites-only VCFs on the reference before
// (CHECK_INPUTS: check_sheets.py known, known_variants.py prepare with the REF check).
//
//  1. left-align and split like the calls (bcftools norm)
//  2. one item per list and sample ('all' = every sample)
//

include { BCFTOOLS_NORM as KNOWN_NORM } from '../../../modules/nf-core/bcftools/norm/main'

workflow PREPARE_KNOWN_VARIANTS {

    take:
    ch_known_lists // channel: [ meta(id: name, sample, mode), vcf, tbi, regions or [] ] checked lists (CHECK_INPUTS)
    ch_samplesheet // channel: samplesheet entries
    ch_reference   // value channel: [ meta, fasta, fai, gzi ]

    main:
    // 1. normalize like the calls (split, left-align, atomize)
    KNOWN_NORM(
        ch_known_lists.map { meta, vcf, tbi, _regions -> [meta, vcf, tbi] },
        ch_reference.map { meta, fasta, _fai, _gzi -> [meta, fasta] }
    )

    // [ meta, vcf, tbi, regions or [] ]
    def ch_lists = KNOWN_NORM.out.vcf
        .join(KNOWN_NORM.out.index)
        .join(ch_known_lists.map { meta, _vcf, _tbi, regions -> [meta, regions] })

    // one item per list and sample ('all' = every sample of the samplesheet)
    def ch_samples = ch_samplesheet.map { entry -> entry[0].sample }.unique()
    def ch_per_sample = ch_lists
        .combine(ch_samples)
        .filter { meta, _vcf, _tbi, _bed, sample -> meta.sample == 'all' || meta.sample == sample }
        .map { meta, vcf, tbi, bed, sample ->
            [meta + [id: "${meta.id}_${sample}".toString(), name: meta.id, sample: sample], vcf, tbi, bed]
        }

    emit:
    per_sample = ch_per_sample                                         // channel: [ meta(id, name, sample, mode), vcf, tbi, regions ]
    sites      = ch_per_sample
        .filter { meta, _vcf, _tbi, _bed -> meta.mode == 'required' }
        .map { meta, vcf, _tbi, _bed -> [meta, vcf] }                   // channel: [ meta(sample), vcf ] positions for the read evidence
}
