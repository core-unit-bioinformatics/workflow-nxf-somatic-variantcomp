//
// Known variant lists: conversion to VCF on the reference, normalization, regions
//

include { KNOWN_PREPARE                            } from '../../../modules/local/known_prepare/main'
include { BCFTOOLS_NORM as KNOWN_NORM              } from '../../../modules/nf-core/bcftools/norm/main'
include { CHECK_REFERENCE as CHECK_REFERENCE_KNOWN } from '../../../modules/local/check_reference/main'

workflow PREPARE_KNOWN_VARIANTS {

    take:
    ch_known       // channel: [ meta(id: name, sample, mode), [ variant files ], regions or [] ]
    ch_samplesheet // channel: samplesheet entries
    ch_reference   // value channel: [ meta, fasta, fai, gzi ]

    main:
    KNOWN_PREPARE(ch_known.map { meta, variants, _regions -> [meta, variants] }, ch_reference)
    KNOWN_NORM(
        KNOWN_PREPARE.out.vcf,
        ch_reference.map { meta, fasta, _fai, _gzi -> [meta, fasta] }
    )
    CHECK_REFERENCE_KNOWN(
        ch_known.filter { _meta, _variants, regions -> regions }.map { meta, _variants, regions -> [meta, regions] },
        'bed',
        ch_reference,
        []
    )

    def ch_lists = KNOWN_NORM.out.vcf
        .join(KNOWN_NORM.out.index)
        .join(CHECK_REFERENCE_KNOWN.out.bed, remainder: true)
        .map { meta, vcf, tbi, bed -> [meta, vcf, tbi, bed ?: []] }

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
        .map { meta, vcf, _tbi, _bed -> [meta.sample, vcf] }            // channel: [ sample, vcf ] positions for the read evidence
}
