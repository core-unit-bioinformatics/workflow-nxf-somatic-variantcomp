//
// Coverage of all alignments in the common region and per-sample callable regions
//
//  1. the common region is translated to the alignment's contig names for mosdepth --by
//     (rename map of the alignment's reference check in CHECK_INPUTS; check_reference.py
//     rename-bed)
//  2. mosdepth per alignment (--quantize: LOW_COVERAGE / CALLABLE bins; --by: mean depth
//     per common-region interval)
//  3. callable intervals per alignment (mosdepth_callable.py: label or mean-depth
//     threshold, contigs renamed to the reference)
//  4. per sample: intersection of the callable intervals of all its alignments
//     (BED_INTERSECT_ALL), restricted to the common region (bedtools intersect -sorted)
//  5. summary per sample (region_summary.py)
//
// A sample gets a callable region if every technology of the sample has a tumor
// alignment; other samples keep the common region.
//

include { RENAME_BED                                } from '../../../modules/local/rename_bed/main'
include { MOSDEPTH                                  } from '../../../modules/nf-core/mosdepth/main'
include { MOSDEPTH_CALLABLE                         } from '../../../modules/local/mosdepth_callable/main'
include { BED_INTERSECT_ALL                         } from '../bed_intersect_all/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_COMMON } from '../../../modules/nf-core/bedtools/intersect/main'
include { REGION_SUMMARY as REGION_SUMMARY_CALLABLE } from '../../../modules/local/region_summary/main'

workflow COVERAGE {

    take:
    ch_samplesheet // channel: samplesheet entries [ meta, vcf, vcf_index, tumor_bam, tumor_index, normal_bam, normal_index, target_bed ] (alignments indexed)
    ch_bam_rename  // channel: [ meta(alignment), rename map ] of the reference check (CHECK_INPUTS)
    ch_reference   // value channel: [ meta, fasta, fai, gzi ]
    common_bed     // value channel: common comparison region (reference contig names; sorted, merged)
    include_normal // boolean: normal alignments must be callable, too

    main:
    def ch_fai = ch_reference.map { meta, _fasta, fai, _gzi -> [meta, fai] }

    // alignments: [ meta(id, entry_id, role, sample, ...), alignment, index ]
    def ch_alignments = ch_samplesheet.flatMap { entry ->
        def meta = entry[0]
        def alignments = []
        if (entry[3]) {
            alignments << [meta + [id: "${meta.id}_tumor".toString(), entry_id: meta.id, role: 'tumor'], entry[3], entry[4]]
        }
        if (entry[5] && include_normal) {
            alignments << [meta + [id: "${meta.id}_normal".toString(), entry_id: meta.id, role: 'normal'], entry[5], entry[6]]
        }
        alignments
    }

    // rename map of each alignment, joined by the alignment id
    def ch_rename = ch_alignments
        .map { meta, _bam, _index -> [meta.id, meta] }
        .join(ch_bam_rename.map { meta, rename -> [meta.id, rename] })
        .map { _id, meta, rename -> [meta, rename] }

    // 1. common region in the alignment's contig names
    RENAME_BED(ch_rename.combine(common_bed))

    // 2. depth per alignment
    MOSDEPTH(
        ch_alignments.join(RENAME_BED.out.bed),
        ch_reference.map { meta, fasta, _fai, _gzi -> [meta, fasta] },
        ['LOW_COVERAGE', 'CALLABLE']
    )

    // 3. callable intervals per alignment (quantized bins: base mode; regions: region mode)
    MOSDEPTH_CALLABLE(
        MOSDEPTH.out.quantized_bed
            .join(MOSDEPTH.out.regions_bed)
            .join(ch_rename),
        ch_fai
    )

    // samples with a tumor alignment for every technology
    def ch_eligibility = ch_samplesheet
        .map { entry -> [entry[0].sample, entry[3] ? true : false] }
        .groupTuple()
        .map { sample, flags -> [sample, flags.every { flag -> flag }] }

    // 4. per eligible sample: bases callable in all its alignments, inside the common region
    //    (BED_INTERSECT_ALL sorts the callable BEDs; they are named by alignment)
    def ch_per_sample = MOSDEPTH_CALLABLE.out.bed
        .map { meta, bed -> [meta.sample, meta.id, bed] }
        .groupTuple()
        .join(ch_eligibility.filter { _sample, eligible -> eligible })
        .map { sample, names, beds, _eligible -> [[id: "${sample}.callable_all".toString(), sample: sample, alignments: beds.size()], beds, names] }
    BED_INTERSECT_ALL(ch_per_sample, ch_fai)
    BEDTOOLS_INTERSECT_COMMON(
        BED_INTERSECT_ALL.out.bed
            .combine(common_bed)
            .map { meta, bed, common -> [[id: meta.sample, sample: meta.sample, alignments: meta.alignments], bed, common] },
        // genome file for -sorted (sorted output in reference order)
        ch_fai
    )
    def ch_callable = BEDTOOLS_INTERSECT_COMMON.out.intersect

    // 5. summary per sample
    REGION_SUMMARY_CALLABLE(
        ch_callable.combine(common_bed).map { meta, bed, common -> [meta, common, [bed], []] },
        'callable',
        []
    )

    // comparison region per samplesheet entry: callable region of its sample, else the common region
    def ch_entries = ch_samplesheet
        .map { entry -> [entry[0].sample, entry[0]] }
        .combine(ch_eligibility, by: 0)
    def ch_regions = ch_entries
        .filter { _sample, _meta, eligible -> eligible }
        .map { sample, meta, _eligible -> [sample, meta] }
        .combine(ch_callable.map { meta, bed -> [meta.sample, bed] }, by: 0)
        .map { _sample, meta, bed -> [meta, bed] }
        .mix(
            ch_entries
                .filter { _sample, _meta, eligible -> !eligible }
                .map { _sample, meta, _eligible -> meta }
                .combine(common_bed)
        )

    emit:
    regions          = ch_regions                             // channel: [ meta, bed ] per samplesheet entry
    callable         = ch_callable                            // channel: [ meta(sample), bed ]
    callable_summary = REGION_SUMMARY_CALLABLE.out.summary    // channel: [ meta(sample), tsv ]
    multiqc_files    = MOSDEPTH.out.global_txt.mix(MOSDEPTH.out.summary_txt, MOSDEPTH.out.regions_txt)
                           .map { _meta, file -> file }       // channel: path
}
