//
// Coverage of all alignments in the common region and per-sample callable regions
//
// A sample gets a callable region if every technology of the sample has a tumor
// alignment; the callable region is covered with >= min_coverage in all of them
// (and in the normals, if include_normal). Other samples keep the common region.
//

include { CHECK_REFERENCE as CHECK_REFERENCE_BAM } from '../../../modules/local/check_reference/main'
include { MOSDEPTH                               } from '../../../modules/nf-core/mosdepth/main'
include { CALLABLE_REGIONS                       } from '../../../modules/local/callable_regions/main'

workflow COVERAGE {

    take:
    ch_samplesheet // channel: samplesheet entries [ meta, vcf, vcf_index, tumor_bam, tumor_bai, normal_bam, normal_bai, target_bed ]
    ch_reference   // value channel: [ meta, fasta, fai, gzi ]
    common_bed     // value channel: common comparison region (reference contig names)
    include_normal // boolean: normal alignments must be callable, too
    callable_mode  // string: 'base' or 'region'
    min_coverage   // integer

    main:
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

    // all alignments are indexed (PREPARE_ALIGNMENTS)
    def ch_indexed = ch_alignments

    // reference check of the @SQ header; common region translated to the alignment's contig names
    CHECK_REFERENCE_BAM(ch_indexed.map { meta, bam, _bai -> [meta, bam] }, 'bam', ch_reference, common_bed)

    MOSDEPTH(
        ch_indexed.join(CHECK_REFERENCE_BAM.out.regions),
        ch_reference.map { meta, fasta, _fai, _gzi -> [meta, fasta] },
        ['LOW_COVERAGE', 'CALLABLE']
    )

    // samples with a tumor alignment for every technology
    def ch_eligibility = ch_samplesheet
        .map { entry -> [entry[0].sample, entry[3] ? true : false] }
        .groupTuple()
        .map { sample, flags -> [sample, flags.every { flag -> flag }] }

    def ch_coverage = callable_mode == 'region' ? MOSDEPTH.out.regions_bed : MOSDEPTH.out.quantized_bed
    def ch_per_sample = ch_coverage
        .join(CHECK_REFERENCE_BAM.out.rename)
        .map { meta, coverage, rename -> [meta.sample, coverage, rename] }
        .groupTuple()
        .join(ch_eligibility.filter { _sample, eligible -> eligible })
        .map { sample, coverage, renames, _eligible -> [[id: sample, sample: sample], coverage, renames] }

    CALLABLE_REGIONS(
        ch_per_sample,
        common_bed,
        ch_reference.map { _meta, _fasta, fai, _gzi -> fai },
        callable_mode,
        min_coverage
    )

    // comparison region per samplesheet entry: callable region of its sample, else the common region
    def ch_entries = ch_samplesheet
        .map { entry -> [entry[0].sample, entry[0]] }
        .combine(ch_eligibility, by: 0)
    def ch_regions = ch_entries
        .filter { _sample, _meta, eligible -> eligible }
        .map { sample, meta, _eligible -> [sample, meta] }
        .combine(CALLABLE_REGIONS.out.bed.map { meta, bed -> [meta.sample, bed] }, by: 0)
        .map { _sample, meta, bed -> [meta, bed] }
        .mix(
            ch_entries
                .filter { _sample, _meta, eligible -> !eligible }
                .map { _sample, meta, _eligible -> meta }
                .combine(common_bed)
        )

    emit:
    regions          = ch_regions                                                             // channel: [ meta, bed ] per samplesheet entry
    callable         = CALLABLE_REGIONS.out.bed                                               // channel: [ meta(sample), bed ]
    callable_summary = CALLABLE_REGIONS.out.summary                                           // channel: [ meta(sample), tsv ]
    reports          = CHECK_REFERENCE_BAM.out.report                                         // channel: [ meta, tsv ]
    multiqc_files    = MOSDEPTH.out.global_txt.mix(MOSDEPTH.out.summary_txt, MOSDEPTH.out.regions_txt)
                           .map { _meta, file -> file }                                       // channel: path
}
