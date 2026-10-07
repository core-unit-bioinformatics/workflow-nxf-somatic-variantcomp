//
// Harmonize all target BEDs to the reference and build the common comparison region
//

include { CHECK_REFERENCE as CHECK_REFERENCE_BED } from '../../../modules/local/check_reference/main'
include { BEDTOOLS_SORT                          } from '../../../modules/nf-core/bedtools/sort/main'
include { COMMON_REGIONS                         } from '../../../modules/local/common_regions/main'

workflow PREPARE_REGIONS {

    take:
    ch_samplesheet // channel: samplesheet entries [ meta, vcf, vcf_index, tumor_bam, tumor_bai, normal_bam, normal_bai, target_bed ]
    ch_reference   // value channel: [ meta, fasta, fai, gzi ]
    regions_bed    // path: optional, restrict the comparison to these regions
    exclude_bed    // path: optional, remove these regions from the comparison

    main:
    // one BED per technology (consistency is checked on samplesheet level)
    def ch_technologies = ch_samplesheet
        .map { entry -> [entry[0].technology, entry[7]] }
        .unique()

    def ch_beds = ch_technologies
        .filter { _technology, bed -> bed }
        .map { technology, bed -> [[id: technology, role: 'target'], bed] }
    if (regions_bed) {
        ch_beds = ch_beds.mix(channel.of([[id: 'regions_bed', role: 'regions'], regions_bed]))
    }
    if (exclude_bed) {
        ch_beds = ch_beds.mix(channel.of([[id: 'exclude_bed', role: 'exclude'], exclude_bed]))
    }

    // reference check: contig renaming, unknown contigs dropped, build mismatch -> error
    CHECK_REFERENCE_BED(ch_beds, 'bed', ch_reference, [])
    BEDTOOLS_SORT(CHECK_REFERENCE_BED.out.bed, ch_reference.map { _meta, _fasta, fai, _gzi -> fai })

    def ch_sorted = BEDTOOLS_SORT.out.sorted.branch { meta, _bed ->
        target: meta.role == 'target'
        regions: meta.role == 'regions'
        exclude: meta.role == 'exclude'
    }
    def ch_targets = ch_sorted.target
        .map { meta, bed -> [meta.id, bed] }
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { pairs -> [[id: 'common_targets'], pairs.collect { pair -> pair[1] }] }
    def ch_genome_wide = ch_technologies
        .filter { _technology, bed -> !bed }
        .map { technology, _bed -> technology }
        .toSortedList()

    COMMON_REGIONS(
        ch_targets,
        ch_genome_wide,
        ch_reference.map { _meta, _fasta, fai, _gzi -> fai },
        ch_sorted.regions.map { _meta, bed -> bed }.collect().ifEmpty([]),
        ch_sorted.exclude.map { _meta, bed -> bed }.collect().ifEmpty([]),
    )

    emit:
    common_bed = COMMON_REGIONS.out.bed.map { _meta, bed -> bed }.first() // value channel: path
    summary    = COMMON_REGIONS.out.summary                                // channel: [ meta, tsv ]
    target_beds = ch_sorted.target                                         // channel: [ meta, bed ]
    reports    = CHECK_REFERENCE_BED.out.report                            // channel: [ meta, tsv ]
}
