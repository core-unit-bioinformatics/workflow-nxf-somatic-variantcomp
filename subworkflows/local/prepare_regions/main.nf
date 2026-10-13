//
// Harmonize all target BEDs to the reference and build the common comparison region
//
//  1. the BEDs come checked and harmonized from CHECK_INPUTS (contig names renamed,
//     unknown and mitochondrial contigs removed)
//  2. sort in reference order (bedtools sort -g <fai>)
//  3. intersection of all target BEDs (BED_INTERSECT_ALL); without any target BED
//     (only WGS) the genome BED is used
//  4. optional: intersection with --regions_bed (bedtools intersect),
//     removal of --exclude_bed (bedtools subtract)
//  5. final merge (common_targets.bed) and summary per technology (region_summary.py)
//

include { BEDTOOLS_SORT                          } from '../../../modules/nf-core/bedtools/sort/main'
include { BED_INTERSECT_ALL                      } from '../bed_intersect_all/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_REGIONS } from '../../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_SUBTRACT                      } from '../../../modules/nf-core/bedtools/subtract/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_COMMON } from '../../../modules/nf-core/bedtools/merge/main'
include { REGION_SUMMARY                         } from '../../../modules/local/region_summary/main'

workflow PREPARE_REGIONS {

    take:
    ch_samplesheet // channel: samplesheet entries [ meta, vcf, vcf_index, tumor_bam, tumor_index, normal_bam, normal_index, target_bed ]
    ch_reference   // value channel: [ meta, fasta, fai, gzi ]
    genome_bed     // value channel: BED of all compared reference contigs
    ch_beds        // channel: [ meta(id, role), bed ] checked BEDs (roles target, regions, exclude)
    use_regions    // boolean: --regions_bed given (restrict the comparison to it)
    use_exclude    // boolean: --exclude_bed given (remove it from the comparison)

    main:
    def ch_fai = ch_reference.map { meta, _fasta, fai, _gzi -> [meta, fai] }

    // technologies with their target BED (or [] for genome-wide), for the summary
    def ch_technologies = ch_samplesheet
        .map { entry -> [entry[0].technology, entry[7]] }
        .unique()

    // 2. sort in reference order
    BEDTOOLS_SORT(
        ch_beds.filter { meta, _bed -> meta.role in ['target', 'regions', 'exclude'] },
        ch_fai.map { _meta, fai -> fai }
    )
    def ch_sorted = BEDTOOLS_SORT.out.sorted.branch { meta, _bed ->
        target: meta.role == 'target'
        regions: meta.role == 'regions'
        exclude: meta.role == 'exclude'
    }

    // 3. intersection of all target BEDs; the genome if there is none
    def ch_target_beds = ch_sorted.target
        .map { _meta, bed -> bed }
        .toSortedList { a, b -> a.name <=> b.name }
    BED_INTERSECT_ALL(
        ch_target_beds
            // wrap the list: combine would otherwise flatten it into the tuple
            .map { beds -> [beds] }
            .combine(genome_bed)
            .map { beds, genome -> [[id: 'targets_intersection'], beds ?: [genome]] }
    )
    def ch_common = BED_INTERSECT_ALL.out.bed

    // 4. optional restriction and exclusion
    if (use_regions) {
        BEDTOOLS_INTERSECT_REGIONS(
            ch_common.combine(ch_sorted.regions.map { _meta, bed -> bed }).map { meta, bed, regions -> [meta, bed, regions] },
            // genome file for -sorted (sorted output in reference order)
            ch_fai
        )
        ch_common = BEDTOOLS_INTERSECT_REGIONS.out.intersect
    }
    if (use_exclude) {
        BEDTOOLS_SUBTRACT(
            ch_common.combine(ch_sorted.exclude.map { _meta, bed -> bed }).map { meta, bed, exclude -> [meta, bed, exclude] }
        )
        ch_common = BEDTOOLS_SUBTRACT.out.bed
    }
    // final merge: one name (common_targets.bed) whatever options were used
    BEDTOOLS_MERGE_COMMON(ch_common)
    def ch_common_bed = BEDTOOLS_MERGE_COMMON.out.bed.map { _meta, bed -> bed }.first()

    // 5. summary per technology
    def ch_genome_wide = ch_technologies
        .filter { _technology, bed -> !bed }
        .map { technology, _bed -> technology }
        .toSortedList()
    REGION_SUMMARY(
        ch_common_bed
            .combine(ch_target_beds.map { beds -> [beds] })
            .combine(genome_bed)
            .map { common, beds, genome -> [[id: 'common_targets'], common, beds, genome] },
        'targets',
        ch_genome_wide
    )

    emit:
    common_bed  = ch_common_bed          // value channel: path
    summary     = REGION_SUMMARY.out.summary // channel: [ meta, tsv ]
    target_beds = ch_sorted.target       // channel: [ meta, bed ]
}
