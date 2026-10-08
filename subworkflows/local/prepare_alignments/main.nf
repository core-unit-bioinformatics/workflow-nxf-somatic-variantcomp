//
// Alignment indices: use the index next to each BAM/CRAM (found when reading the
// samplesheet); create missing ones with samtools index
//

include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'

workflow PREPARE_ALIGNMENTS {

    take:
    ch_samplesheet // channel: samplesheet entries [ meta, vcf, vcf_index, tumor_bam, tumor_index, normal_bam, normal_index, target_bed ]

    main:
    // alignments without index file: [ meta(id, entry_id, role), alignment ]
    def ch_missing = ch_samplesheet.flatMap { entry ->
        def missing = []
        if (entry[3] && !entry[4]) {
            missing << [[id: "${entry[0].id}_tumor".toString(), entry_id: entry[0].id, role: 'tumor'], entry[3]]
        }
        if (entry[5] && !entry[6]) {
            missing << [[id: "${entry[0].id}_normal".toString(), entry_id: entry[0].id, role: 'normal'], entry[5]]
        }
        missing
    }

    SAMTOOLS_INDEX(ch_missing)

    // new indices per entry: [ entry_id, [ role: index ] ]
    def ch_new_index = SAMTOOLS_INDEX.out.index
        .map { meta, index -> [meta.entry_id, meta.role, index] }
        .groupTuple()
        .map { entry_id, roles, indices -> [entry_id, [roles, indices].transpose().collectEntries()] }

    // entries with the created indices filled in (entries without missing index pass unchanged)
    def ch_entries = ch_samplesheet
        .map { entry -> [entry[0].id, entry] }
        .join(ch_new_index, remainder: true)
        .filter { _id, entry, _indices -> entry != null }
        .map { _id, entry, indices ->
            indices
                ? [entry[0], entry[1], entry[2], entry[3], entry[4] ?: indices.tumor, entry[5], entry[6] ?: indices.normal, entry[7]]
                : entry
        }

    emit:
    samplesheet = ch_entries // channel: samplesheet entries with alignment indices
}
