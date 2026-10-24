//
// Alignment indices: use the index next to each BAM/CRAM (found when reading the
// samplesheet); create missing ones with samtools index
//

include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'

workflow PREPARE_ALIGNMENTS {

    take:
    ch_samplesheet // channel: samplesheet entries [ meta, vcf, vcf_index, tumor_bam, tumor_index, normal_bam, normal_index, target_bed ]

    main:
    // entries whose alignments all have an index pass at once; the others wait only for
    // their own new indices (not for the indexing of other entries)
    def ch_by_index = ch_samplesheet.branch { entry ->
        complete: !((entry[3] && !entry[4]) || (entry[5] && !entry[6]))
        missing: true
    }

    // alignments without index file: [ meta(id, entry_id, role, n_missing), alignment ]
    // (n_missing: number of missing indices of the entry, to group them without waiting
    // for the end of the channel)
    def ch_missing = ch_by_index.missing.flatMap { entry ->
        def missing = []
        if (entry[3] && !entry[4]) {
            missing << [[id: "${entry[0].id}_tumor".toString(), entry_id: entry[0].id, role: 'tumor'], entry[3]]
        }
        if (entry[5] && !entry[6]) {
            missing << [[id: "${entry[0].id}_normal".toString(), entry_id: entry[0].id, role: 'normal'], entry[5]]
        }
        missing.collect { item -> [item[0] + [n_missing: missing.size()], item[1]] }
    }

    SAMTOOLS_INDEX(ch_missing)

    // new indices per entry, emitted as soon as all indices of the entry exist:
    // [ entry_id, [ role: index ] ]
    def ch_new_index = SAMTOOLS_INDEX.out.index
        .map { meta, index -> [groupKey(meta.entry_id, meta.n_missing), meta.role, index] }
        .groupTuple()
        .map { key, roles, indices -> [key.getGroupTarget(), [roles, indices].transpose().collectEntries()] }

    // entries with the created indices filled in; entries without missing index unchanged
    def ch_entries = ch_by_index.missing
        .map { entry -> [entry[0].id, entry] }
        .join(ch_new_index)
        .map { _id, entry, indices ->
            [entry[0], entry[1], entry[2], entry[3], entry[4] ?: indices.tumor, entry[5], entry[6] ?: indices.normal, entry[7]]
        }
        .mix(ch_by_index.complete)

    emit:
    samplesheet = ch_entries // channel: samplesheet entries with alignment indices
}
