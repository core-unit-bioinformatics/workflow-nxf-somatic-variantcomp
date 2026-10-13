//
// Intersection of N BEDs: the bases covered by every BED
//
//  1. sort every BED in reference order (bedtools sort -g <fai>; multiinter needs it)
//  2. bedtools multiinter: for every interval, how many (and which) BEDs cover it
//  3. gawk keeps the intervals covered by all BEDs ($4 == NF - 5, conf/modules)
//  4. bedtools merge joins adjacent intervals
// With a single BED, multiinter is skipped and the result is that BED, sorted and merged.
//
// Why not bedtools intersect: `bedtools intersect -a A -b B C` reports the overlap of A
// with *any* of B and C (their union), not with all of them. An intersection of N BEDs
// with bedtools intersect needs N - 1 chained calls, and a fixed chain of modules in
// Nextflow DSL2 can't adapt to the number of BEDs. bedtools multiinter does it in one
// call for any number of BEDs.
//

include { BEDTOOLS_SORT as BEDTOOLS_SORT_INPUT } from '../../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_MULTIINTER                  } from '../../../modules/nf-core/bedtools/multiinter/main'
include { GAWK as GAWK_ALL                     } from '../../../modules/nf-core/gawk/main'
include { BEDTOOLS_MERGE                       } from '../../../modules/nf-core/bedtools/merge/main'

workflow BED_INTERSECT_ALL {

    take:
    ch_beds // channel: [ meta, [ BEDs ], [ names ] ] one item per intersection; one name per BED (unique within the item)
    ch_fai  // value channel: [ meta, fai ] reference index (sort order)

    main:
    // 1. every BED on its own: [ meta(id: name, group: meta of the item, n_beds), bed ]
    def ch_files = ch_beds.flatMap { meta, beds, names ->
        def bed_list = [beds].flatten()
        def name_list = [names].flatten()
        [bed_list, name_list].transpose().collect { bed, name ->
            [[id: name, group: meta, n_beds: bed_list.size()], bed]
        }
    }
    BEDTOOLS_SORT_INPUT(ch_files, ch_fai.map { _meta, fai -> fai })

    // sorted BEDs grouped again per item (the group size is known), in name order
    def ch_sorted = BEDTOOLS_SORT_INPUT.out.sorted
        .map { meta, bed -> [groupKey(meta.group, meta.n_beds), [meta.id, bed]] }
        .groupTuple()
        .map { key, named -> [key.getGroupTarget(), named.sort { a, b -> a[0] <=> b[0] }.collect { item -> item[1] }] }

    // bedtools multiinter needs at least two files; a single BED is only merged
    def ch_input = ch_sorted.branch { _meta, beds ->
        single: beds.size() == 1
        multiple: true
    }

    // 2. per interval: chrom, start, end, number of BEDs, list of BEDs, one 0/1 column per BED
    BEDTOOLS_MULTIINTER(ch_input.multiple, [])

    // 3. keep intervals covered by all BEDs (program in conf/modules: $4 == NF - 5)
    GAWK_ALL(BEDTOOLS_MULTIINTER.out.bed.map { meta, bed -> [meta, bed, 'bed'] }, [], false)

    // 4. join adjacent intervals (of the intersection, or of the single BED)
    BEDTOOLS_MERGE(
        GAWK_ALL.out.output.mix(ch_input.single.map { meta, beds -> [meta, beds[0]] })
    )

    emit:
    bed    = BEDTOOLS_MERGE.out.bed                                                     // channel: [ meta, bed ] the intersection
    sorted = BEDTOOLS_SORT_INPUT.out.sorted.map { meta, bed -> [meta.group, meta.id, bed] } // channel: [ meta of the item, name, sorted BED ]
}
