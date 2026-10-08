//
// Intersection of N BEDs: the bases covered by every BED
//
// bedtools multiinter reports for every interval how many (and which) BEDs cover it;
// gawk keeps intervals covered by all of them; bedtools merge joins adjacent intervals.
// With a single BED, multiinter is skipped and the result is that BED, merged.
//

include { BEDTOOLS_MULTIINTER } from '../../../modules/nf-core/bedtools/multiinter/main'
include { GAWK as GAWK_ALL    } from '../../../modules/nf-core/gawk/main'
include { BEDTOOLS_MERGE      } from '../../../modules/nf-core/bedtools/merge/main'

workflow BED_INTERSECT_ALL {

    take:
    ch_beds // channel: [ meta, [ sorted BEDs ] ] one item per intersection

    main:
    // bedtools multiinter needs at least two files; a single BED is only merged
    def ch_input = ch_beds.branch { _meta, beds ->
        single: [beds].flatten().size() == 1
        multiple: true
    }

    // per interval: chrom, start, end, number of BEDs, list of BEDs, one 0/1 column per BED
    BEDTOOLS_MULTIINTER(ch_input.multiple, [])

    // keep intervals covered by all BEDs (program in conf/modules.config: $4 == NF - 5)
    GAWK_ALL(BEDTOOLS_MULTIINTER.out.bed.map { meta, bed -> [meta, bed, 'bed'] }, [], false)

    // join adjacent intervals (of the intersection, or of the single BED)
    BEDTOOLS_MERGE(
        GAWK_ALL.out.output.mix(ch_input.single.map { meta, beds -> [meta, [beds].flatten()[0]] })
    )

    emit:
    bed = BEDTOOLS_MERGE.out.bed // channel: [ meta, bed ]
}
