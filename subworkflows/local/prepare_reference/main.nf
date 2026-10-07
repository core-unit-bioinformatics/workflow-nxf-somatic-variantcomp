//
// Prepare the reference: FASTA with .fai (and .gzi for bgzipped FASTA)
//

include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_REFERENCE {

    take:
    fasta     // path: reference FASTA (plain or bgzipped)
    fasta_fai // path: optional .fai (ignored for bgzipped FASTA, which also needs a .gzi)

    main:
    def ch_fasta = channel.of([[id: fasta.name.replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')], fasta])

    def ch_reference
    if (fasta_fai && !fasta.name.endsWith('.gz')) {
        ch_reference = ch_fasta.map { meta, fa -> [meta, fa, fasta_fai, []] }
    }
    else {
        SAMTOOLS_FAIDX(ch_fasta.map { meta, fa -> [meta, fa, []] }, false)
        ch_reference = ch_fasta
            .join(SAMTOOLS_FAIDX.out.fai)
            .join(SAMTOOLS_FAIDX.out.gzi, remainder: true)
            .map { meta, fa, fai, gzi -> [meta, fa, fai, gzi ?: []] }
    }

    emit:
    reference = ch_reference.first() // value channel: [ meta, fasta, fai, gzi ]
}
