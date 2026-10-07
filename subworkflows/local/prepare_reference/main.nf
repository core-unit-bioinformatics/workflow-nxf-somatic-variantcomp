//
// Prepare the reference: uncompressed FASTA with .fai, and a BED of all contigs
//

include { GUNZIP         } from '../../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_REFERENCE {

    take:
    fasta     // path: reference FASTA (plain or gzip/bgzip compressed)
    fasta_fai // path: optional .fai of a plain FASTA

    main:
    def ch_fasta = channel.of([[id: fasta.name.replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')], fasta])

    // all tools get an uncompressed FASTA (some nf-core modules only stage the FASTA itself)
    if (fasta.name.endsWith('.gz')) {
        GUNZIP(ch_fasta)
        ch_fasta = GUNZIP.out.gunzip
    }

    def ch_reference
    if (fasta_fai && !fasta.name.endsWith('.gz')) {
        ch_reference = ch_fasta.map { meta, fa -> [meta, fa, fasta_fai, []] }
    }
    else {
        SAMTOOLS_FAIDX(ch_fasta.map { meta, fa -> [meta, fa, []] }, false)
        ch_reference = ch_fasta
            .join(SAMTOOLS_FAIDX.out.fai)
            .map { meta, fa, fai -> [meta, fa, fai, []] }
    }

    // BED with the full length of every reference contig
    def ch_genome_bed = ch_reference
        .map { _meta, _fa, fai, _gzi ->
            fai.readLines().collect { line ->
                def fields = line.tokenize('\t')
                "${fields[0]}\t0\t${fields[1]}"
            }.join('\n') + '\n'
        }
        .collectFile(name: 'genome.bed')

    emit:
    reference  = ch_reference.first()  // value channel: [ meta, fasta, fai, gzi ]
    genome_bed = ch_genome_bed.first() // value channel: path
}
