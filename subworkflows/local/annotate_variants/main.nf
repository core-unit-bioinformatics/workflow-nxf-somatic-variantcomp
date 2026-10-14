//
// Annotation of all compared variants (concordant and discordant), for the
// stratification of the report
//
//  1. consequences (--vep, default): variant table -> sites-only VCF
//     (annotate_variants.py to-vcf) -> Ensembl VEP (nf-core ensemblvep/vep, offline,
//     with a VEP cache or a GTF) -> most severe consequence, gene, impact, variant
//     class, biotype per variant (bcftools +split-vep)
//  2. repeats (--annotate_repeats, default): microsatellites of the reference
//     (nf-core msisensorpro/scan, -> BED: annotate_variants.py repeats-to-bed), or a
//     precomputed --repeats_bed
//  3. variant BED (annotate_variants.py to-bed), all BEDs sorted (bedtools sort);
//     repeat unit and length per variant (bedtools map), number of overlapping
//     intervals per stratification (bedtools intersect -c)
//  4. all results merged into the variant table (annotate_variants.py merge)
//
// The design follows cio-abcd/variantinterpretation (VEP options as parameters,
// offline cache).
//

include { VARIANTS_TO_VCF                                 } from '../../../modules/local/variants_to_vcf/main'
include { ENSEMBLVEP_VEP                                  } from '../../../modules/nf-core/ensemblvep/vep/main'
include { BCFTOOLS_SPLITVEP                               } from '../../../modules/local/bcftools_splitvep/main'
include { GUNZIP as GUNZIP_FASTA                          } from '../../../modules/nf-core/gunzip/main'
include { MSISENSORPRO_SCAN                               } from '../../../modules/nf-core/msisensorpro/scan/main'
include { REPEATS_TO_BED                                  } from '../../../modules/local/repeats_to_bed/main'
include { VARIANTS_TO_BED                                 } from '../../../modules/local/variants_to_bed/main'
include { BEDTOOLS_SORT as BEDTOOLS_SORT_ANNOTATION       } from '../../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_MAP as BEDTOOLS_MAP_REPEATS            } from '../../../modules/nf-core/bedtools/map/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_STRATA } from '../../../modules/nf-core/bedtools/intersect/main'
include { ANNOTATION_MERGE                                } from '../../../modules/local/annotation_merge/main'

workflow ANNOTATE_VARIANTS {

    take:
    ch_variants      // channel: [ meta(id: 'variants'), tsv ] every compared variant of every comparison
    ch_reference     // value channel: [ meta, fasta, fai, gzi ]
    ch_strata        // channel: [ meta(id, role: 'stratum'), bed ] checked stratification BEDs
    ch_repeats_bed   // channel: [ meta(id, role: 'repeats'), bed ] checked --repeats_bed (or empty)
    use_vep          // boolean: annotate with Ensembl VEP (--vep)
    vep_cache        // path or null: VEP cache directory (--vep_cache)
    vep_gtf          // path or null: bgzipped, tabix-indexed GTF instead of a cache (--vep_gtf)
    use_repeats      // boolean: annotate repeat regions (--annotate_repeats)
    repeats_given    // boolean: --repeats_bed given (no msisensor-pro scan)

    main:
    def ch_fai = ch_reference.map { meta, _fasta, fai, _gzi -> [meta, fai] }
    def ch_multiqc = channel.empty()

    //
    // 1. consequences with Ensembl VEP
    //
    def ch_vep = channel.value([])
    if (use_vep) {
        // distinct variants as a sorted sites-only VCF
        VARIANTS_TO_VCF(ch_variants, ch_fai)
        // annotation source: the GTF (with its index) or the cache
        ENSEMBLVEP_VEP(
            VARIANTS_TO_VCF.out.vcf.map { meta, vcf, _tbi -> [meta, vcf, []] },
            params.vep_genome,
            params.vep_species,
            params.vep_cache_version,
            vep_cache ? [[id: 'vep_cache'], vep_cache] : [[], []],
            ch_reference.map { meta, fasta, _fai, _gzi -> [meta, fasta] },
            [],
            vep_gtf ? [vep_gtf, file("${vep_gtf}.tbi", checkIfExists: true)] : [[], []]
        )
        // one line per variant: most severe consequence and its gene, impact, class
        BCFTOOLS_SPLITVEP(ENSEMBLVEP_VEP.out.vcf.join(ENSEMBLVEP_VEP.out.tbi))
        ch_vep = BCFTOOLS_SPLITVEP.out.tsv.map { _meta, tsv -> tsv }
        // VEP summary (html) for MultiQC
        ch_multiqc = ch_multiqc.mix(ENSEMBLVEP_VEP.out.report.map { _meta, _process, _tool, html -> html })
    }

    //
    // 2. repeat regions of the reference: [ meta(role: 'repeats'), bed ]
    //
    def ch_repeats = channel.empty()
    if (use_repeats) {
        ch_repeats = ch_repeats_bed
        if (!repeats_given) {
            // msisensor-pro needs an uncompressed FASTA
            def ch_fasta = ch_reference.map { meta, fasta, _fai, _gzi -> [meta, fasta] }
                .branch { _meta, fasta ->
                    gzipped: fasta.name.endsWith('.gz')
                    plain: true
                }
            GUNZIP_FASTA(ch_fasta.gzipped)
            MSISENSORPRO_SCAN(ch_fasta.plain.mix(GUNZIP_FASTA.out.gunzip))
            REPEATS_TO_BED(MSISENSORPRO_SCAN.out.list)
            ch_repeats = REPEATS_TO_BED.out.bed.map { _meta, bed -> [[id: 'repeats', role: 'repeats'], bed] }
        }
    }

    //
    // 3. variants and annotation BEDs sorted in reference order (bedtools map needs it)
    //
    VARIANTS_TO_BED(ch_variants)
    BEDTOOLS_SORT_ANNOTATION(
        VARIANTS_TO_BED.out.bed
            .map { meta, bed -> [meta + [role: 'variants'], bed] }
            .mix(ch_strata, ch_repeats),
        ch_fai.map { _meta, fai -> fai }
    )
    def ch_sorted = BEDTOOLS_SORT_ANNOTATION.out.sorted.branch { meta, _bed ->
        variants: meta.role == 'variants'
        repeats: meta.role == 'repeats'
        stratum: meta.role == 'stratum'
    }

    // repeat unit and number of repeats overlapping each variant ('.' if none)
    BEDTOOLS_MAP_REPEATS(
        ch_sorted.variants.combine(ch_sorted.repeats.map { _meta, bed -> bed }),
        ch_fai
    )
    // number of stratification intervals overlapping each variant, per stratification
    BEDTOOLS_INTERSECT_STRATA(
        ch_sorted.variants
            .combine(ch_sorted.stratum)
            .map { _vmeta, variants, smeta, stratum -> [[id: smeta.id], variants, stratum] },
        [[:], []]
    )

    //
    // 4. merge the tool results into the variant table
    //
    ANNOTATION_MERGE(
        ch_variants
            .combine(ch_vep.map { tsv -> [tsv] })
            // wrapped in a list: without repeats ([]) the slot must still be filled
            // (combine adds the elements of a list, none for an empty one)
            .combine(BEDTOOLS_MAP_REPEATS.out.mapped.map { _meta, bed -> [bed] }.ifEmpty([[]]))
            // stratification results sorted by name: stable column order in the table
            .combine(BEDTOOLS_INTERSECT_STRATA.out.intersect.map { _meta, bed -> bed }.toSortedList { a, b -> a.name <=> b.name }.map { beds -> [beds] })
            .map { meta, table, vep, repeats, strata -> [meta, table, vep, repeats, strata] }
    )

    emit:
    annotated = ANNOTATION_MERGE.out.tsv // channel: [ meta, tsv ] all compared variants with annotation
    multiqc   = ch_multiqc               // channel: path (VEP summary)
}
