//
// Explain and annotate discordant variants
//
//  1. read evidence: positions of all compared variants of a sample (bcftools query),
//     pileup of every alignment of the sample at these positions (bcftools mpileup |
//     call -m --keep-alts, nf-core module), split and left-aligned like the calls
//     (bcftools norm)
//  2. classification of every FN/FP and power-adjusted recall (classify_discordance.py)
//  3. annotation of all compared variants: variant table -> BED (annotate_variants.py),
//     sorted (bedtools sort), genes (bedtools map), stratifications (bedtools intersect -c),
//     merged back into the table (annotate_variants.py)
//

include { CHECK_REFERENCE as CHECK_REFERENCE_ANNOTATION } from '../../../modules/local/check_reference/main'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_SITES        } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_MPILEUP                              } from '../../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_EVIDENCE       } from '../../../modules/nf-core/bcftools/norm/main'
include { CLASSIFY_DISCORDANCE                          } from '../../../modules/local/classify_discordance/main'
include { AGGREGATE_METRICS as AGGREGATE_POWER          } from '../../../modules/local/aggregate_metrics/main'
include { VARIANTS_TO_BED                               } from '../../../modules/local/variants_to_bed/main'
include { BEDTOOLS_SORT as BEDTOOLS_SORT_ANNOTATION     } from '../../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_MAP as BEDTOOLS_MAP_GENES            } from '../../../modules/nf-core/bedtools/map/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_STRATA } from '../../../modules/nf-core/bedtools/intersect/main'
include { ANNOTATION_MERGE                              } from '../../../modules/local/annotation_merge/main'

workflow DISCORDANCE {

    take:
    ch_samplesheet  // channel: samplesheet entries [ meta, vcf, vcf_index, tumor_bam, tumor_index, normal_bam, normal_index, target_bed ] (alignments indexed)
    ch_normalized   // channel: [ meta, vcf, tbi ] all normalized records per entry
    ch_compare      // channel: [ meta, vcf, tbi ] compared records per entry
    ch_variants     // channel: [ meta(pair), tsv ] per-variant table of each comparison
    ch_reference    // value channel: [ meta, fasta, fai, gzi ]
    gene_bed        // path: optional gene BED (BED4, name = gene)
    ch_strata       // channel: [ meta(id: name), bed ] stratification BEDs
    ch_extra_sites  // channel: [ meta(id, sample), vcf ] further positions for the read evidence (known variants)

    main:
    def ch_fai = ch_reference.map { meta, _fasta, fai, _gzi -> [meta, fai] }

    //
    // 1. read evidence
    //
    // positions (CHROM, POS) of every compared variant and of the extra sites
    BCFTOOLS_QUERY_SITES(
        ch_compare.mix(ch_extra_sites.map { meta, vcf -> [meta, vcf, []] }),
        [],
        [],
        []
    )
    // one sites file per sample (concatenation; mpileup accepts unsorted, duplicate sites)
    def ch_sites = BCFTOOLS_QUERY_SITES.out.output
        .collectFile { meta, sites -> ["${meta.sample}.sites.tsv", sites] }
        .map { sites -> [sites.name.replaceAll(/\.sites\.tsv$/, ''), sites] }

    // every alignment of the samplesheet: [ meta(id, entry_id, role, sample, ...), alignment, index ]
    def ch_alignments = ch_samplesheet.flatMap { entry ->
        def meta = entry[0]
        def alignments = []
        if (entry[3]) {
            alignments << [meta + [id: "${meta.id}_tumor".toString(), entry_id: meta.id, role: 'tumor'], entry[3], entry[4]]
        }
        if (entry[5]) {
            alignments << [meta + [id: "${meta.id}_normal".toString(), entry_id: meta.id, role: 'normal'], entry[5], entry[6]]
        }
        alignments
    }

    // allele counts of each alignment at the positions of its sample
    BCFTOOLS_MPILEUP(
        ch_alignments
            .map { meta, bam, _index -> [meta.sample, meta, bam] }
            .combine(ch_sites, by: 0)
            .map { _sample, meta, bam, sites -> [meta, bam, sites, []] },
        ch_reference.map { meta, fasta, fai, _gzi -> [meta, fasta, fai] },
        false
    )
    // split and left-align like the calls, so alleles match by CHROM/POS/REF/ALT
    BCFTOOLS_NORM_EVIDENCE(
        BCFTOOLS_MPILEUP.out.vcf.join(BCFTOOLS_MPILEUP.out.index),
        ch_reference.map { meta, fasta, _fai, _gzi -> [meta, fasta] }
    )
    def ch_evidence = BCFTOOLS_NORM_EVIDENCE.out.vcf.join(BCFTOOLS_NORM_EVIDENCE.out.index)

    //
    // 2. classification per comparison
    //
    // per entry: [ entry_id, normalized VCF, tumor evidence or [], normal evidence or [] ]
    def ch_tumor_evidence = ch_evidence
        .filter { meta, _vcf, _tbi -> meta.role == 'tumor' }
        .map { meta, vcf, _tbi -> [meta.entry_id, vcf] }
    def ch_normal_evidence = ch_evidence
        .filter { meta, _vcf, _tbi -> meta.role == 'normal' }
        .map { meta, vcf, _tbi -> [meta.entry_id, vcf] }
    def ch_entry_data = ch_samplesheet
        .map { entry -> [entry[0].id] }
        .join(ch_normalized.map { meta, vcf, _tbi -> [meta.id, vcf] })
        .join(ch_tumor_evidence, remainder: true)
        .join(ch_normal_evidence, remainder: true)
        .map { id, norm, tumor, normal -> [id, norm, tumor ?: [], normal ?: []] }

    // inputs of each comparison: truth entry data + query entry data
    def ch_classify_input = ch_variants
        .map { meta, tsv -> ["${meta.sample}_${meta.truth}".toString(), meta, tsv] }
        .combine(ch_entry_data, by: 0)
        .map { _truth_id, meta, tsv, t_norm, t_tumor, t_normal ->
            ["${meta.sample}_${meta.query}".toString(), meta, tsv, t_norm, t_tumor, t_normal]
        }
        .combine(ch_entry_data, by: 0)
        .map { _query_id, meta, tsv, t_norm, t_tumor, t_normal, q_norm, q_tumor, q_normal ->
            [meta, tsv, t_norm, q_norm, t_tumor, q_tumor, [t_normal, q_normal].flatten()]
        }
    CLASSIFY_DISCORDANCE(ch_classify_input)

    AGGREGATE_POWER(
        CLASSIFY_DISCORDANCE.out.power
            .map { _meta, tsv -> tsv }
            .collect()
            .map { files -> [[id: 'power'], files] }
    )

    //
    // 3. annotation of all compared variants of all comparisons
    //
    // one table with every variant of every comparison
    def ch_all_variants = CLASSIFY_DISCORDANCE.out.discordance
        .map { _meta, tsv -> tsv }
        .collectFile(name: 'variants.tsv', keepHeader: true, skip: 1, sort: true)
        .map { tsv -> [[id: 'variants'], tsv] }

    // gene and stratification BEDs: reference check (contig renaming; extra columns kept)
    def ch_annotation_beds = ch_strata.map { meta, bed -> [meta + [role: 'stratum'], bed] }
    if (gene_bed) {
        ch_annotation_beds = ch_annotation_beds.mix(channel.of([[id: 'genes', role: 'genes'], gene_bed]))
    }
    CHECK_REFERENCE_ANNOTATION(ch_annotation_beds, 'bed', ch_fai, [])

    // variants and annotation BEDs sorted in reference order (bedtools map needs both sorted)
    VARIANTS_TO_BED(ch_all_variants)
    BEDTOOLS_SORT_ANNOTATION(
        VARIANTS_TO_BED.out.bed.map { meta, bed -> [meta + [role: 'variants'], bed] }.mix(CHECK_REFERENCE_ANNOTATION.out.bed),
        ch_fai.map { _meta, fai -> fai }
    )
    def ch_sorted = BEDTOOLS_SORT_ANNOTATION.out.sorted.branch { meta, _bed ->
        variants: meta.role == 'variants'
        genes: meta.role == 'genes'
        stratum: meta.role == 'stratum'
    }

    // genes overlapping each variant (comma-separated, '.' if none)
    BEDTOOLS_MAP_GENES(
        ch_sorted.variants.combine(ch_sorted.genes.map { _meta, bed -> bed }).map { meta, variants, genes -> [meta, variants, genes] },
        ch_fai
    )
    // number of stratification intervals overlapping each variant, per stratification
    BEDTOOLS_INTERSECT_STRATA(
        ch_sorted.variants
            .combine(ch_sorted.stratum)
            .map { _vmeta, variants, smeta, stratum -> [[id: smeta.id], variants, stratum] },
        [[:], []]
    )

    // merge the bedtools results into the variant table
    ANNOTATION_MERGE(
        ch_all_variants
            .combine(BEDTOOLS_MAP_GENES.out.mapped.map { _meta, bed -> bed }.ifEmpty([]))
            // stratification results sorted by name: stable column order in the table
            .combine(BEDTOOLS_INTERSECT_STRATA.out.intersect.map { _meta, bed -> bed }.toSortedList { a, b -> a.name <=> b.name }.map { beds -> [beds] })
            .map { meta, table, genes, strata -> [meta, table, genes, strata] }
    )

    emit:
    discordance = CLASSIFY_DISCORDANCE.out.discordance // channel: [ meta(pair), tsv ]
    power       = AGGREGATE_POWER.out.aggregated       // channel: [ meta, tsv ]
    annotated   = ANNOTATION_MERGE.out.tsv             // channel: [ meta, tsv ] all compared variants
    evidence    = ch_evidence                          // channel: [ meta(alignment), vcf, tbi ]
}
