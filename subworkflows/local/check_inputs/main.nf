//
// All input checks, before any long-running step: a wrong input stops the run in the
// first minutes, not after hours. Each check can be switched off with a --check_* param
// (conf/modules/input_checks.config); problems are then reported as warnings and the
// contig renaming is still computed.
//
//  1. sheets: --stratification and --known_variants (check_sheets.py); the samplesheet
//     itself is checked before by CHECK_SAMPLESHEET
//  2. VCFs: contig names and lengths, sample columns, REF alleles against the FASTA
//     (check_reference.py vcf)                                     --check_vcf_contigs, --check_ref_alleles
//  3. alignments: contig names and lengths of the header (check_reference.py bam)
//                                                                  --check_bam_contigs
//  4. BEDs: target BEDs, --regions_bed, --exclude_bed, --repeats_bed, stratification BEDs,
//     GIAB stratifications and known-list regions: format, contig names, interval ends
//     (check_reference.py bed)
//                                                                  --check_bed_contigs
//  5. known-variant lists: contigs, positions, REF alleles (known_variants.py prepare)
//                                                                  --check_known_variants
//  6. 'done' is emitted when all checks above have passed; the workflow starts the
//     other steps only then
//
// The mitochondrial contig is dropped from every input (reported) unless --keep_mito.
//

include { CHECK_SHEETS as CHECK_STRATIFICATION   } from '../../../modules/local/check_sheets/main'
include { CHECK_SHEETS as CHECK_KNOWN            } from '../../../modules/local/check_sheets/main'
include { CHECK_REFERENCE as CHECK_REFERENCE_VCF } from '../../../modules/local/check_reference/main'
include { CHECK_REFERENCE as CHECK_REFERENCE_BAM } from '../../../modules/local/check_reference/main'
include { CHECK_REFERENCE as CHECK_REFERENCE_BED } from '../../../modules/local/check_reference/main'
include { KNOWN_PREPARE                          } from '../../../modules/local/known_prepare/main'

workflow CHECK_INPUTS {

    take:
    ch_samplesheet // channel: samplesheet entries [ meta, vcf, vcf_index, tumor_bam, tumor_index, normal_bam, normal_index, target_bed ]
    samplesheet    // value channel: checked samplesheet CSV (sample names for the known-variant sheet)
    ch_reference   // value channel: [ meta, fasta, fai, gzi ]
    regions_bed    // path or null: --regions_bed
    exclude_bed    // path or null: --exclude_bed
    repeats_bed    // path or null: --repeats_bed (precomputed repeat regions)
    stratification // string or null: --stratification sheet (name, bed)
    giab_beds      // list: [ [ name, path ], ... ] GIAB stratifications (--giab_<name>), or []
    known_variants // string or null: --known_variants sheet (name, sample, mode, variants, regions)

    main:
    def ch_fasta_fai = ch_reference.map { meta, fasta, fai, _gzi -> [meta, fasta, fai] }

    //
    // 1. sheets; the checked sheets (absolute paths) are read into channels
    //
    // stratification BEDs: [ meta(id: name, role: 'stratum'), bed ]
    def ch_strata = channel.empty()
    def ch_sheet_checks = channel.empty()
    if (stratification) {
        def sheet = file(stratification, checkIfExists: true)
        CHECK_STRATIFICATION('stratification', sheet, sheet.parent.toUriString(), [])
        ch_strata = CHECK_STRATIFICATION.out.csv
            .splitCsv(header: true)
            .map { row -> [[id: row.name, role: 'stratum'], file(row.bed, checkIfExists: true)] }
        ch_sheet_checks = ch_sheet_checks.mix(CHECK_STRATIFICATION.out.csv)
    }
    // GIAB stratifications: named giab_<name> (column of the annotated variant table)
    ch_strata = ch_strata.mix(
        channel.fromList(giab_beds)
            .map { name, path -> [[id: "giab_${name}".toString(), role: 'stratum'], file(path, checkIfExists: true)] }
    )
    // known-variant lists: [ meta(id: name, sample, mode), [ variant files ], regions or [] ]
    def ch_known = channel.empty()
    if (known_variants) {
        def sheet = file(known_variants, checkIfExists: true)
        CHECK_KNOWN('known', sheet, sheet.parent.toUriString(), samplesheet)
        ch_known = CHECK_KNOWN.out.csv
            .splitCsv(header: true)
            .map { row ->
                [
                    [id: row.name, sample: row.sample, mode: row.mode],
                    row.variants.tokenize(';').collect { path -> file(path, checkIfExists: true) },
                    row.regions ? file(row.regions, checkIfExists: true) : [],
                ]
            }
        ch_sheet_checks = ch_sheet_checks.mix(CHECK_KNOWN.out.csv)
    }

    //
    // 2. VCFs (the index is only staged next to the VCF)
    //
    CHECK_REFERENCE_VCF(ch_samplesheet.map { entry -> [entry[0], entry[1], entry[2]] }, 'vcf', ch_fasta_fai)

    //
    // 3. alignments: [ meta(id, entry_id, role, ...), alignment, index ]
    //
    CHECK_REFERENCE_BAM(alignmentsOf(ch_samplesheet), 'bam', ch_fasta_fai)

    //
    // 4. BEDs, each with its role: [ meta(id, role), bed ]
    //
    // one target BED per technology (consistency is checked by check_sheets.py)
    def ch_beds = ch_samplesheet
        .map { entry -> [entry[0].technology, entry[7]] }
        .unique()
        .filter { _technology, bed -> bed }
        .map { technology, bed -> [[id: technology, role: 'target'], bed] }
    if (regions_bed) {
        ch_beds = ch_beds.mix(channel.of([[id: 'regions_bed', role: 'regions'], regions_bed]))
    }
    if (exclude_bed) {
        ch_beds = ch_beds.mix(channel.of([[id: 'exclude_bed', role: 'exclude'], exclude_bed]))
    }
    if (repeats_bed) {
        ch_beds = ch_beds.mix(channel.of([[id: 'repeats', role: 'repeats'], repeats_bed]))
    }
    def ch_known_regions = ch_known
        .filter { _meta, _variants, regions -> regions }
        .map { meta, _variants, regions -> [[id: meta.id, role: 'known_regions'], regions] }
    CHECK_REFERENCE_BED(
        ch_beds.mix(ch_strata, ch_known_regions).map { meta, bed -> [meta, bed, []] },
        'bed',
        ch_fasta_fai
    )
    def ch_checked_beds = CHECK_REFERENCE_BED.out.bed

    //
    // 5. known-variant lists -> sites-only VCFs on the reference (REF alleles checked)
    //
    KNOWN_PREPARE(ch_known.map { meta, variants, _regions -> [meta, variants] }, ch_fasta_fai)
    // [ meta(id, sample, mode), vcf, tbi, checked regions or [] ]
    def ch_known_lists = KNOWN_PREPARE.out.vcf
        .map { meta, vcf, tbi -> [meta.id, meta, vcf, tbi] }
        .join(
            ch_checked_beds.filter { meta, _bed -> meta.role == 'known_regions' }.map { meta, bed -> [meta.id, bed] },
            remainder: true
        )
        .map { _id, meta, vcf, tbi, regions -> [meta, vcf, tbi, regions ?: []] }

    //
    // 6. all checks passed (an empty input set passes immediately)
    //
    def ch_reports = CHECK_REFERENCE_VCF.out.report
        .mix(CHECK_REFERENCE_BAM.out.report, CHECK_REFERENCE_BED.out.report, KNOWN_PREPARE.out.report)
    def ch_done = ch_reports
        .map { _meta, report -> report }
        .mix(ch_sheet_checks)
        .collect()
        .ifEmpty([])
        .map { _files -> true }

    emit:
    vcf_rename  = CHECK_REFERENCE_VCF.out.rename                                       // channel: [ meta(entry), rename map ]
    bam_rename  = CHECK_REFERENCE_BAM.out.rename                                       // channel: [ meta(alignment), rename map ]
    beds        = ch_checked_beds                                                      // channel: [ meta(id, role), harmonized BED ]
    known_lists = ch_known_lists                                                       // channel: [ meta(id, sample, mode), vcf, tbi, regions or [] ]
    reports     = ch_reports                                                           // channel: [ meta, reference_check.tsv ]
    done        = ch_done                                                              // value channel: true when every check passed
}

//
// Alignments of the samplesheet entries: [ meta + (id: <entry>_<role>, entry_id, role), alignment, index ]
// (the same metadata is built in COVERAGE and DISCORDANCE, so that channels can be joined by meta)
//
def alignmentsOf(ch_samplesheet) {
    return ch_samplesheet.flatMap { entry ->
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
}
