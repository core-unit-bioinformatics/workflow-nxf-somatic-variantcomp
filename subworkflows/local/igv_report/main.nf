//
// IGV views of every investigated variant (--igv_report, decision D40): one igv-reports
// HTML page per sample and part of at most --igv_max_sites sites, with the tumor and
// normal alignments and the normalized calls of every technology of the sample, and an
// overview page. Sites, track configuration, template and overview are written by
// bin/igv_inputs.py; this subworkflow only connects the channels.
//

include { IGV_SITES    } from '../../../modules/local/igv_sites/main'
include { IGV_TRACKS   } from '../../../modules/local/igv_tracks/main'
include { IGV_TEMPLATE } from '../../../modules/local/igv_template/main'
include { IGV_INDEX    } from '../../../modules/local/igv_index/main'
include { IGVREPORTS   } from '../../../modules/nf-core/igvreports/main'

workflow IGV_REPORT {

    take:
    ch_samplesheet   // channel: samplesheet entries with alignment indices [ meta, vcf, vcf_index, tumor_bam, tumor_index, normal_bam, normal_index, target_bed ]
    ch_normalized    // channel: [ meta, vcf, tbi ] normalized records of every entry (PASS and non-PASS)
    ch_variants      // channel (one item): path of every compared variant (variants.annotated.tsv)
    ch_required      // channel (one item): path of the required known variants, or [] if none
    ch_sheet         // channel (one item): validated samplesheet (CSV)
    ch_reference     // channel (value): [ meta, fasta, fai, gzi ]
    template         // path: report template (assets/igv/variant_template.html)
    igv_js           // path: igv.js (--igv_js)

    main:
    //
    // MODULE: sites tables per sample and part, links table variant_id -> report file
    //
    IGV_SITES(
        ch_variants
            // wrapped, so that combine keeps [] (no required variants) as one element
            .combine(ch_required.map { tsv -> [tsv] })
            .combine(ch_reference.map { _meta, _fasta, fai, _gzi -> fai })
            .map { variants, required, fai -> [[id: 'igv'], variants, required, fai] }
    )

    // sites tables as [ sample, part id (<sample>_<part>), table ]; the folder of a table
    // is its sample (sites/<sample>/<sample>_<part>.sites.tsv)
    def ch_parts = IGV_SITES.out.sites
        .flatMap { _meta, tables -> [tables].flatten() }
        .map { table -> [table.parent.name, table.name - '.sites.tsv', table] }

    // track files of every sample: [ sample, label, file, index ] per alignment and
    // normalized VCF; tumor-only rows have no normal alignment
    def ch_track_files = ch_samplesheet
        .map { entry -> [entry[0].id, entry[0], entry[3], entry[4], entry[5], entry[6]] }
        .join(ch_normalized.map { meta, vcf, tbi -> [meta.id, vcf, tbi] }, failOnMismatch: true, failOnDuplicate: true)
        .flatMap { _id, meta, tumor, tumor_index, normal, normal_index, vcf, tbi ->
            [['tumor', tumor, tumor_index], ['normal', normal, normal_index], ['calls', vcf, tbi]]
                .findAll { item -> item[1] }
                .collect { item -> [meta.sample, "${meta.technology}:${item[0]}".toString(), item[1], item[2]] }
        }
    // per sample: [ meta, labels, files, indices ] in a fixed order (sample, label), so that
    // the staged names (tracks/<n>/) and the task hashes don't depend on the arrival order
    def ch_sample_tracks = ch_track_files
        .toSortedList { a, b -> a[0] <=> b[0] ?: a[1] <=> b[1] }
        .flatMap { items -> items }
        .groupTuple()
        .map { sample, labels, files, indices -> [[id: sample, sample: sample], labels, files, indices] }

    //
    // MODULE: track configuration of every sample (names, order, display filters) and the
    // SAM flags excluded when the reads are embedded
    //
    IGV_TRACKS(
        ch_sample_tracks.map { meta, labels, files, _indices -> [meta, labels, files] },
        ch_sheet
    )

    //
    // MODULE: report template with igv.js inlined (offline)
    //
    IGV_TEMPLATE([[id: 'igv'], template, igv_js])

    //
    // MODULE: igv-reports: one HTML page per sample and part
    //
    // [ sample, meta, files, indices, config, flags ] per sample
    def ch_sample_inputs = ch_sample_tracks
        .map { meta, _labels, files, indices -> [meta.sample, meta, files, indices] }
        .join(IGV_TRACKS.out.config.map { meta, config, flags -> [meta.sample, config, flags] }, failOnMismatch: true)
    IGVREPORTS(
        ch_parts
            .combine(ch_sample_inputs, by: 0)
            .map { sample, part, table, _meta, files, indices, config, flags ->
                [[id: part, sample: sample, exclude_flags: flags], table, files, indices, config]
            },
        ch_reference.map { meta, fasta, fai, _gzi -> [meta, fasta, fai] },
        // a value channel (IGV_TEMPLATE runs once on value inputs): used by every page
        IGV_TEMPLATE.out.template.map { _meta, html -> html }
    )

    //
    // MODULE: overview page report/igv/index.html
    //
    IGV_INDEX(IGV_SITES.out.links)

    emit:
    links   = IGV_SITES.out.links      // channel: [ meta, igv_links.tsv ]
    reports = IGVREPORTS.out.report    // channel: [ meta, <sample>_<part>_report.html ]
    index   = IGV_INDEX.out.index      // channel: [ meta, index.html ]
}
