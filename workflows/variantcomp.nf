/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { CHECK_SHEETS as CHECK_SAMPLESHEET } from '../modules/local/check_sheets/main'
include { PREPARE_REFERENCE      } from '../subworkflows/local/prepare_reference'
include { CHECK_INPUTS           } from '../subworkflows/local/check_inputs'
include { PREPARE_REGIONS        } from '../subworkflows/local/prepare_regions'
include { PREPARE_ALIGNMENTS     } from '../subworkflows/local/prepare_alignments'
include { PREPROCESS_VCF         } from '../subworkflows/local/preprocess_vcf'
include { COVERAGE               } from '../subworkflows/local/coverage'
include { BENCHMARK              } from '../subworkflows/local/benchmark'
include { DISCORDANCE            } from '../subworkflows/local/discordance'
include { ANNOTATE_VARIANTS      } from '../subworkflows/local/annotate_variants'
include { VARIANTCOMP_REPORT     } from '../modules/local/variantcomp_report/main'
include { MQC_TABLES             } from '../modules/local/mqc_tables/main'
include { PREPARE_KNOWN_VARIANTS } from '../subworkflows/local/prepare_known_variants'
include { KNOWN_VARIANTS         } from '../subworkflows/local/known_variants'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_variantcomp_pipeline'
include { giabStratifications    } from '../subworkflows/local/utils_nfcore_variantcomp_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTCOMP {

    take:
    ch_samplesheet_input // channel: samplesheet entries [ meta, vcf, vcf_index, tumor_bam, tumor_index, normal_bam, normal_index, target_bed ]
    input          // string: samplesheet (--input), for the checks across rows
    stratification // string: optional stratification sheet (--stratification: name, bed)
    known_variants // string: optional known-variant sheet (--known_variants: name, sample, mode, variants, regions)
    fasta          // string: reference FASTA
    fasta_fai      // string: optional reference .fai
    regions_bed    // string: optional BED restricting the comparison
    exclude_bed    // string: optional BED of regions excluded from the comparison
    repeats_bed    // string: optional precomputed repeat BED (contig, start, end, unit, times)
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    // alignments are only used for callable regions and read evidence; without both,
    // they are ignored (not checked, not indexed, not read)
    def use_alignments = params.callable_regions || params.read_evidence
    def ch_entries = use_alignments
        ? ch_samplesheet_input
        : ch_samplesheet_input.map { entry -> [entry[0], entry[1], entry[2], [], [], [], [], entry[7]] }

    //
    // MODULE: checks of the samplesheet across rows (bin/check_sheets.py: truth technology
    // per sample, unique rows, one target BED per technology, ...); all problems are
    // reported at once. The columns were already checked by nf-schema.
    //
    def input_file = file(input, checkIfExists: true)
    CHECK_SAMPLESHEET('samplesheet', input_file, input_file.parent.toUriString(), [])
    def ch_samplesheet_validated = CHECK_SAMPLESHEET.out.csv.first()

    //
    // SUBWORKFLOW: reference index and genome BED
    //
    PREPARE_REFERENCE(
        file(fasta, checkIfExists: true),
        fasta_fai ? file(fasta_fai, checkIfExists: true) : null
    )

    //
    // SUBWORKFLOW: all other input checks (sheets, VCFs, alignments, BEDs, known lists)
    // against the reference, before any long-running step
    //
    CHECK_INPUTS(
        ch_entries,
        ch_samplesheet_validated,
        PREPARE_REFERENCE.out.reference,
        regions_bed ? file(regions_bed, checkIfExists: true) : null,
        exclude_bed ? file(exclude_bed, checkIfExists: true) : null,
        repeats_bed ? file(repeats_bed, checkIfExists: true) : null,
        stratification,
        // GIAB genome stratifications given as --giab_* file params: [ [ name, path ], ... ]
        giabStratifications(),
        known_variants
    )

    //
    // SUBWORKFLOW: alignment indices (created where missing); this and every following
    // step starts only after all checks passed
    //
    PREPARE_ALIGNMENTS(
        ch_entries
            .combine(ch_samplesheet_validated)
            .combine(CHECK_INPUTS.out.done)
            .map { entry -> entry[0..-3] }
    )
    def ch_samplesheet = PREPARE_ALIGNMENTS.out.samplesheet

    //
    // SUBWORKFLOW: harmonize target BEDs and build the common comparison region
    //
    PREPARE_REGIONS(
        ch_samplesheet,
        PREPARE_REFERENCE.out.reference,
        PREPARE_REFERENCE.out.genome_bed,
        CHECK_INPUTS.out.beds,
        regions_bed as boolean,
        exclude_bed as boolean
    )

    //
    // SUBWORKFLOW: coverage and per-sample callable regions (--callable_regions); without,
    // every entry is compared in the common region
    //
    // comparison region per samplesheet entry: [ meta, bed ]
    def ch_comparison_regions = ch_samplesheet.map { entry -> entry[0] }.combine(PREPARE_REGIONS.out.common_bed)
    if (params.callable_regions) {
        COVERAGE(
            ch_samplesheet,
            CHECK_INPUTS.out.bam_rename,
            PREPARE_REFERENCE.out.reference,
            PREPARE_REGIONS.out.common_bed,
            params.require_normal_coverage
        )
        ch_comparison_regions = COVERAGE.out.regions
        ch_multiqc_files = ch_multiqc_files.mix(COVERAGE.out.multiqc_files)
        COVERAGE.out.callable_summary
            .map { _meta, tsv -> tsv }
            .collectFile(name: 'callable_summary.tsv', keepHeader: true, skip: 1, sort: { tsv -> tsv.name }, storeDir: "${outdir}/coverage")
    }

    //
    // SUBWORKFLOW: normalization, harmonization and uniform filtering
    //
    PREPROCESS_VCF(
        ch_samplesheet,
        CHECK_INPUTS.out.vcf_rename,
        PREPARE_REFERENCE.out.reference,
        PREPARE_REFERENCE.out.genome_bed,
        ch_comparison_regions,
        params.filter_pass_only
    )

    //
    // SUBWORKFLOW: compare every technology with the truth technology (rtg vcfeval) and compute metrics
    //
    BENCHMARK(
        PREPROCESS_VCF.out.compare,
        ch_comparison_regions,
        PREPARE_REFERENCE.out.reference,
        params.comparison_mode
    )

    //
    // SUBWORKFLOW: known variant lists on the reference, per sample
    //
    PREPARE_KNOWN_VARIANTS(
        CHECK_INPUTS.out.known_lists.combine(CHECK_INPUTS.out.done).map { entry -> entry[0..-2] },
        ch_samplesheet,
        PREPARE_REFERENCE.out.reference
    )

    //
    // SUBWORKFLOW: read evidence and classification of discordant variants
    //
    DISCORDANCE(
        ch_samplesheet,
        PREPROCESS_VCF.out.normalized,
        PREPROCESS_VCF.out.compare,
        BENCHMARK.out.variants,
        PREPARE_REFERENCE.out.reference,
        PREPARE_KNOWN_VARIANTS.out.sites,
        params.read_evidence
    )

    //
    // SUBWORKFLOW: known variants ('required' check, 'benchmark' against an external truth)
    //
    KNOWN_VARIANTS(
        PREPARE_KNOWN_VARIANTS.out.per_sample,
        PREPROCESS_VCF.out.normalized,
        PREPROCESS_VCF.out.compare,
        ch_comparison_regions,
        DISCORDANCE.out.evidence,
        BENCHMARK.out.sdf
    )
    ch_multiqc_files = ch_multiqc_files.mix(KNOWN_VARIANTS.out.multiqc)

    //
    // SUBWORKFLOW: VEP consequences, repeat regions and stratifications of every
    // compared variant
    //
    ANNOTATE_VARIANTS(
        DISCORDANCE.out.variants,
        PREPARE_REFERENCE.out.reference,
        CHECK_INPUTS.out.beds.filter { meta, _bed -> meta.role == 'stratum' },
        CHECK_INPUTS.out.beds.filter { meta, _bed -> meta.role == 'repeats' },
        params.vep,
        params.vep_cache ? file(params.vep_cache, checkIfExists: true) : null,
        params.vep_gtf ? file(params.vep_gtf, checkIfExists: true) : null,
        params.annotate_repeats,
        repeats_bed as boolean
    )
    ch_multiqc_files = ch_multiqc_files.mix(ANNOTATE_VARIANTS.out.multiqc)

    //
    // MODULE: statistics (enrichment, genes, classes), plots and MultiQC sections
    //
    VARIANTCOMP_REPORT(
        BENCHMARK.out.aggregated.map { _meta, tsv -> tsv }
            .combine(BENCHMARK.out.per_pair.map { _meta, tsv -> tsv })
            .combine(DISCORDANCE.out.power.map { _meta, tsv -> tsv })
            .combine(ANNOTATE_VARIANTS.out.annotated.map { _meta, tsv -> tsv })
            .map { metrics, per_pair, power, variants -> [[id: 'variantcomp'], metrics, per_pair, power, variants] }
    )
    ch_multiqc_files = ch_multiqc_files.mix(VARIANTCOMP_REPORT.out.multiqc.flatten())

    //
    // MODULE: MultiQC tables of the step summaries (input checks, normalization,
    // harmonization, regions), one section each
    //
    def ch_summary_sections = CHECK_INPUTS.out.reports
        .map { _meta, tsv -> tsv }
        .collect()
        .map { files -> [[id: 'input_checks', command: 'key-value', title: 'Input checks', description: 'Reference checks of every input at the start of the run (status, renamed and dropped contigs, REF mismatches; input_checks/).'], files] }
        .mix(
            PREPROCESS_VCF.out.norm_logs
                .map { _meta, log -> log }
                .collect()
                .map { files -> [[id: 'normalization', command: 'norm-logs', title: 'Normalization (bcftools norm)', description: 'Records split, left-aligned, removed (REF mismatch, duplicates) per VCF and branch (pass / non_pass).'], files] },
            PREPROCESS_VCF.out.summaries
                .map { _meta, tsv -> tsv }
                .collect()
                .map { files -> [[id: 'harmonization', command: 'key-value', title: 'Harmonization (VAF, depth, alt reads)', description: 'Source of the VAF (FORMAT/AF or AD), depth and dropped records per VCF and branch (harmonize_vcf.py).'], files] },
            PREPARE_REGIONS.out.summary
                .map { _meta, tsv -> tsv }
                .collect()
                .map { files -> [[id: 'regions', command: 'table', title: 'Comparison region', description: 'Size of every target BED and of the common comparison region.'], files] }
        )
    if (params.callable_regions) {
        ch_summary_sections = ch_summary_sections.mix(
            COVERAGE.out.callable_summary
                .map { _meta, tsv -> tsv }
                .collect()
                .map { files -> [[id: 'callable', command: 'table', title: 'Callable regions', description: 'Bases of the common region callable in all alignments of a sample.'], files] }
        )
    }
    MQC_TABLES(ch_summary_sections)
    ch_multiqc_files = ch_multiqc_files.mix(MQC_TABLES.out.multiqc)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'variantcomp_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'variantcomp'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
