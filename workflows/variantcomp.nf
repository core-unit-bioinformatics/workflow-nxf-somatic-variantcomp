/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PREPARE_REFERENCE      } from '../subworkflows/local/prepare_reference'
include { PREPARE_REGIONS        } from '../subworkflows/local/prepare_regions'
include { PREPROCESS_VCF         } from '../subworkflows/local/preprocess_vcf'
include { COVERAGE               } from '../subworkflows/local/coverage'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_variantcomp_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTCOMP {

    take:
    ch_samplesheet // channel: samplesheet entries [ meta, vcf, vcf_index, tumor_bam, tumor_bai, normal_bam, normal_bai, target_bed ]
    fasta          // string: reference FASTA
    fasta_fai      // string: optional reference .fai
    regions_bed    // string: optional BED restricting the comparison
    exclude_bed    // string: optional BED of regions excluded from the comparison
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // Record the validated samplesheet (resolved paths) for provenance
    //
    def samplesheet_columns = ['sample', 'technology', 'library_type', 'tumor_id', 'normal_id', 'caller', 'is_truth',
                               'vcf', 'vcf_index', 'tumor_bam', 'tumor_bai', 'normal_bam', 'normal_bai', 'target_bed']
    ch_samplesheet
        .map { entry ->
            def meta = entry[0]
            def values = ['sample', 'technology', 'library_type', 'tumor_id', 'normal_id', 'caller', 'is_truth']
                .collect { key -> meta[key] == null ? 'n/a' : meta[key].toString() }
            values += entry[1..-1].collect { path -> path ? path.toUriString() : 'n/a' }
            values.join('\t')
        }
        .collectFile(
            name: 'samplesheet_validated.tsv',
            storeDir: "${outdir}/pipeline_info",
            seed: samplesheet_columns.join('\t'),
            newLine: true,
            sort: true
        )

    //
    // SUBWORKFLOW: reference index
    //
    PREPARE_REFERENCE(
        file(fasta, checkIfExists: true),
        fasta_fai ? file(fasta_fai, checkIfExists: true) : null
    )

    //
    // SUBWORKFLOW: harmonize target BEDs and build the common comparison region
    //
    PREPARE_REGIONS(
        ch_samplesheet,
        PREPARE_REFERENCE.out.reference,
        regions_bed ? file(regions_bed, checkIfExists: true) : null,
        exclude_bed ? file(exclude_bed, checkIfExists: true) : null
    )

    //
    // SUBWORKFLOW: coverage and per-sample callable regions (alignments are optional)
    //
    COVERAGE(
        ch_samplesheet,
        PREPARE_REFERENCE.out.reference,
        PREPARE_REGIONS.out.common_bed,
        params.require_normal_coverage,
        params.callable_mode,
        params.min_coverage
    )
    ch_multiqc_files = ch_multiqc_files.mix(COVERAGE.out.multiqc_files)
    COVERAGE.out.callable_summary
        .map { _meta, tsv -> tsv }
        .collectFile(name: 'callable_summary.tsv', keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/coverage")

    //
    // SUBWORKFLOW: reference check, normalization, harmonization and uniform filtering
    //
    PREPROCESS_VCF(
        ch_samplesheet,
        PREPARE_REFERENCE.out.reference,
        PREPARE_REFERENCE.out.genome_bed,
        COVERAGE.out.regions
    )

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
