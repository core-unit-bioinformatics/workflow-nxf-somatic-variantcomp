/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PREPARE_REFERENCE      } from '../subworkflows/local/prepare_reference'
include { PREPARE_REGIONS        } from '../subworkflows/local/prepare_regions'
include { PREPARE_ALIGNMENTS     } from '../subworkflows/local/prepare_alignments'
include { PREPROCESS_VCF         } from '../subworkflows/local/preprocess_vcf'
include { COVERAGE               } from '../subworkflows/local/coverage'
include { BENCHMARK              } from '../subworkflows/local/benchmark'
include { DISCORDANCE            } from '../subworkflows/local/discordance'
include { VARIANTCOMP_REPORT     } from '../modules/local/variantcomp_report/main'
include { PREPARE_KNOWN_VARIANTS } from '../subworkflows/local/prepare_known_variants'
include { KNOWN_VARIANTS         } from '../subworkflows/local/known_variants'
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
    ch_samplesheet_input // channel: samplesheet entries [ meta, vcf, vcf_index, tumor_bam, tumor_index, normal_bam, normal_index, target_bed ]
    ch_strata      // channel: stratification BEDs [ meta(id: name), bed ]
    ch_known       // channel: known variant lists [ meta(id: name, sample, mode), [ files ], regions ]
    fasta          // string: reference FASTA
    fasta_fai      // string: optional reference .fai
    regions_bed    // string: optional BED restricting the comparison
    exclude_bed    // string: optional BED of regions excluded from the comparison
    gene_bed       // string: optional gene BED (BED4)
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // SUBWORKFLOW: alignment indices (created where missing)
    //
    PREPARE_ALIGNMENTS(ch_samplesheet_input)
    def ch_samplesheet = PREPARE_ALIGNMENTS.out.samplesheet

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
    // SUBWORKFLOW: compare every technology with the truth technology (rtg vcfeval) and compute metrics
    //
    BENCHMARK(
        PREPROCESS_VCF.out.compare,
        COVERAGE.out.regions,
        PREPARE_REFERENCE.out.reference,
        params.comparison_mode
    )

    //
    // SUBWORKFLOW: known variant lists on the reference, per sample
    //
    PREPARE_KNOWN_VARIANTS(ch_known, ch_samplesheet, PREPARE_REFERENCE.out.reference)

    //
    // SUBWORKFLOW: read evidence, classification and annotation of discordant variants
    //
    DISCORDANCE(
        ch_samplesheet,
        PREPROCESS_VCF.out.normalized,
        PREPROCESS_VCF.out.compare,
        BENCHMARK.out.variants,
        PREPARE_REFERENCE.out.reference,
        gene_bed ? file(gene_bed, checkIfExists: true) : null,
        ch_strata,
        PREPARE_KNOWN_VARIANTS.out.sites
    )

    //
    // SUBWORKFLOW: known variants ('required' check, 'benchmark' against an external truth)
    //
    KNOWN_VARIANTS(
        PREPARE_KNOWN_VARIANTS.out.per_sample,
        PREPROCESS_VCF.out.normalized,
        PREPROCESS_VCF.out.compare,
        COVERAGE.out.regions,
        DISCORDANCE.out.evidence,
        BENCHMARK.out.sdf
    )
    def ch_required = KNOWN_VARIANTS.out.required_variants
        .map { _meta, tsv -> tsv }
        .collectFile(name: 'required_variants.tsv', keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/known_variants")
    def ch_required_summary = KNOWN_VARIANTS.out.required_summary
        .map { meta, tsv -> tsv.readLines().drop(1).collect { line -> "${meta.name}\t${line}" }.join('\n') }
        .collectFile(
            name: 'required_summary.tsv', sort: true, newLine: true, storeDir: "${outdir}/known_variants",
            seed: 'list\tsample\ttechnology\trequired\tpassed\tfailed\tall_passed'
        )
    ch_required
        .map { tsv -> tsv.readLines() }
        .subscribe { lines ->
            def header = lines[0].tokenize('\t')
            def failed = lines.drop(1).collect { line -> [header, line.split('\t', -1).toList()].transpose().collectEntries() }
                .findAll { row -> row.status != 'PASS' }
            failed.each { row ->
                log.warn("Required variant ${row.label} ${row.chrom}:${row.pos} ${row.ref}>${row.alt} (${row.gene}) not found in ${row.sample} ${row.technology}: ${row.status}")
            }
            if (failed && params.fail_on_missing_required) {
                error("${failed.size()} required variant(s) not found, see ${outdir}/known_variants/required_variants.tsv (--fail_on_missing_required)")
            }
        }
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_required_summary.map { tsv ->
            def header = "# id: 'variantcomp_required'\n# section_name: 'Required variants'\n" +
                "# description: 'Known variants that must be found in every technology (--known_variants, mode required).'\n" +
                "# plot_type: 'table'\n"
            def lines = tsv.readLines()
            def table = lines.withIndex().collect { line, i -> i == 0 ? 'id\t' + line : "${line.tokenize('\t')[0..2].join(' ')}\t${line}" }
            def mqc = file("${workflow.workDir}/tmp/variantcomp_required_mqc.tsv")
            mqc.parent.mkdirs()
            mqc.text = header + table.join('\n') + '\n'
            mqc
        }
    )


    //
    // MODULE: statistics (enrichment, genes, classes), plots and MultiQC sections
    //
    VARIANTCOMP_REPORT(
        BENCHMARK.out.aggregated.map { _meta, tsv -> tsv }
            .combine(BENCHMARK.out.per_pair.map { _meta, tsv -> tsv })
            .combine(DISCORDANCE.out.power.map { _meta, tsv -> tsv })
            .combine(DISCORDANCE.out.annotated.map { _meta, tsv -> tsv })
            .map { metrics, per_pair, power, variants -> [[id: 'variantcomp'], metrics, per_pair, power, variants] }
    )
    ch_multiqc_files = ch_multiqc_files.mix(VARIANTCOMP_REPORT.out.multiqc.flatten())

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
