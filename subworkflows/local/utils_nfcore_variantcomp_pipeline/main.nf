//
// Subworkflow with functionality specific to the core-unit-bioinformatics/workflow-nxf-somatic-variantcomp pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    truth_technology  //  string: Technology used as ground truth
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.tsv --truth_technology <technology> --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //

    // Relative paths are resolved against the samplesheet's directory, 'n/a' marks
    // an empty optional field. Rows are checked across the whole samplesheet.
    //
    def entries = buildSamplesheetEntries(
        samplesheetToList(input, "${projectDir}/assets/schema_input.json"),
        file(input).parent,
        truth_technology
    )
    def errors = samplesheetErrors(entries, truth_technology) + missingSamplesheetFiles(entries)
    if (errors) {
        error("Invalid samplesheet '${input}':\n  - " + errors.join("\n  - "))
    }
    samplesheetWarnings(entries).each { msg -> log.warn(msg) }

    ch_samplesheet = channel.fromList(entries)

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    monochrome_logs // boolean: Disable ANSI colour codes in log output

    main:

    //
    // Completion email and summary
    //
    workflow.onComplete {

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Samplesheet handling
//
// A samplesheet entry is a list:
//   [ meta, vcf, vcf_index, tumor_bam, tumor_bai, normal_bam, normal_bai, target_bed ]
// with meta = [ id, sample, technology, library_type, tumor_id, normal_id, caller, is_truth ].
// Empty optional fields are represented by [].
//

//
// Resolve a samplesheet path: 'n/a'/empty -> [], URLs and absolute paths as given,
// relative paths against the samplesheet's directory
//
def resolveSamplesheetPath(value, base_dir) {
    def str = value instanceof List ? '' : (value == null ? '' : value.toString().trim())
    if (str == '' || str == 'n/a') {
        return []
    }
    if (str ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*/ || str.startsWith('/')) {
        return file(str)
    }
    return file("${base_dir}/${str}")
}

//
// Convert the rows returned by samplesheetToList into samplesheet entries
//
def buildSamplesheetEntries(rows, base_dir, truth_technology) {
    return rows.collect { row ->
        def meta = row[0]
        def paths = row[1..-1].collect { value -> resolveSamplesheetPath(value, base_dir) }
        def normal_id = meta.normal_id in [null, '', 'n/a'] ? null : meta.normal_id
        def new_meta = meta + [
            id       : "${meta.sample}_${meta.technology}".toString(),
            normal_id: normal_id,
            caller   : meta.caller ?: 'unknown',
            is_truth : meta.technology == truth_technology,
        ]
        [new_meta] + paths
    }
}

//
// Cross-row checks; returns a list of error messages (empty if valid)
//
def samplesheetErrors(entries, truth_technology) {
    def errors = []
    def metas = entries.collect { entry -> entry[0] }

    metas.groupBy { meta -> meta.id }.each { _id, group ->
        if (group.size() > 1) {
            errors << "Duplicate rows for sample '${group[0].sample}' and technology '${group[0].technology}'"
        }
    }

    def technologies = metas.collect { meta -> meta.technology }.unique()
    if (technologies.size() < 2) {
        errors << "At least two technologies are required, found: ${technologies.join(', ')}"
    }
    if (!truth_technology) {
        errors << "No truth technology given (--truth_technology)"
    }
    else if (!(truth_technology in technologies)) {
        errors << "Truth technology '${truth_technology}' is not in the samplesheet (technologies: ${technologies.join(', ')})"
    }
    else {
        metas.groupBy { meta -> meta.sample }.each { sample, group ->
            def sample_techs = group.collect { meta -> meta.technology }
            if (!(truth_technology in sample_techs)) {
                errors << "Sample '${sample}' has no row for the truth technology '${truth_technology}'"
            }
            else if (sample_techs.unique().size() < 2) {
                errors << "Sample '${sample}' has no technology to compare with the truth technology '${truth_technology}'"
            }
        }
    }

    entries.groupBy { entry -> entry[0].technology }.each { technology, group ->
        def beds = group.collect { entry -> entry[7] ? entry[7].toUriString() : 'n/a' }.unique()
        if (beds.size() > 1) {
            errors << "Technology '${technology}' must use one target BED for all samples, found: ${beds.join(', ')}"
        }
        def library_types = group.collect { entry -> entry[0].library_type }.unique()
        if (library_types.size() > 1) {
            errors << "Technology '${technology}' must have one library_type, found: ${library_types.join(', ')}"
        }
        if (library_types.any { type -> type != 'wgs' } && group.any { entry -> !entry[7] }) {
            errors << "Technology '${technology}' (library_type ${library_types.join('/')}) requires a target_bed"
        }
    }

    entries.each { entry ->
        def meta = entry[0]
        if (meta.normal_id && meta.tumor_id == meta.normal_id) {
            errors << "Row '${meta.id}': tumor_id and normal_id are identical ('${meta.tumor_id}')"
        }
        if (entry[4] && !entry[3]) {
            errors << "Row '${meta.id}': tumor_bai given without tumor_bam"
        }
        if (entry[6] && !entry[5]) {
            errors << "Row '${meta.id}': normal_bai given without normal_bam"
        }
        if (entry[5] && !entry[3]) {
            errors << "Row '${meta.id}': normal_bam given without tumor_bam"
        }
        if (entry[5] && !meta.normal_id) {
            errors << "Row '${meta.id}': normal_bam given but normal_id is 'n/a'"
        }
    }
    return errors.collect { msg -> msg.toString() }
}

//
// Check that all files referenced in the samplesheet exist
//
def missingSamplesheetFiles(entries) {
    def columns = ['vcf', 'vcf_index', 'tumor_bam', 'tumor_bai', 'normal_bam', 'normal_bai', 'target_bed']
    def errors = []
    entries.each { entry ->
        columns.eachWithIndex { column, idx ->
            def path = entry[idx + 1]
            if (path && !path.exists()) {
                errors << "Row '${entry[0].id}': ${column} not found: ${path.toUriString()}"
            }
        }
    }
    return errors.collect { msg -> msg.toString() }
}

//
// Non-fatal samplesheet issues
//
def samplesheetWarnings(entries) {
    def warnings = []
    def with_bam = entries.findAll { entry -> entry[3] }
    if (with_bam && with_bam.size() < entries.size()) {
        warnings << "Alignments are given for only ${with_bam.size()} of ${entries.size()} rows. Coverage-based callable regions and read-evidence checks need alignments for all technologies of a sample."
    }
    entries.findAll { entry -> !entry[0].normal_id }.each { entry ->
        warnings << "Row '${entry[0].id}': no normal sample (tumor-only). Germline variants may appear as false positives."
    }
    return warnings.collect { msg -> msg.toString() }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
