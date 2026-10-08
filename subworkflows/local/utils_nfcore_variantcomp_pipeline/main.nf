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
include { paramsHelp                } from 'plugin/nf-schema'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { CHECK_SHEETS as CHECK_SAMPLESHEET    } from '../../../modules/local/check_sheets/main'
include { CHECK_SHEETS as CHECK_STRATIFICATION } from '../../../modules/local/check_sheets/main'
include { CHECK_SHEETS as CHECK_KNOWN          } from '../../../modules/local/check_sheets/main'

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
    stratification    //  string: Optional stratification sheet (name, bed)
    known_variants    //  string: Optional known-variant sheet (name, sample, mode, variants, regions)
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

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --fasta genome.fa --truth_technology <technology> --outdir <OUTDIR>"

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
    // Sheets: content checks in bin/check_sheets.py (all problems reported at once);
    // here only the validated CSVs are read into channels.
    //
    def input_file = file(input, checkIfExists: true)
    CHECK_SAMPLESHEET('samplesheet', input_file, input_file.parent.toUriString(), [])
    // samplesheet entry: [ meta, vcf, vcf_index, tumor_bam, tumor_index, normal_bam, normal_index, target_bed ]
    ch_samplesheet = CHECK_SAMPLESHEET.out.csv
        .splitCsv(header: true)
        .map { row -> samplesheetEntry(row) }

    def ch_strata = channel.empty()
    if (stratification) {
        def strat_file = file(stratification, checkIfExists: true)
        CHECK_STRATIFICATION('stratification', strat_file, strat_file.parent.toUriString(), [])
        // stratum: [ meta(id: name), bed ]
        ch_strata = CHECK_STRATIFICATION.out.csv
            .splitCsv(header: true)
            .map { row -> [[id: row.name], file(row.bed, checkIfExists: true)] }
    }

    def ch_known = channel.empty()
    if (known_variants) {
        def known_file = file(known_variants, checkIfExists: true)
        CHECK_KNOWN('known', known_file, known_file.parent.toUriString(), CHECK_SAMPLESHEET.out.csv)
        // known-variant list: [ meta(id: name, sample, mode), [ variant files ], regions or [] ]
        ch_known = CHECK_KNOWN.out.csv
            .splitCsv(header: true)
            .map { row ->
                [
                    [id: row.name, sample: row.sample, mode: row.mode],
                    row.variants.tokenize(';').collect { path -> file(path, checkIfExists: true) },
                    row.regions ? file(row.regions, checkIfExists: true) : [],
                ]
            }
    }

    emit:
    samplesheet = ch_samplesheet
    strata      = ch_strata
    known       = ch_known
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
// Samplesheet entry from a row of the validated samplesheet (bin/check_sheets.py)
//
// Content was checked by the script; here the files are only staged (checkIfExists
// fails with the missing path) and index files next to them are looked up.
//
def samplesheetEntry(row) {
    def meta = [
        id          : "${row.sample}_${row.technology}".toString(),
        sample      : row.sample,
        technology  : row.technology,
        library_type: row.library_type,
        tumor_id    : row.tumor_id,
        normal_id   : row.normal_id ?: null,
        caller      : row.caller,
        is_truth    : row.is_truth == 'true',
    ]
    def vcf = file(row.vcf, checkIfExists: true)
    def tumor_bam = row.tumor_bam ? file(row.tumor_bam, checkIfExists: true) : []
    def normal_bam = row.normal_bam ? file(row.normal_bam, checkIfExists: true) : []
    return [
        meta,
        vcf,
        findIndex(vcf, ['.tbi', '.csi'], []),
        tumor_bam,
        tumor_bam ? findIndex(tumor_bam, ['.bai', '.csi', '.crai'], [/\.bam$/: '.bai', /\.cram$/: '.crai']) : [],
        normal_bam,
        normal_bam ? findIndex(normal_bam, ['.bai', '.csi', '.crai'], [/\.bam$/: '.bai', /\.cram$/: '.crai']) : [],
        row.target_bed ? file(row.target_bed, checkIfExists: true) : [],
    ]
}

//
// Index file next to a data file: <file><suffix> (e.g. x.bam.bai) or, with a
// replacement, <file without extension><new extension> (e.g. x.bai); [] if none exists
// (the pipeline creates missing indices where they are needed)
//
def findIndex(data_file, suffixes, replacements) {
    def name = data_file.toUriString()
    def candidates = suffixes.collect { suffix -> name + suffix }
    replacements.each { pattern, extension -> candidates << name.replaceAll(pattern, extension) }
    def found = candidates.findAll { candidate -> candidate != name }.collect { candidate -> file(candidate) }.find { candidate -> candidate.exists() }
    return found ?: []
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // in-text citations of every tool of the pipeline (full list: CITATIONS.md)
    def citation_text = [
            "Tools used in the workflow included:",
            "BCFtools, SAMtools and HTSlib (Danecek et al. 2021) for VCF and alignment processing,",
            "BEDTools (Quinlan and Hall 2010) for region operations,",
            "mosdepth (Pedersen and Quinlan 2018) for coverage (if alignments were given),",
            "RTG Tools vcfeval (Cleary et al. 2015) for haplotype-aware variant matching,",
            "GNU Awk, pysam, pandas (McKinney 2010), NumPy (Harris et al. 2020), SciPy (Virtanen et al. 2020) and Matplotlib (Hunter 2007) in the pipeline scripts,",
            "the Benjamini-Hochberg procedure (Benjamini and Hochberg 1995) for multiple-testing correction",
            "and MultiQC (Ewels et al. 2016) for the report."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // bibliography entries matching toolCitationText()
    def reference_text = [
            "<li>Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. doi: 10.1093/gigascience/giab008</li>",
            "<li>Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841–842. doi: 10.1093/bioinformatics/btq033</li>",
            "<li>Pedersen, B. S., & Quinlan, A. R. (2018). Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics, 34(5), 867–868. doi: 10.1093/bioinformatics/btx699</li>",
            "<li>Cleary, J. G., Braithwaite, R., Gaastra, K., Hilbush, B. S., Inglis, S., Irvine, S. A., Jackson, A., Littin, R., Rathod, M., Ware, D., Zook, J. M., Trigg, L., & De La Vega, F. M. (2015). Comparing variant call files for performance benchmarking of next-generation sequencing variant calling pipelines. bioRxiv. doi: 10.1101/023754</li>",
            "<li>McKinney, W. (2010). Data structures for statistical computing in Python. Proceedings of the 9th Python in Science Conference, 56–61. doi: 10.25080/Majora-92bf1922-00a</li>",
            "<li>Harris, C. R., Millman, K. J., van der Walt, S. J., et al. (2020). Array programming with NumPy. Nature, 585, 357–362. doi: 10.1038/s41586-020-2649-2</li>",
            "<li>Virtanen, P., Gommers, R., Oliphant, T. E., et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods, 17, 261–272. doi: 10.1038/s41592-019-0686-2</li>",
            "<li>Hunter, J. D. (2007). Matplotlib: a 2D graphics environment. Computing in Science & Engineering, 9(3), 90–95. doi: 10.1109/MCSE.2007.55</li>",
            "<li>Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society: Series B, 57(1), 289–300. doi: 10.1111/j.2517-6161.1995.tb02031.x</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. doi: 10.1093/bioinformatics/btw354</li>"
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

    meta["tool_citations"] = toolCitationText()
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
