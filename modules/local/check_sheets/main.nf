process CHECK_SHEETS {
    tag "${sheet.name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    val command                    // 'samplesheet', 'stratification' or 'known'
    path sheet                     // the CSV sheet given by the user
    val base_dir                   // directory of the sheet: relative paths are resolved against it
    path samplesheet               // validated samplesheet (only for 'known'), else []

    output:
    path "${prefix}.validated.csv", emit: csv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: e.g. '--truth-technology wgs' for the samplesheet (set in conf/modules.config)
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${command}"
    def samplesheet_arg = samplesheet ? "--samplesheet ${samplesheet}" : ''
    """
    # validate the sheet; all problems are reported at once, exit 1 if invalid
    check_sheets.py ${command} \\
        --sheet ${sheet} \\
        --base-dir ${base_dir} \\
        --output ${prefix}.validated.csv \\
        ${samplesheet_arg} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${command}"
    """
    # header row only: the stub run reads it with splitCsv (no rows -> empty channels)
    head -n 1 ${sheet} > ${prefix}.validated.csv
    """
}
