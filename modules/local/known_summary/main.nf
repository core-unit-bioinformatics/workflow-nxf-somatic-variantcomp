process KNOWN_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    // per list and sample: list name, known_variants.tsv and known_summary.tsv (same order)
    tuple val(meta), val(lists), path(variants, stageAs: 'variants?/*'), path(summaries, stageAs: 'summary?/*')

    output:
    tuple val(meta), path("${prefix}.required_variants.tsv"), emit: variants
    tuple val(meta), path("${prefix}.required_summary.tsv"), emit: summary
    path "${prefix}_required_mqc.tsv", emit: multiqc
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: '--fail-on-missing' with --fail_on_missing_required (conf/modules.config)
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def variant_files = [variants].flatten()
    def summary_files = [summaries].flatten()
    // LIST=FILE arguments in the order of the lists
    def variant_args = [lists, variant_files].transpose().collect { name, file -> "${name}=${file}" }.join(' ')
    def summary_args = [lists, summary_files].transpose().collect { name, file -> "${name}=${file}" }.join(' ')
    """
    # one table of all required checks; warnings for every missing variant
    known_summary.py \\
        --variants ${variant_args} \\
        --summary ${summary_args} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.required_variants.tsv ${prefix}.required_summary.tsv ${prefix}_required_mqc.tsv
    """
}
