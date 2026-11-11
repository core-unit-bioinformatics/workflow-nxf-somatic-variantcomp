process IGV_INDEX {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:e96f382ca64afe53'
        : 'community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:ddf4905651532171'}"

    input:
    // links table of IGV_SITES (variant_id, sample, part, igv_file)
    tuple val(meta), path(links)

    output:
    tuple val(meta), path("index.html"), emit: index
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args (conf/modules/igv_report.config): --max-sites (shown), links to the reports
    def args = task.ext.args ?: ''
    """
    # overview page of the IGV views: every sample and part with its number of sites
    igv_inputs.py index \\
        --links ${links} \\
        --output index.html \\
        ${args}
    """

    stub:
    """
    touch index.html
    """
}
