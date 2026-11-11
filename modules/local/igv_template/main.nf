process IGV_TEMPLATE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:e96f382ca64afe53'
        : 'community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:ddf4905651532171'}"

    input:
    // the report template of the pipeline (assets/igv/variant_template.html) and igv.js
    // (--igv_js)
    tuple val(meta), path(template, stageAs: 'in/template.html'), path(igv_js)

    output:
    tuple val(meta), path("${prefix}.html"), emit: template
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "variant_template"
    """
    # the template with igv.js inlined instead of the CDN script tag (reports work offline)
    igv_inputs.py template \\
        --template ${template} \\
        --igv-js ${igv_js} \\
        --output ${prefix}.html \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "variant_template"
    """
    touch ${prefix}.html
    """
}
