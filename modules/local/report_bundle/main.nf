process REPORT_BUNDLE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:e96f382ca64afe53'
        : 'community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:ddf4905651532171'}"

    input:
    // the MultiQC report and the datavzrd output folder
    tuple val(meta), path(multiqc, stageAs: 'in/multiqc_report.html'), path(tables, stageAs: 'in/tables')

    output:
    // versions are not sent to the `versions` topic: this process waits for the MultiQC
    // report, which waits for the collated versions, so it would keep the topic open forever
    tuple val(meta), path("${prefix}.zip"), emit: zip
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "variantcomp_report"
    """
    # one zip with the MultiQC report and the table browser in the layout of the output
    # folder (the reports link to each other with relative paths)
    report_bundle.py \\
        --multiqc ${multiqc} \\
        --tables ${tables} \\
        --name ${prefix} \\
        --output ${prefix}.zip \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "variantcomp_report"
    """
    touch ${prefix}.zip
    """
}
