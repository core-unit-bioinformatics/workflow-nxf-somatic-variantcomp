process CHECKS_GATE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:e96f382ca64afe53'
        : 'community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:ddf4905651532171'}"

    input:
    // the consistency checks; the reports are only staged, so that the gate runs after
    // they were written
    tuple val(meta), path(checks), path(reports, stageAs: 'reports/*')

    output:
    // versions are not sent to the `versions` topic: this process waits for the collated
    // versions (directly or through the reports), so it would keep the topic open forever
    tuple val(meta), path(checks), emit: checks
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # stop the run (exit status 1) if a consistency check failed
    consistency_checks.py gate \\
        --checks ${checks} \\
        ${args}
    """

    stub:
    """
    true
    """
}
