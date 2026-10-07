process VARIANTCOMP_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:e96f382ca64afe53'
        : 'community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:ddf4905651532171'}"

    input:
    tuple val(meta), path(metrics), path(per_pair), path(power), path(variants)

    output:
    tuple val(meta), path("${prefix}.enrichment.tsv"), emit: enrichment
    tuple val(meta), path("${prefix}.genes.tsv"), emit: genes
    tuple val(meta), path("${prefix}.classes.tsv"), emit: classes
    tuple val(meta), path("plots/*.png"), emit: plots, optional: true
    path "*_mqc.{tsv,yaml}", emit: multiqc
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pandas'), eval("python -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas
    tuple val("${task.process}"), val('scipy'), eval("python -c 'import scipy; print(scipy.__version__)'"), topic: versions, emit: versions_scipy
    tuple val("${task.process}"), val('matplotlib'), eval("python -c 'import matplotlib; print(matplotlib.__version__)'"), topic: versions, emit: versions_matplotlib

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "variantcomp"
    def power_arg = power ? "--power ${power}" : ''
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    variantcomp_report.py \\
        --metrics ${metrics} \\
        --per-pair ${per_pair} \\
        ${power_arg} \\
        --variants ${variants} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "variantcomp"
    """
    touch ${prefix}.enrichment.tsv ${prefix}.genes.tsv ${prefix}.classes.tsv ${prefix}_metrics_mqc.tsv
    """
}
