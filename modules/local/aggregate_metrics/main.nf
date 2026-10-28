process AGGREGATE_METRICS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:2.2.1'
        : 'quay.io/biocontainers/pandas:2.2.1'}"

    input:
    tuple val(meta), path(metrics, stageAs: 'metrics/*')

    output:
    tuple val(meta), path("${prefix}.per_pair.tsv"), emit: per_pair
    tuple val(meta), path("${prefix}.aggregated.tsv"), emit: aggregated
    // the same tables split by matching method (<prefix>.<method>.per_pair/aggregated.tsv)
    tuple val(meta), path("${prefix}.*.{per_pair,aggregated}.tsv"), emit: by_method
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pandas'), eval("python -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # per-pair table of all inputs and micro/macro aggregation over samples
    aggregate_metrics.py \\
        metrics/* \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.per_pair.tsv ${prefix}.aggregated.tsv ${prefix}.rtgtools_vcfeval.per_pair.tsv ${prefix}.rtgtools_vcfeval.aggregated.tsv
    """
}
