process BENCHMARK_METRICS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    tuple val(meta), path(tp_baseline), path(fn), path(tp), path(fp), path(truth_vcf, stageAs: 'truth/*'), path(query_vcf, stageAs: 'query/*'), path(regions, stageAs: 'regions/*')

    output:
    tuple val(meta), path("${prefix}.metrics.tsv"), emit: metrics
    tuple val(meta), path("${prefix}.variants.tsv"), emit: variants
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions_arg = regions ? "--regions ${[regions].flatten().join(' ')}" : ''
    """
    benchmark_metrics.py \\
        --tp-baseline ${tp_baseline} \\
        --fn ${fn} \\
        --tp ${tp} \\
        --fp ${fp} \\
        --truth-vcf ${truth_vcf} \\
        --query-vcf ${query_vcf} \\
        --sample ${meta.sample} \\
        --truth ${meta.truth} \\
        --query ${meta.query} \\
        --prefix ${prefix} \\
        ${regions_arg} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\ttruth\\tquery\\tstratum_type\\tstratum\\ttp_baseline\\tfn\\ttp_call\\tfp\\trecall\\tprecision\\tf1\\n' > ${prefix}.metrics.tsv
    touch ${prefix}.variants.tsv
    """
}
