process CLASSIFY_DISCORDANCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    tuple val(meta), path(variants), path(truth_normalized, stageAs: 'truth/*'), path(query_normalized, stageAs: 'query/*'),
        path(truth_tumor_evidence, stageAs: 'evidence_truth/*'), path(query_tumor_evidence, stageAs: 'evidence_query/*'),
        path(normal_evidence, stageAs: 'evidence_normal/*')

    output:
    tuple val(meta), path("${prefix}.discordance.tsv"), emit: discordance
    tuple val(meta), path("${prefix}.power.tsv"), emit: power
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def truth_tumor = truth_tumor_evidence ? "--truth-tumor-evidence ${truth_tumor_evidence}" : ''
    def query_tumor = query_tumor_evidence ? "--query-tumor-evidence ${query_tumor_evidence}" : ''
    def normals = normal_evidence ? "--normal-evidence ${normal_evidence.join(' ')}" : ''
    """
    classify_discordance.py \\
        --variants ${variants} \\
        --truth-normalized ${truth_normalized} \\
        --query-normalized ${query_normalized} \\
        ${truth_tumor} \\
        ${query_tumor} \\
        ${normals} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.discordance.tsv ${prefix}.power.tsv
    """
}
