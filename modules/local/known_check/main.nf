process KNOWN_CHECK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    // technologies, normalized, compare, regions and evidence are lists in the same order;
    // evidence entries are [] if a technology has no tumor alignment
    tuple val(meta), path(known), val(technologies), path(normalized, stageAs: 'normalized?/*'),
        path(compare, stageAs: 'compare?/*'), path(regions, stageAs: 'regions?/*'), path(evidence, stageAs: 'evidence?/*')

    output:
    tuple val(meta), path("${prefix}.known_variants.tsv"), emit: variants
    tuple val(meta), path("${prefix}.known_summary.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def evidence_files = [evidence].flatten()
    def evidence_args = technologies.withIndex().collect { tech, _i ->
        def match = evidence_files.find { f -> f.name.startsWith("${meta.sample}_${tech}_tumor.") }
        match ? "${tech}=${match}" : null
    }.findAll { arg -> arg }
    """
    # PASS/FAIL of every required known variant in every technology of the sample
    # (compared calls, filtered calls, comparison region, tumor read evidence)
    known_variants.py check \\
        --known ${known} \\
        --sample ${meta.sample} \\
        --technology ${technologies.join(' ')} \\
        --normalized ${[normalized].flatten().join(' ')} \\
        --compare ${[compare].flatten().join(' ')} \\
        --regions ${[regions].flatten().join(' ')} \\
        ${evidence_args ? '--evidence ' + evidence_args.join(' ') : ''} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.known_variants.tsv ${prefix}.known_summary.tsv
    """
}
