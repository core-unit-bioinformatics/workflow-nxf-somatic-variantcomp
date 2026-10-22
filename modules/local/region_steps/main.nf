process REGION_STEPS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    // the region after every step of PREPARE_REGIONS, in step order; names: one per BED
    tuple val(meta), val(names), path(beds, stageAs: 'steps/?/*')

    output:
    tuple val(meta), path("${prefix}.steps.tsv"), emit: steps
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // one NAME=BED argument per step, in step order
    def step_args = [names, [beds].flatten()].transpose().collect { name, bed -> "${name}=${bed}" }.join(' ')
    """
    # bp of the comparison region after every step (genome, intersection of the target
    # BEDs, --regions_bed, --exclude_bed)
    region_summary.py steps \\
        --step ${step_args} \\
        --output ${prefix}.steps.tsv \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'step\\tname\\tintervals\\tbp\\tfraction_of_genome\\tfraction_of_previous\\tremoved_bp\\n' > ${prefix}.steps.tsv
    """
}
