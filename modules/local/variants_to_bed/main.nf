process VARIANTS_TO_BED {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    tuple val(meta), path(table) // variant table (TSV with chrom, pos, ref)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.variants"
    """
    # REF span of every variant as BED4, row number as name
    annotate_variants.py to-bed \\
        --table ${table} \\
        --output ${prefix}.bed \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.variants"
    """
    touch ${prefix}.bed
    """
}
