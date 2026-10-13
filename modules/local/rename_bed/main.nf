process RENAME_BED {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    tuple val(meta), path(rename), path(bed) // rename map of an input (input name, reference name), reference-named BED

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.renamed"
    """
    # translate the BED to the contig names of the input (e.g. an alignment with '1'
    # instead of 'chr1'), using the rename map of its reference check
    check_reference.py rename-bed \\
        --bed ${bed} \\
        --rename ${rename} \\
        --output ${prefix}.bed \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.renamed"
    """
    touch ${prefix}.bed
    """
}
