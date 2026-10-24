process MOSDEPTH_CALLABLE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    tuple val(meta), path(quantized), path(regions), path(rename) // mosdepth outputs and the contig rename map of the alignment
    tuple val(meta2), path(fai)
    path(clip)                                                     // BED the callable intervals are restricted to (reference names), or []

    output:
    tuple val(meta), path("${prefix}.callable.bed"), emit: bed
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: '--mode base|region --min-coverage <n>' (conf/modules.config)
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // only the parts of the callable intervals inside the clip BED are written
    def clip_arg = clip ? "--clip ${clip}" : ''
    """
    # callable intervals of the alignment, contigs renamed to the reference names,
    # restricted to the clip BED if one is given
    mosdepth_callable.py \\
        --quantized ${quantized} \\
        --regions ${regions} \\
        --rename ${rename} \\
        --fai ${fai} \\
        ${clip_arg} \\
        --output ${prefix}.callable.bed \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.callable.bed
    """
}
