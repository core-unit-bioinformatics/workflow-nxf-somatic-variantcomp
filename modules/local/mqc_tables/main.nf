process MQC_TABLES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    // meta: id (MultiQC section id), command (key-value, norm-logs, table), title,
    // description; summaries: the files of the section (staged in in1/, in2/, ... so that
    // equal names don't collide)
    tuple val(meta), path(summaries, stageAs: 'in?/*')

    output:
    path "${meta.id}_mqc.tsv", emit: multiqc
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: e.g. '--strip .reference_check.tsv' or '--id-columns 2'
    def args = task.ext.args ?: ''
    """
    # one MultiQC table of all summaries of this section
    mqc_tables.py ${meta.command} \\
        ${[summaries].flatten().join(' ')} \\
        --section ${meta.id} \\
        --title '${meta.title}' \\
        --description '${meta.description}' \\
        --prefix ${meta.id} \\
        ${args}
    """

    stub:
    """
    touch ${meta.id}_mqc.tsv
    """
}
