process REGION_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    // command 'targets': beds = target BEDs named <technology>.bed, other = genome BED
    // command 'callable': beds = [ callable BED ], other = []
    tuple val(meta), path(common), path(beds, stageAs: 'beds/*'), path(other, stageAs: 'other/*')
    val command
    val genome_wide // 'targets': technologies without target BED (WGS); 'callable': []

    output:
    tuple val(meta), path("${prefix}.summary.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: e.g. '--sample S1 --alignments 4 --mode base --min-coverage 10' for 'callable'
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def bed_list = [beds].flatten()
    // 'targets': one technology=BED argument per target BED; the technology is the base
    // name of the staged file ('beds/<technology>.bed') without directory and suffix
    def target_args = bed_list.collect { bed -> "${bed.toString().tokenize('/').last().replaceAll(/\.bed$/, '')}=${bed}" }.join(' ')
    def command_args = command == 'targets'
        ? "--genome ${other} ${target_args ? '--target ' + target_args : ''} ${genome_wide ? '--genome-wide ' + genome_wide.join(' ') : ''}"
        : "--callable-bed ${bed_list[0]}"
    """
    # size of the regions and their overlap with the common comparison region
    region_summary.py ${command} \\
        --common ${common} \\
        ${command_args} \\
        --output ${prefix}.summary.tsv \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.summary.tsv
    """
}
