process CHECK_REFERENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    tuple val(meta), path(input)
    val type // 'vcf', 'bed' or 'bam'
    tuple val(meta2), path(fasta), path(fai), path(gzi)
    path regions // optional (bam): BED in reference names, translated to the alignment's names

    output:
    tuple val(meta), path("*.reference_check.tsv"), emit: report
    tuple val(meta), path("*.rename_chrs.tsv"), emit: rename
    tuple val(meta), path("*.keep_contigs.txt"), emit: keep
    tuple val(meta), path("*.harmonized.bed"), emit: bed, optional: true
    tuple val(meta), path("*.regions.bed"), emit: regions, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${type}"
    def regions_arg = regions ? "--regions ${regions}" : ''
    """
    check_reference.py \\
        --input ${input} \\
        --type ${type} \\
        --fai ${fai} \\
        --fasta ${fasta} \\
        --prefix ${prefix} \\
        ${regions_arg} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${type}"
    """
    printf 'key\\tvalue\\nstatus\\tok\\n' > ${prefix}.reference_check.tsv
    touch ${prefix}.rename_chrs.tsv ${prefix}.keep_contigs.txt
    if [ "${type}" == "bed" ]; then touch ${prefix}.harmonized.bed; fi
    if [ "${type}" == "bam" ] && [ -n "${regions}" ]; then touch ${prefix}.regions.bed; fi
    """
}
