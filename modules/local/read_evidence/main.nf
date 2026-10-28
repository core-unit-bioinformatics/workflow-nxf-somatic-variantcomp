process READ_EVIDENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    // one alignment (BAM/CRAM with index) and the variants of its sample
    // (CHROM POS REF ALT tables)
    tuple val(meta), path(alignment), path(index), path(sites, stageAs: 'sites/*')
    // reference FASTA with its index (and .gzi if bgzipped)
    tuple val(meta2), path(fasta), path(fai), path(gzi)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args (conf/modules/discordance.config): read filters, depth cap, indel window
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.evidence"
    """
    # reads supporting REF, ALT or another allele of every variant of the sample
    read_evidence.py \\
        --alignment ${alignment} \\
        --fasta ${fasta} \\
        --sites ${sites} \\
        --sample ${meta.id} \\
        --output ${prefix}.vcf.gz \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.evidence"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
