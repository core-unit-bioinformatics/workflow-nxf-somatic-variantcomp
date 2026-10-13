process KNOWN_PREPARE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    tuple val(meta), path(variants, stageAs: 'input?/*')
    tuple val(meta2), path(fasta), path(fai) // reference: contig names/lengths and REF alleles

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${prefix}.reference_check.tsv"), emit: report
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.known"
    """
    # list(s) -> sorted sites-only VCF on the reference contig names; contigs, positions
    # and REF alleles are checked against the reference (ext.args: --keep-mito, --warn-only)
    known_variants.py prepare \\
        --input ${[variants].flatten().join(' ')} \\
        --name ${meta.id} \\
        --fai ${fai} \\
        --fasta ${fasta} \\
        --output ${prefix}.vcf.gz \\
        --report ${prefix}.reference_check.tsv \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.known"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    printf 'key\\tvalue\\nstatus\\tok\\n' > ${prefix}.reference_check.tsv
    """
}
