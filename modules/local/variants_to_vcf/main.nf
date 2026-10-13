process VARIANTS_TO_VCF {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    tuple val(meta), path(table)  // variant table (TSV with chrom, pos, ref, alt)
    tuple val(meta2), path(fai)   // reference index (contig order of the VCF)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.sites"
    """
    # distinct variants of the table as a sorted sites-only VCF (input of VEP)
    annotate_variants.py to-vcf \\
        --table ${table} \\
        --fai ${fai} \\
        --output ${prefix}.vcf.gz \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.sites"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
