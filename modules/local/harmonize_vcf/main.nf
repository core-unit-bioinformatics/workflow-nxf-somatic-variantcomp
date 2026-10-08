process HARMONIZE_VCF {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${prefix}.harmonize_summary.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.harmonized"
    def normal = meta.normal_id ? "--normal ${meta.normal_id}" : ''
    if ("${vcf}" == "${prefix}.vcf.gz") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    # add the VC_* INFO tags (VAF, depth, alt reads and their source, variant type),
    # set GT to 0/1 (tumor) and 0/0 (normal), drop symbolic alleles; the rules for
    # picking AD / AF / DP are documented in bin/harmonize_vcf.py and docs/usage.md
    harmonize_vcf.py \\
        --input ${vcf} \\
        --output ${prefix}.vcf.gz \\
        --tumor ${meta.tumor_id} \\
        ${normal} \\
        --index \\
        --summary ${prefix}.harmonize_summary.tsv \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.harmonized"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    printf 'key\\tvalue\\n' > ${prefix}.harmonize_summary.tsv
    """
}
