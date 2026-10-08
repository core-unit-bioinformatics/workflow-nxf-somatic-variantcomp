process CHECK_REF_ALLELES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    tuple val(meta), path(norm_log) // stderr of bcftools norm --check-ref x (patched nf-core module)

    output:
    tuple val(meta), path("*.ref_alleles.tsv"), emit: report
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: '--max-ref-mismatch <fraction>' (conf/modules.config, --max_ref_mismatch)
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # REF mismatch rate from the bcftools norm summary; exits 1 above the threshold
    check_reference.py ref-alleles \\
        --name ${meta.id} \\
        --norm-log ${norm_log} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'key\\tvalue\\nstatus\\tok\\n' > ${prefix}.ref_alleles.tsv
    """
}
