process CHECK_REFERENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    // the input to check: a VCF/BCF, an alignment (BAM/CRAM) or a BED; its index is
    // optional ([]) and only staged next to it
    tuple val(meta), path(input), path(input_index)
    val type                                  // 'vcf', 'bam' or 'bed'
    tuple val(meta2), path(fasta), path(fai)  // reference (FASTA: REF alleles of VCFs, CRAM)

    output:
    tuple val(meta), path("*.reference_check.tsv"), emit: report
    tuple val(meta), path("*.rename_chrs.tsv"), emit: rename
    tuple val(meta), path("*.harmonized.bed"), emit: bed, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args (conf/modules/input_checks.config): --samples, --max-ref-mismatch,
    // --skip-ref-alleles, --keep-mito, --warn-only
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${type}"
    // the option that names the input for each type
    def input_arg = type == 'vcf' ? "--vcf ${input} --fasta ${fasta}"
        : type == 'bam' ? "--alignment ${input} --fasta ${fasta}"
        : "--bed ${input}"
    """
    # check contig names and lengths (and samples and REF alleles of VCFs) against the
    # reference; writes the rename map, the report and (BED) the harmonized BED;
    # exits 1 on an error
    check_reference.py ${type} \\
        --name ${meta.id} \\
        --fai ${fai} \\
        --prefix ${prefix} \\
        ${input_arg} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${type}"
    """
    printf 'key\\tvalue\\nstatus\\tok\\n' > ${prefix}.reference_check.tsv
    touch ${prefix}.rename_chrs.tsv
    if [ "${type}" == "bed" ]; then touch ${prefix}.harmonized.bed; fi
    """
}
