process VARIANT_ACCOUNTING {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    // input VCF (+ index, only staged), rename map, norm log and harmonization summary of
    // the compared branch, harmonized VCF (before the filters), compared VCF, region BED
    tuple val(meta), path(input_vcf, stageAs: 'input/*'), path(input_index, stageAs: 'input/*'), path(rename), path(norm_log), path(summary), path(harmonized), path(harmonized_tbi), path(compare), path(compare_tbi), path(regions)
    path genome_bed
    val pass_only // boolean: only PASS records were compared (the PREPROCESS_VCF input)

    output:
    tuple val(meta), path("${prefix}.accounting.tsv"), emit: accounting
    tuple val(meta), path("${prefix}.removed.tsv"), emit: removed
    tuple val(meta), path("${prefix}.checks.tsv"), emit: checks
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args (conf/modules/preprocess_vcf.config): the uniform filters
    def args = task.ext.args ?: ''
    def pass_arg = pass_only ? '--pass-only' : ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # records kept and removed per preprocessing step and reason, and checks that the
    # steps add up (the checks are evaluated at the end of the pipeline)
    variant_accounting.py \\
        --sample ${meta.sample} \\
        --technology ${meta.technology} \\
        --input-vcf ${input_vcf} \\
        --rename ${rename} \\
        --genome-bed ${genome_bed} \\
        --norm-log ${norm_log} \\
        --harmonize-summary ${summary} \\
        --harmonized ${harmonized} \\
        --compare ${compare} \\
        --regions ${regions} \\
        --prefix ${prefix} \\
        ${pass_arg} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\ttechnology\\tstep\\tstep_name\\toutcome\\treason\\trecords\\tsnv\\tindel\\tother\\n' > ${prefix}.accounting.tsv
    printf 'sample\\ttechnology\\tstep_name\\treason\\tchrom\\tpos\\tref\\talt\\tfilter\\n' > ${prefix}.removed.tsv
    printf 'sample\\ttechnology\\tcheck\\texpected\\tobserved\\tstatus\\n' > ${prefix}.checks.tsv
    """
}
