process BENCHMARK_METRICS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    // vcfeval outputs, bcftools isec output folders (exact, position-only) and optional BEDs
    // that restrict the isec counts
    tuple val(meta), path(tp_baseline), path(fn), path(tp), path(fp), path(isec_exact, stageAs: 'isec_exact/*'), path(isec_position, stageAs: 'isec_position/*'), path(regions, stageAs: 'regions/*')

    output:
    tuple val(meta), path("${prefix}.metrics.tsv"), emit: metrics
    tuple val(meta), path("${prefix}.variants.tsv"), emit: variants
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions_arg = regions ? "--regions ${[regions].flatten().join(' ')}" : ''
    // ext.args: '--vaf-bins ...' (conf/modules.config)
    """
    # recall, precision and F1 per stratum and the per-variant table
    benchmark_metrics.py \\
        --tp-baseline ${tp_baseline} \\
        --fn ${fn} \\
        --tp ${tp} \\
        --fp ${fp} \\
        --isec-exact ${isec_exact} \\
        --isec-position ${isec_position} \\
        --sample ${meta.sample} \\
        --truth ${meta.truth} \\
        --query ${meta.query} \\
        --prefix ${prefix} \\
        ${regions_arg} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\ttruth\\tquery\\tstratum_type\\tstratum\\ttp_baseline\\tfn\\ttp_call\\tfp\\trecall\\tprecision\\tf1\\n' > ${prefix}.metrics.tsv
    touch ${prefix}.variants.tsv
    """
}
