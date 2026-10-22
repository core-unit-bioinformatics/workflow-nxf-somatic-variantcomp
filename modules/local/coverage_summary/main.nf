process COVERAGE_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    // alignment table (alignment, sample, technology, role) and the mosdepth outputs of
    // all alignments, named by the alignment id (common region, target BED, genome)
    tuple val(meta), path(alignments), path(mosdepth_files, stageAs: 'mosdepth/*')

    output:
    tuple val(meta), path("${prefix}.coverage.tsv"), emit: coverage
    path "${prefix}_generalstats_mqc.tsv", emit: multiqc
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args (conf/modules/coverage.config): --thresholds (= --coverage_thresholds)
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # mean depth and % bases >= each threshold per alignment, in the common region and in
    # its target region (genome without target BED)
    coverage_summary.py \\
        --alignments ${alignments} \\
        --dir mosdepth \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'alignment\\tsample\\ttechnology\\trole\\tregion\\tbp\\tmean_depth\\n' > ${prefix}.coverage.tsv
    touch ${prefix}_generalstats_mqc.tsv
    """
}
