process IGV_TRACKS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:e96f382ca64afe53'
        : 'community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:ddf4905651532171'}"

    input:
    // one sample: a label <technology>:<tumor|normal|calls> per track file, and the track
    // files staged exactly as IGVREPORTS stages them (tracks/<n>/<name>), so that the
    // configuration names the files of the IGVREPORTS work directory; validated samplesheet
    tuple val(meta), val(labels), path(tracks, stageAs: 'tracks/??/*')
    path(samplesheet)

    output:
    // track configuration (igv.js track objects) and the SAM flags excluded when the reads
    // are embedded (create_report --exclude-flags), as a value
    tuple val(meta), path("${prefix}.tracks.json"), eval("cat ${prefix}.exclude_flags.txt"), emit: config
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args (conf/modules/igv_report.config): --height, --sampling-depth,
    // --amplicon-keep-duplicates
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # igv.js track objects of the tumor and normal alignments and the normalized calls of
    # every technology of the sample (truth technology first)
    igv_inputs.py tracks \\
        --sample ${meta.sample} \\
        --samplesheet ${samplesheet} \\
        --labels ${labels.join(' ')} \\
        --files ${[tracks].flatten().join(' ')} \\
        --output ${prefix}.tracks.json \\
        --flags-output ${prefix}.exclude_flags.txt \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '[]' > ${prefix}.tracks.json
    echo 1536 > ${prefix}.exclude_flags.txt
    """
}
