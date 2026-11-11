process IGV_SITES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:e96f382ca64afe53'
        : 'community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:ddf4905651532171'}"

    input:
    // every compared variant (discordance/variants.annotated.tsv), the required known
    // variants ([] if none) and the reference index (contig order)
    tuple val(meta), path(variants, stageAs: 'in/variants.tsv'), path(required, stageAs: 'in/required_variants.tsv'), path(fai)

    output:
    // sites tables per sample and part: sites/<sample>/<sample>_<part>.sites.tsv (none if
    // there are no sites)
    tuple val(meta), path("sites/*/*.sites.tsv"), emit: sites, optional: true
    tuple val(meta), path("igv_links.tsv"), emit: links
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pandas'), eval("python -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args (conf/modules/igv_report.config): --max-sites, --flanking, --window
    def args = task.ext.args ?: ''
    def required_arg = required ? "--required ${required}" : ''
    """
    # one row per variant_id (compared and required variants), sorted by reference
    # position, split per sample into parts of at most --max-sites sites; links table
    # variant_id -> IGV report file
    igv_inputs.py sites \\
        --variants ${variants} \\
        ${required_arg} \\
        --fai ${fai} \\
        --outdir sites \\
        --links igv_links.tsv \\
        ${args}
    """

    stub:
    """
    mkdir -p sites/stub
    touch sites/stub/stub_001.sites.tsv
    printf 'variant_id\\tsample\\tpart\\tigv_file\\n' > igv_links.tsv
    """
}
