process CONSISTENCY_CHECKS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:e96f382ca64afe53'
        : 'community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:ddf4905651532171'}"

    input:
    // accounting and check tables of every VCF (variant_accounting.py), metrics, variant
    // table, discordance classes, region steps and summary; callable summary and required
    // known-variant summary are optional ([])
    tuple val(meta), path(accounting, stageAs: 'accounting/*'), path(per_pair, stageAs: 'in/per_pair.tsv'), path(aggregated, stageAs: 'in/aggregated.tsv'), path(variants, stageAs: 'in/variants.tsv'), path(classes, stageAs: 'in/classes.tsv'), path(region_steps, stageAs: 'in/region_steps.tsv'), path(targets, stageAs: 'in/targets.tsv'), path(callable, stageAs: 'in/callable.tsv'), path(required, stageAs: 'in/required.tsv')

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: checks
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "consistency_checks"
    def callable_arg = callable ? "--callable ${callable}" : ''
    def required_arg = required ? "--required ${required}" : ''
    """
    # do the counts of the steps add up? (accounting, metrics, variant table, classes,
    # regions, required variants); failed checks are listed, the run is stopped later
    # by CHECKS_GATE, after the report was written
    consistency_checks.py check \\
        --accounting ${accounting} \\
        --per-pair ${per_pair} \\
        --aggregated ${aggregated} \\
        --variants ${variants} \\
        --classes ${classes} \\
        --region-steps ${region_steps} \\
        --targets ${targets} \\
        ${callable_arg} \\
        ${required_arg} \\
        --output ${prefix}.tsv \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "consistency_checks"
    """
    printf 'category\\tscope\\tcheck\\texpected\\tobserved\\tstatus\\n' > ${prefix}.tsv
    """
}
