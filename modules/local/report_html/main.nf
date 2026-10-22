process REPORT_HTML {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:e96f382ca64afe53'
        : 'community.wave.seqera.io/library/python_pandas_scipy_matplotlib-base_upsetplot:ddf4905651532171'}"

    input:
    // every table of the run; each is staged under a fixed name (the accounting and removed
    // tables of all VCFs in folders); optional tables are [] (power, callable, coverage,
    // known variants)
    tuple val(meta), path(per_pair, stageAs: 'in/per_pair.tsv'), path(metrics, stageAs: 'in/metrics.tsv'), path(power, stageAs: 'in/power.tsv'), path(variants, stageAs: 'in/variants.tsv'), path(truth_variants, stageAs: 'in/truth_variants.tsv'), path(classes, stageAs: 'in/classes.tsv'), path(enrichment, stageAs: 'in/enrichment.tsv'), path(genes, stageAs: 'in/genes.tsv'), path(annotation_counts, stageAs: 'in/annotation_counts.tsv')
    tuple path(accounting, stageAs: 'accounting/*'), path(removed, stageAs: 'removed/*'), path(region_steps, stageAs: 'in/region_steps.tsv'), path(targets, stageAs: 'in/targets.tsv'), path(callable, stageAs: 'in/callable.tsv'), path(coverage, stageAs: 'in/coverage.tsv')
    tuple path(required_summary, stageAs: 'in/required_summary.tsv'), path(required_variants, stageAs: 'in/required_variants.tsv'), path(known_benchmark, stageAs: 'in/known_benchmark.tsv')
    tuple path(checks, stageAs: 'in/checks.tsv'), path(samplesheet, stageAs: 'in/samplesheet.csv'), path(params_table, stageAs: 'in/params.tsv'), path(versions, stageAs: 'in/versions.yml')

    output:
    // versions are not sent to the `versions` topic: this process waits for the collated
    // versions (directly or through the reports), so it would keep the topic open forever
    tuple val(meta), path("${prefix}.html"), emit: html
    tuple val(meta), path("datavzrd/datavzrd.yaml"), path("datavzrd/*.tsv"), emit: datavzrd
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), emit: versions_python
    tuple val("${task.process}"), val('pandas'), eval("python -c 'import pandas; print(pandas.__version__)'"), emit: versions_pandas
    tuple val("${task.process}"), val('matplotlib'), eval("python -c 'import matplotlib; print(matplotlib.__version__)'"), emit: versions_matplotlib

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args (conf/modules/report.config): run information and links to the other reports
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "variantcomp_report"
    def power_arg = power ? "--power ${power}" : ''
    def callable_arg = callable ? "--callable ${callable}" : ''
    def coverage_arg = coverage ? "--coverage ${coverage}" : ''
    def required_arg = required_summary ? "--required-summary ${required_summary} --required-variants ${required_variants}" : ''
    def known_arg = known_benchmark ? "--known-benchmark ${known_benchmark}" : ''
    """
    # matplotlib needs a writable config/cache directory (the home may be read-only on HPC)
    export MPLCONFIGDIR=\$PWD/.matplotlib
    # narrative report (Summary, Introduction, Results, Methods) and the tables and
    # configuration of the datavzrd table browser
    report_html.py \\
        --per-pair ${per_pair} \\
        --metrics ${metrics} \\
        ${power_arg} \\
        --variants ${variants} \\
        --truth-variants ${truth_variants} \\
        --classes ${classes} \\
        --enrichment ${enrichment} \\
        --genes ${genes} \\
        --annotation-counts ${annotation_counts} \\
        --accounting ${accounting} \\
        --removed ${removed} \\
        --region-steps ${region_steps} \\
        --targets ${targets} \\
        ${callable_arg} \\
        ${coverage_arg} \\
        ${required_arg} \\
        ${known_arg} \\
        --checks ${checks} \\
        --samplesheet ${samplesheet} \\
        --params ${params_table} \\
        --versions ${versions} \\
        --html ${prefix}.html \\
        --datavzrd datavzrd \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "variantcomp_report"
    """
    mkdir datavzrd
    touch ${prefix}.html datavzrd/datavzrd.yaml datavzrd/truth-variants.tsv
    """
}
