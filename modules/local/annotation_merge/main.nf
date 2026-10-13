process ANNOTATION_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    tuple val(meta), path(table), path(vep), path(repeats), path(strata, stageAs: 'strata/*')
    // vep: bcftools +split-vep TSV (or []); repeats: bedtools map result with the repeat BED
    // (or []); strata: bedtools intersect -c results named <stratum>.counts.bed

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.annotated"
    def vep_arg = vep ? "--vep ${vep}" : ''
    def repeats_arg = repeats ? "--repeats ${repeats}" : ''
    // one NAME=FILE argument per stratification result; the name is the base name of the
    // staged file ('strata/<stratum>.counts.bed') without directory and suffix
    def strata_args = [strata].flatten().collect { file ->
        "${file.toString().tokenize('/').last().replaceAll(/\.counts\.bed$/, '')}=${file}"
    }
    """
    # add the VEP consequence, gene, impact, variant class, the repeat context and one
    # yes/no column per stratification to the variant table
    annotate_variants.py merge \\
        --table ${table} \\
        ${vep_arg} \\
        ${repeats_arg} \\
        ${strata_args ? '--stratum ' + strata_args.join(' ') : ''} \\
        --output ${prefix}.tsv \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.annotated"
    """
    cp ${table} ${prefix}.tsv
    """
}
