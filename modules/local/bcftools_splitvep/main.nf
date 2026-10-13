process BCFTOOLS_SPLITVEP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(vcf), path(tbi) // VCF annotated by VEP (INFO/CSQ)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: the fields and the transcript selection, e.g.
    // -H -f '%CHROM\t%POS\t...\n' -s worst (conf/modules/annotation.config)
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.vep"
    """
    # one line per variant with the VEP fields of the selected consequence
    bcftools +split-vep \\
        ${args} \\
        ${vcf} \\
        > ${prefix}.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.vep"
    """
    touch ${prefix}.tsv
    """
}
