process READ_EVIDENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(alignment), path(index), path(sites, stageAs: 'sites/*') // VCFs whose positions are piled up
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''   // bcftools mpileup options, e.g. --ff for amplicon data
    prefix = task.ext.prefix ?: "${meta.id}.evidence"
    """
    # positions (and alleles) of all variants to look at
    for vcf in sites/*; do
        bcftools query -f '%CHROM\\t%POS\\n' "\$vcf"
    done | sort -u -k1,1 -k2,2n > sites.tsv

    # raw allele counts at these positions, split and normalized like the calls
    bcftools mpileup \\
        --fasta-ref ${fasta} \\
        --targets-file sites.tsv \\
        --annotate FORMAT/AD,FORMAT/DP \\
        --max-depth 100000 \\
        --min-MQ 20 \\
        --min-BQ 13 \\
        ${args} \\
        --output-type u \\
        ${alignment} \\
        | bcftools norm \\
            --fasta-ref ${fasta} \\
            --multiallelics -any \\
            --output-type z \\
            --write-index=tbi \\
            --output ${prefix}.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.evidence"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
