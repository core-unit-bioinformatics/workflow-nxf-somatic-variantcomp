process CHECK_REFERENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    // facts extracted by tools before this step:
    //   type 'vcf': [ header (bcftools view -h), record contigs (bcftools query -f '%CHROM\n') ]
    //   type 'bed': the BED itself
    //   type 'bam': samtools idxstats output
    tuple val(meta), path(facts, stageAs: 'facts/*')
    val type
    tuple val(meta2), path(fai)
    path regions // optional ('bam'): reference-named BED to translate to the alignment's contig names

    output:
    tuple val(meta), path("*.reference_check.tsv"), emit: report
    tuple val(meta), path("*.rename_chrs.tsv"), emit: rename
    tuple val(meta), path("*.keep_contigs.txt"), emit: keep
    tuple val(meta), path("*.harmonized.bed"), emit: bed, optional: true
    tuple val(meta), path("*.regions.bed"), emit: regions, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: e.g. '--samples TUMOR,NORMAL' for VCFs (set in conf/modules.config)
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${type}"
    def files = [facts].flatten()
    // the script option for each kind of extracted fact
    def facts_args = type == 'vcf'
        ? "--vcf-header ${files[0]} --record-contigs ${files[1]}"
        : type == 'bed' ? "--bed ${files[0]}" : "--idxstats ${files[0]}"
    def regions_arg = regions ? "--regions ${regions}" : ''
    """
    # decide on contig names/lengths (and sample names for VCFs); writes the rename map,
    # the contigs to keep and the report; exits 1 on a build mismatch or missing samples
    check_reference.py contigs \\
        --type ${type} \\
        --name ${meta.id} \\
        --fai ${fai} \\
        --prefix ${prefix} \\
        ${facts_args} \\
        ${regions_arg} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${type}"
    """
    printf 'key\\tvalue\\nstatus\\tok\\n' > ${prefix}.reference_check.tsv
    touch ${prefix}.rename_chrs.tsv ${prefix}.keep_contigs.txt
    if [ "${type}" == "bed" ]; then touch ${prefix}.harmonized.bed; fi
    if [ "${type}" == "bam" ] && [ -n "${regions}" ]; then touch ${prefix}.regions.bed; fi
    """
}
