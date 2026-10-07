process COMMON_REGIONS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0'
        : 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0'}"

    input:
    tuple val(meta), path(target_beds, stageAs: 'targets/*') // sorted, contig-harmonized BEDs named <technology>.bed
    val genome_wide                                          // list of technologies without target BED
    path fai
    path regions_bed                                         // optional: restrict to these regions
    path exclude_bed                                         // optional: remove these regions

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val(meta), path("${prefix}.summary.tsv"), emit: summary
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def wgs_list = genome_wide.join(' ')
    """
    set -euo pipefail

    # genome-wide intervals of all reference contigs
    awk -v OFS='\\t' '{ print \$1, 0, \$2 }' ${fai} > genome.bed

    mkdir merged
    targets=()
    for bed in targets/*.bed; do
        [ -e "\$bed" ] || continue
        bedtools merge -i "\$bed" > "merged/\$(basename \$bed)"
        targets+=("merged/\$(basename \$bed)")
    done

    # N-way intersection of all target BEDs (WGS technologies are genome-wide)
    if [ \${#targets[@]} -eq 0 ]; then
        cp genome.bed common.raw.bed
    elif [ \${#targets[@]} -eq 1 ]; then
        cp "\${targets[0]}" common.raw.bed
    else
        bedtools multiinter -i "\${targets[@]}" \\
            | awk -v OFS='\\t' -v n=\${#targets[@]} '\$4 == n { print \$1, \$2, \$3 }' \\
            | bedtools merge -i stdin > common.raw.bed
    fi

    if [ -s "${regions_bed}" ] && [ "${regions_bed}" != "" ]; then
        bedtools intersect -a common.raw.bed -b ${regions_bed} \\
            | bedtools sort -g ${fai} -i stdin | bedtools merge -i stdin > common.tmp.bed
        mv common.tmp.bed common.raw.bed
    fi
    if [ -s "${exclude_bed}" ] && [ "${exclude_bed}" != "" ]; then
        bedtools subtract -a common.raw.bed -b ${exclude_bed} > common.tmp.bed
        mv common.tmp.bed common.raw.bed
    fi
    bedtools sort -g ${fai} -i common.raw.bed | bedtools merge -i stdin > ${prefix}.bed

    # summary: intervals and bp per technology and of the common region
    bp() { awk '{ s += \$3 - \$2 } END { printf "%d", s }' "\$1"; }
    common_bp=\$(bp ${prefix}.bed)
    {
        printf 'region\\ttype\\tintervals\\tbp\\tbp_in_common\\tfraction_in_common\\n'
        for bed in "\${targets[@]}"; do
            name=\$(basename "\$bed" .bed)
            total=\$(bp "\$bed")
            inside=\$(bedtools intersect -a "\$bed" -b ${prefix}.bed | awk '{ s += \$3 - \$2 } END { printf "%d", s }')
            awk -v n="\$name" -v t="\$total" -v i="\$inside" 'END { printf "%s\\ttarget\\t%d\\t%d\\t%d\\t%.4f\\n", n, NR, t, i, (t ? i / t : 0) }' "\$bed"
        done
        for name in ${wgs_list}; do
            awk -v n="\$name" -v c="\$common_bp" '{ t += \$3 - \$2 } END { printf "%s\\tgenome\\t%d\\t%d\\t%d\\t%.4f\\n", n, NR, t, c, (t ? c / t : 0) }' genome.bed
        done
        awk -v OFS='\\t' '{ t += \$3 - \$2 } END { printf "common\\tcommon\\t%d\\t%d\\t%d\\t1.0000\\n", NR, t, t }' ${prefix}.bed
    } > ${prefix}.summary.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    printf 'region\\ttype\\tintervals\\tbp\\tbp_in_common\\tfraction_in_common\\n' > ${prefix}.summary.tsv
    """
}
