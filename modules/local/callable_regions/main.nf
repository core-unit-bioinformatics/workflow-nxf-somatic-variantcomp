process CALLABLE_REGIONS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0'
        : 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0'}"

    input:
    tuple val(meta), path(coverage, stageAs: 'coverage/*'), path(renames, stageAs: 'rename/*') // mosdepth <id>.{quantized,regions}.bed.gz and <id>.bam.rename_chrs.tsv per alignment
    path common_bed
    path fai
    val mode         // 'base': mosdepth quantized CALLABLE bins; 'region': mean depth per region
    val min_coverage

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val(meta), path("${prefix}.summary.tsv"), emit: summary
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.callable"
    """
    set -euo pipefail
    cut -f1 ${fai} > reference_contigs.txt
    mkdir per_alignment
    files=()
    for cov in coverage/*.bed.gz; do
        id=\$(basename "\$cov")
        id=\${id%.quantized.bed.gz}
        id=\${id%.regions.bed.gz}
        if [ "${mode}" == "base" ]; then
            zcat "\$cov" | awk -v OFS='\\t' '\$4 == "CALLABLE" { print \$1, \$2, \$3 }'
        else
            zcat "\$cov" | awk -v OFS='\\t' -v min=${min_coverage} '\$4 >= min { print \$1, \$2, \$3 }'
        fi \\
            | awk -v OFS='\\t' 'FILENAME == ARGV[1] { map[\$1] = \$2; next }
                                FILENAME == ARGV[2] { ref[\$1] = 1; next }
                                { if (\$1 in map) \$1 = map[\$1]; if (\$1 in ref) print }' \\
                "rename/\${id}.bam.rename_chrs.tsv" reference_contigs.txt - \\
            | bedtools sort -g ${fai} -i stdin \\
            | bedtools merge -i stdin > "per_alignment/\${id}.bed"
        files+=("per_alignment/\${id}.bed")
    done

    # callable in every alignment of the sample
    if [ \${#files[@]} -eq 1 ]; then
        cp "\${files[0]}" all.bed
    else
        bedtools multiinter -i "\${files[@]}" \\
            | awk -v OFS='\\t' -v n=\${#files[@]} '\$4 == n { print \$1, \$2, \$3 }' \\
            | bedtools merge -i stdin > all.bed
    fi
    bedtools intersect -a all.bed -b ${common_bed} \\
        | bedtools sort -g ${fai} -i stdin \\
        | bedtools merge -i stdin > ${prefix}.bed

    bp() { awk '{ s += \$3 - \$2 } END { printf "%d", s }' "\$1"; }
    {
        printf 'sample\\talignments\\tmode\\tmin_coverage\\tcommon_bp\\tcallable_bp\\tfraction_callable\\n'
        awk -v OFS='\\t' -v s="${meta.id}" -v n=\${#files[@]} -v m="${mode}" -v c=${min_coverage} \\
            -v total=\$(bp ${common_bed}) -v callable=\$(bp ${prefix}.bed) \\
            'BEGIN { printf "%s\\t%d\\t%s\\t%s\\t%d\\t%d\\t%.4f\\n", s, n, m, c, total, callable, (total ? callable / total : 0) }'
    } > ${prefix}.summary.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.callable"
    """
    touch ${prefix}.bed
    printf 'sample\\talignments\\tmode\\tmin_coverage\\tcommon_bp\\tcallable_bp\\tfraction_callable\\n' > ${prefix}.summary.tsv
    """
}
