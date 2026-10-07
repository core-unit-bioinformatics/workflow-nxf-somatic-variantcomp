process ANNOTATE_VARIANTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0'
        : 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0'}"

    input:
    tuple val(meta), path(table)                     // TSV with columns chrom, pos (1-based), ref
    path gene_bed                                    // optional: BED4 (name = gene), reference contig names
    path strata, stageAs: 'strata/*'                 // optional: one BED per stratification, named <stratum>[.harmonized].bed
    path fai

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.annotated"
    if ("${table}" == "${prefix}.tsv") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    set -euo pipefail

    # variants as BED (REF span), row number as name
    awk -F'\\t' -v OFS='\\t' 'NR == 1 { for (i = 1; i <= NF; i++) col[\$i] = i; next }
        { print \$col["chrom"], \$col["pos"] - 1, \$col["pos"] - 1 + length(\$col["ref"]), NR - 1 }' ${table} \\
        | bedtools sort -g ${fai} -i stdin > variants.bed

    # gene names (comma-separated if several overlap, '.' if none)
    if [ -s "${gene_bed}" ]; then
        cut -f1-4 ${gene_bed} | bedtools sort -g ${fai} -i stdin > genes.sorted.bed
        bedtools map -a variants.bed -b genes.sorted.bed -g ${fai} -c 4 -o distinct \\
            | awk -F'\\t' -v OFS='\\t' '{ print \$4, \$5 }' > genes.tsv
    else
        touch genes.tsv
    fi

    # one yes/no column per stratification
    header=""
    files=""
    for bed in strata/*.bed; do
        [ -e "\$bed" ] || continue
        name=\$(basename "\$bed")
        name=\${name%.harmonized.bed}
        name=\${name%.bed}
        bedtools intersect -a variants.bed -b "\$bed" -c \\
            | awk -F'\\t' -v OFS='\\t' '{ print \$4, (\$5 > 0 ? "yes" : "no") }' > "stratum_\${name}.tsv"
        header="\${header}\\t\${name}"
        files="\${files} stratum_\${name}.tsv"
    done

    awk -F'\\t' -v OFS='\\t' -v has_genes=\$( [ -s genes.tsv ] && echo 1 || echo 0 ) -v extra="\$header" -v files="\$files" '
        BEGIN {
            n = split(files, f, " ")
            while ((getline line < "genes.tsv") > 0) { split(line, a, "\\t"); gene[a[1]] = a[2] }
            for (i = 1; i <= n; i++) {
                while ((getline line < f[i]) > 0) { split(line, a, "\\t"); flag[i, a[1]] = a[2] }
            }
        }
        NR == 1 { print \$0 (has_genes ? "\\tgene" : "") extra; next }
        {
            row = \$0
            if (has_genes) row = row "\\t" ((NR - 1) in gene ? gene[NR - 1] : ".")
            for (i = 1; i <= n; i++) row = row "\\t" flag[i, NR - 1]
            print row
        }' ${table} > ${prefix}.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.annotated"
    """
    cp ${table} ${prefix}.tsv
    """
}
