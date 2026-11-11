// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process IGVREPORTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/igv-reports:1.16.0--pyh7e72e81_0':
        'quay.io/biocontainers/igv-reports:1.16.0--pyh7e72e81_0' }"

    input:
    // variantcomp patch: every track is staged in its own folder (tracks/<n>/), next to its
    // index (same position in tracks_indices), so tracks of the same name don't collide;
    // optional track configuration (--track-config, [] if none) and template ([] if none)
    tuple val(meta), path(sites), path(tracks, stageAs: 'tracks/??/*'), path(tracks_indices, stageAs: 'tracks/??/*'), path(track_config)
    tuple val(meta2), path(fasta), path(fai)
    path(template)

    output:
    tuple val(meta), path("*.html") , emit: report
    // variantcomp patch: igv_reports.__version__ is not maintained (0.1.0); the package version is
    tuple val("${task.process}"), val("igvreports"), eval("python -c \"import importlib.metadata; print(importlib.metadata.version('igv-reports'))\""), topic: versions, emit: versions_igvreports

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_opt = fasta ? "--fasta ${fasta}" : ""
    // variantcomp patch: the tracks are named in the track configuration, if one is given
    def config_opt = track_config ? "--track-config ${track_config}" : ""
    def template_opt = template ? "--template ${template}" : ""
    // If tracks is not null, create a string of the track paths
    // variantcomp patch: [ ].flatten(): a single staged track is a path, not a list
    def track_arg = tracks && !track_config ? "--tracks "+ [tracks].flatten().collect { track -> track.toString() }.join(' ') : ""
    // if "--tracks" is in the args, then add track_string immediately after it in
    // the args string and set the track_arg to ""
    if (args.contains("--tracks") && track_arg) {
        args = args.replace("--tracks", track_arg)
        track_arg = ""
    }

    """
    create_report ${sites} \\
        ${args} \\
        ${fasta_opt} \\
        ${track_arg} \\
        ${config_opt} \\
        ${template_opt} \\
        --output ${prefix}_report.html
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_report.html
    """
}
