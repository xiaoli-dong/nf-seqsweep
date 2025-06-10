process KRAKENTOOLS_GETTOPMATCHES {
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"
    input:
    tuple val(meta), path(mpareport_count_perc)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    extrac_top_organism_from_mpa.py ${mpareport_count_perc} ${prefix}.top_matches.csv ${meta.id}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combine_kreports.py: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.top_matches.csv

    cat <<-END_VERSIONS > versions.yml
     "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
