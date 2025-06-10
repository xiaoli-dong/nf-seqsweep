process REPORT_STATS{
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'biocontainers/pandas:1.4.3' }"

    input:
    tuple val(meta), path(input_stats), path(trimmed_stats), path(dehost_stats)

    output:
    tuple val(meta), path("${prefix}_qc_report.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    prefix = task.ext.prefix ?: "${meta.id}"

    def input_tsv = input_stats ? "--input-stats  ${input_stats}" : ""
    def trimmed_tsv = trimmed_stats ? " --trim-stats ${trimmed_stats}" : ""
    def dehost_tsv = dehost_stats ? "--dehost-stats ${dehost_stats}" : ""

    """
    qc_summary.py ${input_tsv} ${trimmed_tsv} ${dehost_tsv} --output ${prefix}_qc_report.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
