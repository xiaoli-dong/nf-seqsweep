process KRAKENTOOLS_KREPORT2MPA {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakentools:1.2--pyh5e36f6f_0':
        'biocontainers/krakentools:1.2--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(kreport)

    output:
    tuple val(meta), path("*_count.txt"), emit: count
    tuple val(meta), path("*_perc.txt"), emit: percent
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    kreport2mpa.py ${args} --read_count -r ${kreport} -o ${prefix}.kraken2_mpareport_count.txt
    kreport2mpa.py ${args} --percentages -r ${kreport} -o ${prefix}.kraken2_mpareport_perc.txt

    awk '
        NR==1 {
        printf "#Classification"
        for (i=2; i<=NF; i++) {
            gsub(".kraken2_kreport.txt", "", \$i)
            printf "\\t%s_count", \$i
        }
        print ""
        next
        }
        { print }
    ' ${prefix}.kraken2_mpareport_count.txt > count.txt

    mv count.txt ${prefix}.kraken2_mpareport_count.txt

     awk '
        NR==1 {
        printf "#Classification"
        for (i=2; i<=NF; i++) {
            gsub(".kraken2_kreport.txt", "", \$i)
            printf "\\t%s_perc", \$i
        }
        print ""
        next
        }
        { print }
    ' ${prefix}.kraken2_mpareport_perc.txt > perc.txt

    mv perc.txt ${prefix}.kraken2_mpareport_perc.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combine_kreports.py: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kreport2mpa.py: ${VERSION}
    END_VERSIONS
    """
}
