process KRAKEN2_KRAKEN2 {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kraken2%3A2.1.5--pl5321h077b44d_0' :
        'biocontainers/kraken2:2.1.5--pl5321h077b44d_0' }"

    input:
    tuple val(meta), path(reads)
    path  db
    val save_output_fastqs
    val save_reads_assignment

    output:
    tuple val(meta), path('*.classified{.,_}*')     , optional:true, emit: classified_reads_fastq
    tuple val(meta), path('*.unclassified{.,_}*')   , optional:true, emit: unclassified_reads_fastq
    tuple val(meta), path('*classifiedreads.txt')   , optional:true, emit: classified_reads_assignment
    tuple val(meta), path('*report.txt')                           , emit: report
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired       = meta.single_end ? "" : "--paired"
    def classified   = meta.single_end ? "${prefix}.classified.fastq"   : "${prefix}.classified#.fastq"
    def unclassified = meta.single_end ? "${prefix}.unclassified.fastq" : "${prefix}.unclassified#.fastq"
    def classified_option = save_output_fastqs ? "--classified-out ${classified}" : ""
    def unclassified_option = save_output_fastqs ? "--unclassified-out ${unclassified}" : ""
    def readclassification_option = save_reads_assignment ? "--output ${prefix}.kraken2.classifiedreads.txt" : "--output /dev/null"
    //def compress_reads_command = save_output_fastqs ? "pigz -p $task.cpus *.fastq" : ""

    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${prefix}.kraken2_kreport.txt \\
        --gzip-compressed \\
        $unclassified_option \\
        $classified_option \\
        $readclassification_option \\
        $paired \\
        $args \\
        $reads



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')

    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired       = meta.single_end ? "" : "--paired"
    def classified   = meta.single_end ? "${prefix}.classified.fastq.gz"   : "${prefix}.classified_1.fastq.gz ${prefix}.classified_2.fastq.gz"
    def unclassified = meta.single_end ? "${prefix}.unclassified.fastq.gz" : "${prefix}.unclassified_1.fastq.gz ${prefix}.unclassified_2.fastq.gz"
    def readclassification_option = save_reads_assignment ? "--output ${prefix}.kraken2.classifiedreads.txt" : "--output /dev/null"


    """
    touch ${prefix}.kraken2.report.txt
    if [ "$save_output_fastqs" == "true" ]; then
        touch $classified
        touch $unclassified
    fi
    if [ "$save_reads_assignment" == "true" ]; then
        touch ${prefix}.kraken2.classifiedreads.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')

    END_VERSIONS
    """

}
