include {
    FASTQC as fastqc_input ;
    FASTQC as fastqc_trimmed_with_bbduk ;
    FASTQC as fastqc_trimmed_with_fastp ;
    FASTQC as fastqc_dehost
} from '../../modules/nf-core/fastqc/main'

include { BBMAP_BBDUK } from '../../modules/nf-core/bbmap/bbduk/main'

include { FASTP } from '../../modules/nf-core/fastp/main.nf'

include {
    SEQKIT_STATS as stats_input ;
    SEQKIT_STATS as stats_trimmed_with_bbduk ;
    SEQKIT_STATS as stats_trimmed_with_fastp ;
    SEQKIT_STATS as stats_dehost
} from '../../modules/local/seqkit/stats/main.nf'

include {
    HTML_COPYDIR
} from '../../modules/local/html/copydir/main.nf'

include {
    HTML_DATATABLE_CSV2JSON as qcreport_to_json_for_datatable;
    HTML_DATATABLE_CSV2JSON as topmatches_to_json_for_datatable;
} from '../../modules/local/html/datatable/csv2json/main.nf'

include {
    HTML_JSTREE_KRAKEN2_MPA2JSON as kraken2mpa_to_json_for_jstree
} from '../../modules/local/html/jstree/kraken2/mpa2json/main.nf'


include {
    CSVTK_CONCAT as concat_stats_report;
    CSVTK_CONCAT as concat_topmatches;
} from '../../modules/nf-core/csvtk/concat/main.nf'

include {
    HOSTILE_CLEAN
} from '../../modules/local/hostile/clean/main.nf'

include {
    KRAKEN2_KRAKEN2
} from '../../modules/local/kraken2/kraken2/main.nf'

include {
    KRAKENTOOLS_COMBINEKREPORTS
} from '../../modules/nf-core/krakentools/combinekreports/main.nf'

include {
    KRAKENTOOLS_KREPORT2MPA
} from '../../modules/local/krakentools/kreport2mpa/main.nf'

include {
    KRAKENTOOLS_COMBINEMPA as KRAKENTOOLS_COMBINEMPA_COUNT ;
    KRAKENTOOLS_COMBINEMPA as KRAKENTOOLS_COMBINEMPA_PERC
} from '../../modules/local/krakentools/combinempa/main.nf'

include {
    KRAKENTOOLS_GETTOPMATCHES
} from '../../modules/local/krakentools/gettopmatches/main.nf'

include {
    CSVTK_JOIN as CSVTK_JOIN_SUMMARY_REPORT ;
    CSVTK_JOIN as CSVTK_JOIN_MPA_COUNT_PERC
} from '../../modules/nf-core/csvtk/join/main.nf'

include {
    REPORT_STATS
} from '../../modules/local/report/stats/main.nf'

workflow QC_ILLUMINA {
    take:
    ch_input_reads
    adapter_fasta
    host_ref_dir
    host_ref_name

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    trimmed_reads = ch_input_reads
    qc_reads = ch_input_reads
    ch_input_to_qc_summary = Channel.empty()

    ch_input_reads.filter { it[1][0].size() > 0 && it[1][0].countFastq() > 0 }
    fastqc_input(ch_input_reads)
    ch_versions = ch_versions.mix(fastqc_input.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(fastqc_input.out.zip)

    stats_input(ch_input_reads)
    ch_input_to_qc_summary = stats_input.out.stats
    ch_versions = ch_versions.mix(stats_input.out.versions.first())

    ch_multiqc_files.mix(stats_input.out.stats)
    qc_stats = Channel.empty()

    //default
    if (params.illumina_reads_qc_tool == 'bbduk') {
        BBMAP_BBDUK(ch_input_reads, adapter_fasta)
        ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
        //get rid of zero size contig file and avoid the downstream crash
        BBMAP_BBDUK.out.reads
            .filter { it[1][0].size() > 0 && it[1][0].countFastq() > 0 }
            .set { trimmed_reads }

        fastqc_trimmed_with_bbduk(trimmed_reads)
        stats_trimmed_with_bbduk(trimmed_reads)
        ch_multiqc_files = ch_multiqc_files.mix(fastqc_trimmed_with_bbduk.out.zip)

        qc_reads = trimmed_reads
        qc_stats = stats_trimmed_with_bbduk.out.stats
        ch_input_to_qc_summary = ch_input_to_qc_summary.join(qc_stats)
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBDUK.out.log)
    }
    else if (params.illumina_reads_qc_tool == 'fastp') {
        save_trimmed_fail = false
        discard_trimmed_pass = false
        save_merged = false
        FASTP(ch_input_reads, adapter_fasta, discard_trimmed_pass, save_trimmed_fail, save_merged)
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        //get rid of zero size contig file and avoid the downstream crash

        FASTP.out.reads
            .filter { it[1][0].size() > 0 && it[1][0].countFastq() > 0 }
            .set { trimmed_reads }

        //qc_reads = FASTP.out.reads
        fastqc_trimmed_with_fastp(trimmed_reads)
        stats_trimmed_with_fastp(trimmed_reads)
        ch_multiqc_files = ch_multiqc_files.mix(fastqc_trimmed_with_fastp.out.zip)

        qc_reads = trimmed_reads
        qc_stats = stats_trimmed_with_fastp.out.stats
        ch_input_to_qc_summary = ch_input_to_qc_summary.join(qc_stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json)
    }


    //FASTQC_QC(qc_reads)
    if (!params.skip_illumina_dehost) {
        HOSTILE_CLEAN(trimmed_reads, "auto", host_ref_dir, host_ref_name)
        ch_versions = ch_versions.mix(HOSTILE_CLEAN.out.versions.first())
        HOSTILE_CLEAN.out.fastq
            .filter { it[1][0].size() > 0 && it[1][0].countFastq() > 0 }
            .set { qc_reads }

        fastqc_dehost(qc_reads)
        stats_dehost(qc_reads)

        //qc_reads = HOSTILE_ILLUMINA.out.reads
        qc_stats = stats_dehost.out.stats
        ch_input_to_qc_summary = ch_input_to_qc_summary.join(qc_stats)
    }
    ch_input_to_qc_summary.view()

    KRAKEN2_KRAKEN2(
        qc_reads,
        Channel.value(file(params.kraken2_db)),
        false,
        true,
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)

    //produce two reports: count and percent
    KRAKENTOOLS_KREPORT2MPA(
        KRAKEN2_KRAKEN2.out.report
    )
    ch_versions = ch_versions.mix(KRAKENTOOLS_KREPORT2MPA.out.versions)
    ch_input_tsv = KRAKENTOOLS_KREPORT2MPA.out.count
        .join(KRAKENTOOLS_KREPORT2MPA.out.percent)
        .map { it ->
            [it[0], it[1..-1]]
        }
        .view()
    CSVTK_JOIN_MPA_COUNT_PERC(ch_input_tsv)

    KRAKENTOOLS_GETTOPMATCHES(CSVTK_JOIN_MPA_COUNT_PERC.out.csv)
    KRAKENTOOLS_GETTOPMATCHES.out.csv.view()
    concat_topmatches(
        KRAKENTOOLS_GETTOPMATCHES.out.csv.map { it ->
            it[1]
        }.collect().map { files ->
            tuple([id: "reads_illumina.topmatches"], files)
        },
        'csv',
        'csv',
    )
    KRAKENTOOLS_COMBINEKREPORTS(
        KRAKEN2_KRAKEN2.out.report.map { cfg, tsv ->
            tsv
        }.collect().map {
            it.sort()
        }.map { files ->
            [[id: "reads_illumina.kraken2_kreport"], files]
        }.view()
    )
    ch_versions = ch_versions.mix(KRAKENTOOLS_COMBINEKREPORTS.out.versions)

    KRAKENTOOLS_COMBINEMPA_COUNT(
        KRAKENTOOLS_KREPORT2MPA.out.count.map { cfg, tsv ->
            tsv
        }.collect().map { files ->
            [[id: "reads_illumina.kraken2_mpareport_count"], files.sort()]
        }
    )
    KRAKENTOOLS_COMBINEMPA_PERC(
        KRAKENTOOLS_KREPORT2MPA.out.percent.map { cfg, tsv ->
            tsv
        }.collect().map { files ->
            [[id: "reads_illumina.kraken2_mpareport_perc"], files.sort()]
        }
    )
    ch_versions = ch_versions.mix(KRAKENTOOLS_COMBINEMPA_COUNT.out.versions)

    ch_mpa_report = Channel.empty()
    ch_mpa_report = KRAKENTOOLS_COMBINEMPA_COUNT.out.txt
        .map { it ->
            it[1]
        }
        .combine(
            KRAKENTOOLS_COMBINEMPA_PERC.out.txt.map { it ->
                it[1]
            }
        )
        .map { it ->
            [[id: "reads_illumina.kraken2_mpareport_count_perc"], it]
        }
    CSVTK_JOIN_SUMMARY_REPORT(ch_mpa_report)
    kraken2mpa_to_json_for_jstree(CSVTK_JOIN_SUMMARY_REPORT.out.csv)
    ch_versions = ch_versions.mix(kraken2mpa_to_json_for_jstree.out.versions)
    REPORT_STATS(ch_input_to_qc_summary)
    ch_versions = ch_versions.mix(REPORT_STATS.out.versions)

    concat_stats_report(
        REPORT_STATS.out.csv.map { it ->
            it[1]
        }.collect().map { files ->
            tuple([id: "reads_illumina.qc_report"], files)
        },
        'csv',
        'csv',
    )
    ch_versions = ch_versions.mix(concat_stats_report.out.versions)

    ch_html_report_template = Channel.fromPath("${projectDir}/assets/html_report_template", checkIfExists: true)

    HTML_COPYDIR(ch_html_report_template, "html_report")

    qcreport_to_json_for_datatable(concat_stats_report.out.csv)
    ch_versions = ch_versions.mix(qcreport_to_json_for_datatable.out.versions)

    topmatches_to_json_for_datatable(concat_topmatches.out.csv)
    ch_versions = ch_versions.mix(qcreport_to_json_for_datatable.out.versions)


    emit:
    input_stats = stats_input.out.stats
    qc_reads
    qc_stats
    multiqc_files = ch_multiqc_files
    versions = ch_versions
}
