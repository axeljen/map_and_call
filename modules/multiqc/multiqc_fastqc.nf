process multiqc_fastqc {
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    // publishDir "${params.outdir}/00_reports/00_fastqc_${stage}", mode: 'copy'

    input:
    tuple path(fastqc_zips), val(stage)

    output:
    path "${stage}_multiqc_report.html", emit: report
    path "${stage}_multiqc_report_data", emit: data

    script:
    """
    multiqc . --filename ${stage}_multiqc_report.html

    """

    stub:
    """
    mkdir -p ${stage}_multiqc_report_data
    touch ${stage}_multiqc_report.html
    """
}
