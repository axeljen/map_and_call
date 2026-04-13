process multiqc_cram {
    tag "cram_multiqc"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    // publishDir "${params.outdir}/02_cramfiles", mode: 'copy'

    input:
    path(cram_qc_reports)

    output:
    path "multiqc_cram_report.html", emit: report
    path "multiqc_cram_report_data", emit: data

    script:
    """
    multiqc . --filename multiqc_cram_report.html
    
    """

    stub:
    """
    touch multiqc_cram_report.html
    mkdir -p multiqc_cram_report_data
    """
}
