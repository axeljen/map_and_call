// Module for running some qc on the final processed bam/cram files

process cramqc {
    tag "$sample_id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    
    input:
    tuple val(sample_id), path(cram), path(crai)

    output:
    tuple val(sample_id), path("${sample_id}_mapqc.txt"), emit: qc_report

    script:
    """
    samtools stats ${cram} > ${sample_id}_mapqc.txt
    
    """

    stub:
    """
    touch ${sample_id}_mapqc.txt
    """
}
