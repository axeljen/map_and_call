// Module for running some qc on the final processed bam/cram files

process cramqc {
    scratch params.use_scratch
    tag "$sample_id"
    label 'medium_medium'
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
