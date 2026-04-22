// Process to concatenate readpairs stemming from the same library, prior to deduplication.

process concat_reads {
    tag "$sample_id"

    input:
    tuple val(sample_id), val(library), val(datatype), path(forward_reads), path(reverse_reads)
    
    output:
    tuple val(sample_id), val(library), val(datatype), path("${sample_id}_${library}_concat_R1.fastq.gz"), path("${sample_id}_${library}_concat_R2.fastq.gz"), emit: reads_concat

    script:
    """
    cat ${forward_reads} > ${sample_id}_${library}_concat_R1.fastq.gz
    cat ${reverse_reads} > ${sample_id}_${library}_concat_R2.fastq.gz
    """

    stub:
    """
    touch ${sample_id}_${library}_concat_R1.fastq.gz
    touch ${sample_id}_${library}_concat_R2.fastq.gz
    """
}

process concat_collapsed {
    tag "$sample_id"

    input:
    tuple val(sample_id), val(library), val(datatype), path(collapsed_reads)
    
    output:
    tuple val(sample_id), val(library), val(datatype), path("${sample_id}_${library}_concat_collapsed.fastq.gz"), emit: collapsed_concat

    script:
    """
    cat ${collapsed_reads} > ${sample_id}_${library}_concat_collapsed.fastq.gz
    """

    stub:
    """
    touch ${sample_id}_${library}_concat_collapsed.fastq.gz
    """
}