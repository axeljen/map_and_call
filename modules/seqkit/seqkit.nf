/*
 * seqkit split fastq files by readlength
 */

process split_fq_by_length {
    tag "$sample_id"
    conda "${moduleDir}/environment.yml"
    label 'medium_short'

    input:
    tuple val(sample_id), val(library), val(datatype), path(reads)
    

    output:
    tuple val(sample_id), val(library), val(datatype), path("${sample_id}_${library}_shortreads.fq.gz"), emit: shortreads
    tuple val(sample_id), val(library), val(datatype), path("${sample_id}_${library}_longreads.fq.gz"), emit: longreads

    script:
    def minlen_longreads = params.short_reads_threshold + 1
    """
    # extract short reads
    seqkit seq -M ${params.short_reads_threshold} -o ${sample_id}_${library}_shortreads.fq.gz ${reads}
    # extract long reads
    seqkit seq -m ${minlen_longreads} -o ${sample_id}_${library}_longreads.fq.gz ${reads}

    """

    stub:
    """
    touch ${sample_id}_${library}_longreads.fq.gz
    touch ${sample_id}_${library}_shortreads.fq.gz
    """
}
