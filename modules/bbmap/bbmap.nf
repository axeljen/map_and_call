/*
 * Clumpify - Read deduplication
 */

process clumpify_single {
    tag "$sample_id"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), val(lane), val(library), path(reads)
    
    output:
    tuple val(sample_id), val(lane), val(library), path("${sample_id}_${lane}_${library}.dedup.fastq.gz"), emit: dedup_reads

    script:
    """
    clumpify.sh \
        in=${reads} \
        out=${sample_id}_${lane}_${library}.dedup.fastq.gz \
        dedupe
        
    """

    stub:
    """
    touch ${sample_id}_${lane}_${library}.dedup.fastq.gz
    """
}

process clumpify_paired {
    tag "$sample_id"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), val(lane), val(library), path(reads1), path(reads2)
    
    output:
    tuple val(sample_id), val(lane), val(library), path("${sample_id}_${lane}_${library}.dedup_R1.fastq.gz"), path("${sample_id}_${lane}_${library}.dedup_R2.fastq.gz"), emit: dedup_reads

    script:
    """
    clumpify.sh \
        in=${reads1} \
        in2=${reads2} \
        out=${sample_id}_${lane}_${library}.dedup_R1.fastq.gz \
        out2=${sample_id}_${lane}_${library}.dedup_R2.fastq.gz \
        dedupe
        
    """

    stub:
    """
    touch ${sample_id}_${lane}_${library}.dedup_R1.fastq.gz
    touch ${sample_id}_${lane}_${library}.dedup_R2.fastq.gz
    """
}
