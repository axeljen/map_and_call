/*
 * Fastp - Adapter trimming and quality filtering
 */

process fastp {
    tag "$sample_id"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), val(lane), val(datatype), val(library), path(read1), path(read2)

    output:
    tuple val(sample_id), val(lane), val(datatype), val(library), path("${sample_id}_${lane}_trimmed_R1.fastq.gz"), path("${sample_id}_${lane}_trimmed_R2.fastq.gz"), emit: reads
    tuple val(sample_id), val(lane), val(datatype), val(library), path("${sample_id}_${lane}.fastp.json"), emit: json
    tuple val(sample_id), val(lane), val(datatype), val(library), path("${sample_id}_${lane}.fastp.html"), emit: html

    script:
    """
    fastp \
        --in1 ${read1} \
        --in2 ${read2} \
        --out1 ${sample_id}_${lane}_trimmed_R1.fastq.gz \
        --out2 ${sample_id}_${lane}_trimmed_R2.fastq.gz \
        --json ${sample_id}_${lane}.fastp.json \
        --html ${sample_id}_${lane}.fastp.html \
        --trim_front1 ${params.trim_front} \
        --trim_front2 ${params.trim_front} \
        --trim_tail1 ${params.trim_tail} \
        --trim_tail2 ${params.trim_tail} \
        --length_required ${params.min_readlength} \
        --detect_adapter_for_pe \
        --thread ${task.cpus}

    """

    stub:
    """
    touch ${sample_id}_${lane}_trimmed_R1.fastq.gz
    touch ${sample_id}_${lane}_trimmed_R2.fastq.gz
    touch ${sample_id}_${lane}.fastp.json
    touch ${sample_id}_${lane}.fastp.html
    """
}
