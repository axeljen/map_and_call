/*
 * FastQC - Quality control for raw and trimmed reads
 */

process fastqc {
    tag "$sample_id"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), val(lane), val(library), path(reads), val(stage)

    output:
    tuple val(sample_id), val(lane), val(library), path("${sample_id}*_${library}_${stage}/*.html"), emit: html
    tuple val(sample_id), val(lane), val(library), path("${sample_id}*_${library}_${stage}/*.zip"),  emit: zip

    script:
    def lanelabel = lane == 'combined' ? '' : "_${lane}"
    """
    mkdir -p ${sample_id}${lanelabel}_${library}_${stage}
    fastqc --threads ${task.cpus} --quiet ${reads} -o ${sample_id}${lanelabel}_${library}_${stage}
    """

    stub:
    """
    mkdir -p ${sample_id}${lanelabel}_${library}_${stage}
    touch ${sample_id}${lanelabel}_${library}_${stage}/fastqc_${sample_id}${lanelabel}_${library}_${stage}.html
    touch ${sample_id}${lanelabel}_${library}_${stage}/fastqc_${sample_id}${lanelabel}_${library}_${stage}.zip
    """
}