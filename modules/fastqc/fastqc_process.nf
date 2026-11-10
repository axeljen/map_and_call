/*
 * FastQC - Quality control for raw and trimmed reads
 */

process fastqc {
    tag "$sample_id"
    label 'medium_short'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), val(lane), val(library), path(reads), val(stage)

    output:
    tuple val(sample_id), val(lane), val(library), path("${sample_id}/*.html"), emit: html
    tuple val(sample_id), val(lane), val(library), path("${sample_id}/*.zip"),  emit: zip

    script:
    def lanelabel = lane == 'combined' ? '' : "_${lane}"
    """
    mkdir -p ${sample_id}
    fastqc --threads ${task.cpus} --quiet ${reads} -o ${sample_id}
    """

    stub:
    """
    mkdir -p ${sample_id}_${reads.baseName}
    touch ${sample_id}_${reads.baseName}/fastqc_${sample_id}_${reads.baseName}.html
    touch ${sample_id}_${reads.baseName}/fastqc_${sample_id}_${reads.baseName}.zip
    """
}