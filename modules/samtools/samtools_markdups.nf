process samtools_markdups {
    tag "$sample_id"
    conda "${moduleDir}/environment.yml"

    // publishDir "${params.outdir}/02_cramfiles", mode: 'copy'
    input:
    tuple val(sample_id), val(datatype), path(bam), path(bam_idx)

    output:
    tuple val(sample_id), val(datatype), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.dedup_metrics.txt"), emit: metrics

    script:
    """
    samtools collate -@ ${task.cpus} -O -u ${bam} | \
    samtools fixmate -@ ${task.cpus} -m -u - - | \
    samtools sort -@ ${task.cpus} -u - | \
    samtools markdup -@ ${task.cpus} -r -f ${sample_id}.dedup_metrics.txt -O BAM - ${sample_id}.dedup.bam
    #--use-read-groups 
    samtools index ${sample_id}.dedup.bam

    """

    stub:
    """
    touch ${sample_id}.dedup.bam
    touch ${sample_id}.dedup.bam.bai
    touch ${sample_id}.dedup_metrics.txt
    """
}
