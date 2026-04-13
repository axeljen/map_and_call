process samtools_markdups {
    tag "$sample_id"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    // publishDir "${params.outdir}/02_cramfiles", mode: 'copy'
    input:
    tuple val(sample_id), path(bam), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.cram"), path("${sample_id}.dedup.cram.crai"), emit: bam
    tuple val(sample_id), path("${sample_id}.dedup_metrics.txt"), emit: metrics

    script:
    """
    samtools collate -@ ${task.cpus} -O -u ${bam} | \
    samtools fixmate -@ ${task.cpus} -m -u - - | \
    samtools sort -@ ${task.cpus} -u - | \
    samtools markdup -@ ${task.cpus} -r -f ${sample_id}.dedup_metrics.txt --reference ${reference} -O CRAM --output-fmt-option version=3.0 - ${sample_id}.dedup.cram
    #--use-read-groups 
    samtools index -c ${sample_id}.dedup.cram

    """

    stub:
    """
    touch ${sample_id}.dedup.cram
    touch ${sample_id}.dedup.cram.crai
    touch ${sample_id}.dedup_metrics.txt
    """
}
