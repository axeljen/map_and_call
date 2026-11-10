process samtools_stats {
    tag "$sample_id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    // publishDir "${params.outdir}/stats/mapping", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.stats.txt"), emit: stats
    tuple val(sample_id), path("${sample_id}.flagstat.txt"), emit: flagstat

    script:
    """
    // samtools stats ${bam} > ${sample_id}.stats.txt
    // samtools flagstat ${bam} > ${sample_id}.flagstat.txt
    echo "Simulated stats content for sample ${sample_id}" > ${sample_id}.stats.txt
    echo "Simulated flagstat content for sample ${sample_id}" > ${sample_id}.flagstat.txt
    """

    stub:
    """
    touch ${sample_id}.stats.txt
    touch ${sample_id}.flagstat.txt
    """
}
