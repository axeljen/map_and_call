

process samtools_downsample {
    tag "samtools_downsample_${sample}"
    label 'medium_medium'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bamindex), val(fraction)
    
    output:
    tuple val(sample), path("${sample}_downsampled.bam"), path("${sample}_downsampled.bam.bai"), emit: downsampled_bam

    script:
    """
    samtools view -s ${fraction} -o ${sample}_downsampled.bam ${bam}
    samtools index ${sample}_downsampled.bam
    """

    stub:
    """
    touch ${sample}_downsampled.bam
    touch ${sample}_downsampled.bam.bai
    """
}
