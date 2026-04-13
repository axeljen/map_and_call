

process samtools_downsample {
    tag "samtools_downsample_${sample}"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bamindex), val(fraction)
    
    output:
    tuple val(sample), path("${sample}_downsampled.cram"), path("${sample}_downsampled.cram.crai"), emit: downsampled_bam

    script:
    """
    samtools view -s ${fraction} -C -o ${sample}_downsampled.cram ${bam}
    samtools index ${sample}_downsampled.cram
    """

    stub:
    """
    touch ${sample}_downsampled.cram
    touch ${sample}_downsampled.cram.crai
    """
}
