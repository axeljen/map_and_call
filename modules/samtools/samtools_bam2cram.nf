

process samtools_bam2cram {
    tag "samtools_bam2cram_${sample}"
    label 'medium_medium'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bamindex), path(reference_fasta), path(reference_fai), path(reference_gzi)
    
    output:
    tuple val(sample), path("${bam.baseName}.cram"), path("${bam.baseName}.cram.crai"), emit: cram

    script:
    """
    samtools view -C -T ${reference_fasta} -o ${bam.baseName}.cram ${bam}
    samtools index ${bam.baseName}.cram
    """

    stub:
    """
    touch ${bam.baseName}.cram
    touch ${bam.baseName}.cram.crai
    """
}
