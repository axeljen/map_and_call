
process bcftools_concat {
    scratch params.use_scratch
    tag "concat_vcfs"
    label 'thin_medium'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/03_genotypes", mode: 'copy'

    input:
    tuple path(vcf_files), path(vcf_indices)
    val category

    output:
    tuple path("${category}.vcf.gz"), path("${category}.vcf.gz.*"), emit: vcf

    script:
    """
    bcftools concat -Oz -o ${category}.vcf.gz ${vcf_files.join(' ')}
    bcftools index ${category}.vcf.gz
    """

    stub:
    """
    touch ${category}.vcf.gz
    touch ${category}.vcf.gz.csi
    """
}