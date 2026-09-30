/*
 * Compress (if needed) and index a user-supplied VCF file
 */
process bcftools_index_vcf {
    scratch params.use_scratch
    tag "index_vcf"
    conda "${moduleDir}/environment.yml"

    input:
    path vcf

    output:
    tuple path("input.vcf.gz"), path("input.vcf.gz.csi"), emit: indexed_vcf

    script:
    """
    bcftools view -Oz -o input.vcf.gz ${vcf}
    bcftools index -f input.vcf.gz
    """

    stub:
    """
    touch input.vcf.gz
    touch input.vcf.gz.csi
    """
}
