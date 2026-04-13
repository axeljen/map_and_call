process bcftools_filter {
    tag "filter_vcf"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(region_id), path(vcf), path(csi)
    val filter_expression
    val category

    output:
    tuple val(region_id), path("region-${region_id}.${category}.filtered.vcf.gz"), path("region-${region_id}.${category}.filtered.vcf.gz.*"), emit: vcf

    script:
    """
    bcftools view -e ' ${filter_expression} ' -Oz -o region-${region_id}.${category}.filtered.vcf.gz ${vcf}
    bcftools index region-${region_id}.${category}.filtered.vcf.gz
    """

    stub:
    """
    touch region-${region_id}.${category}.filtered.vcf.gz
    touch region-${region_id}.${category}.filtered.vcf.gz.csi
    """
}
