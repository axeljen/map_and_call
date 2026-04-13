process bcftools_filter_indels {
    tag "filter_indels"    
    label 'process_wide'
    conda "${moduleDir}/environment.yml"


    input:
    tuple val(region_id), path(vcf), path(csi)

    output:
    tuple val(region_id), path("region-${region_id}.indels.filtered.vcf.gz"), path("region-${region_id}.indels.filtered.vcf.gz.*"), emit: vcf
    val(filter_expression)
    val(category)

    script:
    """
    bcftools filter -e ' ${params.indel_filter_expression} ' -Oz -o region-${region_id}.indels.filtered.vcf.gz ${vcf}
    bcftools index region-${region_id}.indels.filtered.vcf.gz
    """

    stub:
    """
    touch region-${region_id}.indels.filtered.vcf.gz
    touch region-${region_id}.indels.filtered.vcf.gz.csi
    """
}
