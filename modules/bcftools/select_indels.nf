process select_indels {
    tag "select_indels"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(region_id), path(vcf), path(csi)

    output:
    tuple val(region_id), path("region-${region_id}.indels.vcf.gz"), path("region-${region_id}.indels.vcf.gz.*"), emit: vcf

    script:
    """
    bcftools view -v indels -Oz -o region-${region_id}.indels.vcf.gz ${vcf}
    bcftools index region-${region_id}.indels.vcf.gz
    """

    stub:
    """
    touch region-${region_id}.indels.vcf.gz
    touch region-${region_id}.indels.vcf.gz.csi
    """
}
