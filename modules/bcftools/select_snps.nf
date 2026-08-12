process select_snps {
    scratch params.use_scratch
    tag "select_snps"
    label 'thin_short'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(region_id), path(vcf), path(csi)

    output:
    tuple val(region_id), path("region-${region_id}.snps.vcf.gz"), path("region-${region_id}.snps.vcf.gz.*"), emit: vcf

    script:
    """
    bcftools view -v snps -Oz -o region-${region_id}.snps.vcf.gz ${vcf}
    bcftools index region-${region_id}.snps.vcf.gz
    """

    stub:
    """
    touch region-${region_id}.snps.vcf.gz
    touch region-${region_id}.snps.vcf.gz.csi
    """
}
