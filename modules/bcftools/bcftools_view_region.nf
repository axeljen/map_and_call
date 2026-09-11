/*
 * Subset a VCF file to a given region, to scatter a pre-existing VCF
 * across the same region chunks used elsewhere in the pipeline
 */
process bcftools_view_region {
    scratch params.use_scratch
    tag "region-${region_id}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(region_id), val(regions), path(vcf), path(idx)

    output:
    tuple val(region_id), path("region-${region_id}.vcf.gz"), path("region-${region_id}.vcf.gz.csi"), emit: vcf

    script:
    def regions_list = regions.join(',')
    """
    bcftools view -r ${regions_list} -Oz -o region-${region_id}.vcf.gz ${vcf}
    bcftools index -f region-${region_id}.vcf.gz
    """

    stub:
    """
    touch region-${region_id}.vcf.gz
    touch region-${region_id}.vcf.gz.csi
    """
}
