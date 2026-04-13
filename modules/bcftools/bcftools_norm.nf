
process bcftools_norm {
    tag "normalize_vcfs"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(region_id), path(vcf), path(csi), path(reference), path(reference_fai), path(reference_gzi)

    output:
    tuple val(region_id), path("region-${region_id}.norm.vcf.gz"), path("region-${region_id}.norm.vcf.gz.*"), emit: vcf

    script:
    """
    bcftools norm --check-ref w -f ${reference} -a --atom-overlaps '*' -m -both -Oz -o region-${region_id}.norm.vcf.gz ${vcf}
    bcftools index region-${region_id}.norm.vcf.gz
    """

    stub:
    """
    touch region-${region_id}.norm.vcf.gz
    touch region-${region_id}.norm.vcf.gz.csi
    """
}