process bcftools_filter {
    scratch params.use_scratch
    tag "filter_vcf"
    label 'thin_medium'
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
    # Try to index; if it fails (e.g., empty after filtering), create a valid empty VCF with header
    if ! bcftools index region-${region_id}.${category}.filtered.vcf.gz 2>/dev/null; then
        echo "Warning: Could not index region-${region_id}.${category}.filtered.vcf.gz (likely empty after filtering)"
        bcftools view -h ${vcf} | bgzip > region-${region_id}.${category}.filtered.vcf.gz
        bcftools index region-${region_id}.${category}.filtered.vcf.gz
    fi
    """

    stub:
    """
    touch region-${region_id}.${category}.filtered.vcf.gz
    touch region-${region_id}.${category}.filtered.vcf.gz.csi
    """
}
