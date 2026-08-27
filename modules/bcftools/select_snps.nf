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
    # Try to index; if it fails (e.g., no snps), create a valid empty VCF with header
    if ! bcftools index region-${region_id}.snps.vcf.gz 2>/dev/null; then
        echo "Warning: Could not index region-${region_id}.snps.vcf.gz (likely no snps found)"
        bcftools view -h ${vcf} | bgzip > region-${region_id}.snps.vcf.gz
        bcftools index region-${region_id}.snps.vcf.gz
    fi
    """

    stub:
    """
    touch region-${region_id}.snps.vcf.gz
    touch region-${region_id}.snps.vcf.gz.csi
    """
}
