process select_indels {
    scratch params.use_scratch
    tag "select_indels"
    label 'thin_short'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(region_id), path(vcf), path(csi)

    output:
    tuple val(region_id), path("region-${region_id}.indels.vcf.gz"), path("region-${region_id}.indels.vcf.gz.*"), emit: vcf

    script:
    """
    bcftools view -v indels -Oz -o region-${region_id}.indels.vcf.gz ${vcf}
    # Try to index; if it fails (e.g., no indels), create a valid empty VCF with header
    if ! bcftools index region-${region_id}.indels.vcf.gz 2>/dev/null; then
        echo "Warning: Could not index region-${region_id}.indels.vcf.gz (likely no indels found)"
        bcftools view -h ${vcf} | bgzip > region-${region_id}.indels.vcf.gz
        bcftools index region-${region_id}.indels.vcf.gz
    fi
    """

    stub:
    """
    touch region-${region_id}.indels.vcf.gz
    touch region-${region_id}.indels.vcf.gz.csi
    """
}
