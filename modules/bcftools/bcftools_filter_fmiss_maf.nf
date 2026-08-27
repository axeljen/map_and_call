process bcftools_filter_fmiss_maf {
    scratch params.use_scratch
    tag "fmiss_maf_filter"
    label 'thin_medium'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(region_id), path(vcf), path(csi)
    val category

    output:
    tuple val(region_id), path("region-${region_id}.${category}.filtered.vcf.gz"), path("region-${region_id}.${category}.filtered.vcf.gz.*"), emit: vcf

    script:
    def filter_expression = "F_MISSING > ${params.fmiss_threshold} || MAF < ${params.maf_threshold}"
    // if min and max global dp are set in the params, add them to the filter expression
    if (params.min_global_dp != null) {
        filter_expression += " || DP < ${params.min_global_dp}"
    }
    if (params.max_global_dp != null) {
        filter_expression += " || DP > ${params.max_global_dp}"
    }
    """
    # add MAF and F_MISSING annotations to the vcf
    bcftools +fill-tags ${vcf} -Oz -- -t F_MISSING,MAF | \
        bcftools view -e "${filter_expression}" -Oz - | \
        # remove the GL field
        bcftools annotate -x 'FORMAT/GL' -Oz -o region-${region_id}.${category}.filtered.vcf.gz -
    # Try to index; if it fails (e.g., empty after filtering), create a valid empty VCF with header
    if ! bcftools index region-${region_id}.${category}.filtered.vcf.gz 2>/dev/null; then
        echo "Warning: Could not index region-${region_id}.${category}.filtered.vcf.gz (likely empty after filtering)"
        bcftools view -h ${vcf} | bcftools annotate -x 'FORMAT/GL' -Oz -o region-${region_id}.${category}.filtered.vcf.gz -
        bcftools index region-${region_id}.${category}.filtered.vcf.gz
    fi

    """
    stub:
    """
    touch region-${region_id}.${category}.filtered.vcf.gz
    touch region-${region_id}.${category}.filtered.vcf.gz.csi
    """
}
