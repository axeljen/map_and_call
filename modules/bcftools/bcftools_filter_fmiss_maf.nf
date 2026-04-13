process bcftools_filter_fmiss_maf {
    tag "fmiss_maf_filter"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(region_id), path(vcf), path(csi)
    val category

    output:
    tuple val(region_id), path("region-${region_id}.${category}.filtered.vcf.gz"), path("region-${region_id}.${category}.filtered.vcf.gz.*"), emit: vcf

    script:
    def filter_expression = "F_MISSING >= ${params.fmiss_threshold} || MAF <= ${params.maf_threshold}"
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
        bcftools view -e "${filter_expression}" -Oz -o region-${region_id}.${category}.filtered.vcf.gz -
    bcftools index region-${region_id}.${category}.filtered.vcf.gz

    """
    stub:
    """
    touch region-${region_id}.${category}.filtered.vcf.gz
    touch region-${region_id}.${category}.filtered.vcf.gz.csi
    """
}
