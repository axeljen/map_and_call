process vcf_stats {
    tag "vcf_stats"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(region_id), path(vcf), path(csi)
    val category

    output:
    path("${category}_summary_statistics/${category}_${region_id}_ab_dp.txt"), emit: ab_dp
    path("${category}_summary_statistics/${category}_${region_id}_qual_fmiss_maf_dp.txt"), emit: qual_fmiss_maf
    path("${category}_summary_statistics/${category}_${region_id}_rec_counts.txt"), emit: rec_counts
    path("${category}_summary_statistics/${category}_${region_id}_sample_sumstats.txt"), emit: sample_stats

    script:
    """
    mkdir -p ${category}_summary_statistics
    vcfstats.py -i ${vcf} -o ${category}_summary_statistics/${category}_${region_id}    
    """
    stub:
    """
    mkdir -p ${category}_summary_statistics
    touch ${category}_summary_statistics/${category}_${region_id}_ab_dp.txt
    touch ${category}_summary_statistics/${category}_${region_id}_qual_fmiss_maf_dp.txt
    touch ${category}_summary_statistics/${category}_${region_id}_rec_counts.txt
    touch ${category}_summary_statistics/${category}_${region_id}_sample_sumstats.txt
    """
}
