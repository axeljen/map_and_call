process combine_stats {
    tag "combine_stats"
    conda "${moduleDir}/environment.yml"

    // publishDir "${params.outdir}/04_variant_stats", mode: 'copy'

    input:
    path ab_dp_files
    path gt_files 
    path sample_stats_files 
    path rec_count_files 
    val category

    output:
    path("${category}*"), emit: combined_summary_statistics

    script:
    """
    #mkdir -p ${category}_summary_statistics
    combine_stats.py --ab-dp ${ab_dp_files.join(',')} --gt ${gt_files.join(',')} --sample-stats ${sample_stats_files.join(',')} --rec-counts ${rec_count_files.join(',')} -o ${category}
    """
    stub:
    """
    #mkdir -p ${category}_summary_statistics
    touch ${category}_qual_fmiss_maf_dp.tsv
    touch ${category}_sample_sumstats.tsv
    touch ${category}_record_counts.tsv
    touch ${category}_ab_dp.tsv
    """
}
