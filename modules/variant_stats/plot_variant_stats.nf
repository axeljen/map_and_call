process plot_variant_stats {
    tag "plot_variant_stats_${category}"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/${params.variant_caller}/00_reports/02_variantstats", mode: 'copy'

    input:
    path(stat_files)
    val category

    output:
    path("*_report.pdf"), emit: variant_stats_pdf

    script:
    """
    # Create a temporary directory and copy stat files there
    mkdir -p stats_input
    cp ${stat_files} stats_input/
    
    # Run the R script to generate PDF report
    Rscript ${projectDir}/bin/plot_variantstats.R stats_input
    
    # The output PDF will be created as variant_stats_report.pdf or similar
    # Rename it with the category name for clarity  
    if [ -f "stats_input/variant_stats_report.pdf" ]; then
        cp stats_input/variant_stats_report.pdf ${category}_report.pdf
    elif [ -f "stats_input/"*"_report.pdf" ]; then
        cp stats_input/*_report.pdf ${category}_report.pdf
    else
        echo "Warning: No PDF report generated"
    fi
    """

    stub:
    """
    touch ${category}_report.pdf
    """
}

