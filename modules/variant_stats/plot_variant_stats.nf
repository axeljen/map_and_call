process plot_variant_stats {
    tag "plot_variant_stats_${category}"
    conda "${moduleDir}/environment.yml"

    input:
    path(stat_files)
    val category

    output:
    path "${category}_report.html", emit: report
    script:
    """
    # Create a temporary directory and copy stat files there
    mkdir -p stats_input
    cp ${stat_files} stats_input/

    # Run the R script to generate HTML report
    plot_variantstats.R stats_input

    # The R script creates a file based on the input directory name (stats_input_report.html)
    # Copy it to the expected output name
    OUTPUT_FILE=\$(find stats_input -name '*_report.html' -type f | head -1)
    if [ -z "\$OUTPUT_FILE" ]; then
        echo "Error: No report HTML file was generated"
        exit 1
    fi

    cp "\$OUTPUT_FILE" ${category}_report.html
    
    """

    stub:
    """
    touch ${category}_report.html
    """
}

