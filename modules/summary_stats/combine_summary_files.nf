process combine_summary_tables {
    tag "combine_summary_tables"
    label 'process_small'
    // conda "${moduleDir}/environment.yml"

    input:
        path(summary_files)

    output:
        path("summary_statistics.tsv"), emit: table

    script:
    """
    awk 'FNR==1 && NR!=1{next}1' ${summary_files} > summary_statistics.tsv
    """

    stub:
    """
    touch summary_statistics.tsv
    """
}