process combine_summary_tables {
    scratch params.use_scratch
    tag "combine_summary_tables"
    label 'thin_short'
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