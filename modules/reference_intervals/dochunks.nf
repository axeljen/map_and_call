// Simple process to split input files into chunks for parallel processing

process dochunks {
    tag "refintervals"
    label 'thin_short'
    conda "${moduleDir}/environment.yml"

    // publishDir "${params.outdir}/01_reference_genome", mode: 'copy'

    input:
    path refindex
    val chunk_size
    val scaffolds

    output:
    stdout emit: reference_intervals

    script:
    def scaffold_list = scaffolds.join(' ')
    """
    # make a temporary file listing the scaffolds to be chunked
    for scaffold in '${scaffold_list}';
        do 
        echo \$scaffold >> scaffolds_to_chunk.txt
    done

    # do the chunking
    dochunks.py $refindex ${chunk_size} scaffolds_to_chunk.txt
    """

    stub:
    """
    # create an empty file to emit something
    touch reference_intervals.txt
    """
}