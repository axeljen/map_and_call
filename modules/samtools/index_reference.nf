
process samtools_index {
    tag "$reference"
    label 'thin_short'
    conda "${moduleDir}/environment.yml"

    // publishDir "${params.outdir}/01_reference_genome", mode: 'copy'

    input:
    path reference

    output:
    path "${reference}", emit: reference_fasta
    path "${reference}.fai", emit: reference_fai
    path "${reference}.gzi", emit: reference_gzi

    script:
    """
    samtools faidx ${reference}
    # .gzi is only created for bgzipped FASTA; create a placeholder for plain FASTA
    if [ ! -f ${reference}.gzi ]; then
        touch ${reference}.gzi
    fi
    
    """

    stub:
    """
    touch ${reference}.fai
    touch ${reference}.gzi
    echo "chr1:1-10000"
    """
}

