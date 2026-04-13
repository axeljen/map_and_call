process samtools_faidx {
    tag "$reference"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path reference

    output:
    path "${reference}.fai", emit: fai

    script:
    """
    samtools faidx ${reference}
    """

    stub:
    """
    touch ${reference}.fai
    """
}
