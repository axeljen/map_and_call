/*
 * BWA index - Index reference genome for BWA
 */

process bwa_index {
    tag "$reference"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/01_reference_genome", mode: 'copy'

    input:
    path reference

    output:
    tuple path("${reference}"), path("${reference}.*"), emit: reference

    script:
    """
    bwa index ${reference}
    """

    stub:
    """
    touch ${reference}.amb
    touch ${reference}.ann
    touch ${reference}.bwt
    touch ${reference}.pac
    touch ${reference}.sa
    """
}
