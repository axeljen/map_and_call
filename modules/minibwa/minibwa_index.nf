/*
 * minibwa index - index reference genome for downstream minibwa mapping
 */

process minibwa_index {
    scratch params.use_scratch
    tag "$reference"
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/01_reference_genome", mode: 'copy'

    input:
    path reference

    output:
    tuple path("${reference}"), path("${reference}.*"), emit: reference

    script:
    """
    minibwa index ${reference}
    """

    stub:
    """
    touch ${reference}.l2b
    touch ${reference}.mbw
    """
}