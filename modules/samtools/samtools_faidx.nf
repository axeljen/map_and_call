process samtools_faidx {
    scratch params.use_scratch
    tag "$reference"
    label 'thin_short'
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
