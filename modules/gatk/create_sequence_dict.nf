/*
 * GATK CreateSequenceDictionary - Create sequence dictionary for reference
 */

process create_sequence_dict {
    tag "$reference"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    // publishDir "${params.outdir}/00_reference", mode: 'copy'

    input:
    path reference

    output:
    path "${reference.baseName}.dict", emit: dict

    script:
    """
    gatk CreateSequenceDictionary \\
        -R ${reference} \\
        -O ${reference.baseName}.dict
    """

    stub:
    """
    touch ${reference.baseName}.dict
    """
}
