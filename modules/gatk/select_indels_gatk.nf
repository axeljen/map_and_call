/*
 * GATK SelectVariants - Select indels
 */

process select_indels_gatk {
    tag "select_indels"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"

    input:
    tuple path(vcf), path(tbi)
    path reference
    path ref_index
    path ref_dict

    output:
    tuple path("cohort.indels.vcf.gz"), path("cohort.indels.vcf.gz.tbi"), emit: vcf

    script:
    def avail_mem = (task.memory.giga * 0.8).intValue()
    """
    gatk --java-options "-Xmx${avail_mem}g" SelectVariants \\
        -R ${reference} \\
        -V ${vcf} \\
        --select-type-to-include INDEL \\
        -O cohort.indels.vcf.gz
    """

    stub:
    """
    touch cohort.indels.vcf.gz
    touch cohort.indels.vcf.gz.tbi
    """
}
