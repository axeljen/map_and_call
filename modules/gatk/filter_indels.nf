/*
 * GATK VariantFiltration - Filter indels
 */

process filter_indels {
    tag "filter_indels"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"

    input:
    tuple path(vcf), path(tbi)
    path reference
    path ref_index
    path ref_dict

    output:
    tuple path("cohort.indels.filtered.vcf.gz"), path("cohort.indels.filtered.vcf.gz.tbi"), emit: vcf

    script:
    def avail_mem = (task.memory.giga * 0.8).intValue()
    """
    gatk --java-options "-Xmx${avail_mem}g" VariantFiltration \\
        -R ${reference} \\
        -V ${vcf} \\
        --filter-expression "${params.indel_filter_expression}" \\
        --filter-name "INDEL_FILTER" \\
        -O cohort.indels.filtered.vcf.gz
    """

    stub:
    """
    touch cohort.indels.filtered.vcf.gz
    touch cohort.indels.filtered.vcf.gz.tbi
    """
}
