/*
 * GATK VariantFiltration - Filter SNPs
 */

process filter_snps {
    tag "filter_snps"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"

    input:
    tuple path(vcf), path(tbi)
    path reference
    path ref_index
    path ref_dict

    output:
    tuple path("cohort.snps.filtered.vcf.gz"), path("cohort.snps.filtered.vcf.gz.tbi"), emit: vcf

    script:
    def avail_mem = (task.memory.giga * 0.8).intValue()
    """
    gatk --java-options "-Xmx${avail_mem}g" VariantFiltration \\
        -R ${reference} \\
        -V ${vcf} \\
        --filter-expression "${params.snp_filter_expression}" \\
        --filter-name "SNP_FILTER" \\
        -O cohort.snps.filtered.vcf.gz
    """

    stub:
    """
    touch cohort.snps.filtered.vcf.gz
    touch cohort.snps.filtered.vcf.gz.tbi
    """
}
