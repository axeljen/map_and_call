/*
 * GATK SelectVariants - Select SNPs
 */

process select_snps_gatk {
    tag "select_snps"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"

    input:
    tuple path(vcf), path(tbi)
    path reference
    path ref_index
    path ref_dict

    output:
    tuple path("cohort.snps.vcf.gz"), path("cohort.snps.vcf.gz.tbi"), emit: vcf

    script:
    def avail_mem = (task.memory.giga * 0.8).intValue()
    """
    gatk --java-options "-Xmx${avail_mem}g" SelectVariants \\
        -R ${reference} \\
        -V ${vcf} \\
        --select-type-to-include SNP \\
        -O cohort.snps.vcf.gz
    """

    stub:
    """
    touch cohort.snps.vcf.gz
    touch cohort.snps.vcf.gz.tbi
    """
}
