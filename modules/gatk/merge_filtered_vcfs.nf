/*
 * GATK MergeVcfs - Merge filtered VCFs
 */

process merge_filtered_vcfs {
    tag "merge_vcfs"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"

    // publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    tuple path(snps_vcf), path(snps_tbi)
    tuple path(indels_vcf), path(indels_tbi)
    path reference
    path ref_index
    path ref_dict

    output:
    tuple path("cohort.filtered.vcf.gz"), path("cohort.filtered.vcf.gz.tbi"), emit: vcf

    script:
    def avail_mem = (task.memory.giga * 0.8).intValue()
    """
    gatk --java-options "-Xmx${avail_mem}g" MergeVcfs \\
        -I ${snps_vcf} \\
        -I ${indels_vcf} \\
        -O cohort.filtered.vcf.gz
    """

    stub:
    """
    touch cohort.filtered.vcf.gz
    touch cohort.filtered.vcf.gz.tbi
    """
}
