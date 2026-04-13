/*
 * GATK GenotypeGVCFs - Joint genotyping of cohort GVCF
 */

process genotype_gvcfs {
    tag "genotype_gvcfs"
    label 'process_high'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    tuple path(gvcf), path(tbi)
    path reference
    path ref_index
    path ref_dict

    output:
    tuple path("cohort.vcf.gz"), path("cohort.vcf.gz.tbi"), emit: vcf

    script:
    def avail_mem = (task.memory.giga * 0.8).intValue()
    """
    gatk --java-options "-Xmx${avail_mem}g" GenotypeGVCFs \\
        -R ${reference} \\
        -V ${gvcf} \\
        -O cohort.vcf.gz
    """
}
