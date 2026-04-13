/*
 * GATK CombineGVCFs - Combine per-sample GVCFs into a cohort GVCF
 */

process combine_gvcfs {
    tag "combine_gvcfs"
    label 'process_high'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(region_id), val(region), path(reference), path(gvcfs), path(gvcf_indexes), path(seqdict), path(reference_fai)

    output:
    tuple val(region_id), val(region), path("region_${region_id}.cohort.g.vcf.gz"), path("region_${region_id}.cohort.g.vcf.gz.tbi"), emit: gvcf

    script:
    def avail_mem = (task.memory.giga * 0.8).intValue()
    def gvcf_args = gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk --java-options "-Xmx${avail_mem}g" CombineGVCFs \
        -R ${reference} \
        ${gvcf_args} \
        -L ${region} \
        -O region_${region_id}.cohort.g.vcf.gz
    """

    stub:
    """
    touch region_${region_id}.cohort.g.vcf.gz
    touch region_${region_id}.cohort.g.vcf.gz.tbi
    """
}
