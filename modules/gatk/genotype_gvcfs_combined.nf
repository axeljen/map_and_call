process genotype_gvcfs_combined {
    tag "combine_gvcfs"
    label 'process_high'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(region_id), val(region), path(reference), path(gvcf), path(gvcf_index), path(seqdict), path(reference_fai)

    output:
    tuple val(region_id), val(region), path("region_${region_id}.cohort.vcf.gz"), path("region_${region_id}.cohort.vcf.gz.tbi"), emit: vcf

    script:
    def avail_mem = (task.memory.giga * 0.8).intValue()
    """
    gatk --java-options "-Xmx${avail_mem}g" GenotypeGVCFs \
        -R ${reference} \
        -V ${gvcf} \
        -L ${region} \
        -O region_${region_id}.cohort.vcf.gz

    bcftools index region_${region_id}.cohort.vcf.gz
    """

    stub:
    """
    touch region_${region_id}.cohort.vcf.gz
    touch region_${region_id}.cohort.vcf.gz.tbi
    """
}
