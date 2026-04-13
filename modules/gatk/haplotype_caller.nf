/*
 * GATK HaplotypeCaller - Standard variant calling mode
 */

process haplotype_caller {
    tag "$sample_id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), path(cram), path(crai), path(reference), path(reference_indices), val(region_id), val(regions)

    output:
    tuple val(region_id), val(sample_id), path("${sample_id}.region-${region_id}.vcf.gz*"), emit: vcf

    script:
    def avail_mem = (task.memory.giga * 0.8).intValue()
    """
    # put the regions in an intervals file that we'll pass to HaplotypeCaller
    for region in \$(echo "${regions}" | tr ',' ' '); do
        echo \$region >> ${sample_id}.region-${region_id}.intervals
    done

    gatk --java-options "-Xmx${avail_mem}g" HaplotypeCaller \
        -R ${reference} \
        -I ${cram} \
        -L ${sample_id}.region-${region_id}.intervals \
        -ip 100 \
        -O ${sample_id}.region-${region_id}.vcf.gz \
        --sample-ploidy ${params.ploidy}

    # get rid of the padding to avoid overlaps on concatenation
    bcftools view -r ${regions} -Oz -o tmp.vcf.gz ${sample_id}.region-${region_id}.vcf.gz && mv tmp.vcf.gz ${sample_id}.region-${region_id}.vcf.gz
    bcftools index ${sample_id}.region-${region_id}.vcf.gz

    """

    stub:
    """
    touch ${sample_id}.region-${region_id}.vcf.gz
    touch ${sample_id}.region-${region_id}.vcf.gz.tbi
    """
}
