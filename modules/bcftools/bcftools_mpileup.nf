
process mpileup{
    tag "$reference"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(pop_id), path(cram), path(crai), path(reference), path(reference_indices), val(region_id), val(regions)

    output:
    tuple val(region_id), val(pop_id), path("${pop_id}.region-${region_id}.vcf.gz"), path("${pop_id}.region-${region_id}.vcf.gz.csi"), emit: vcf

    script:
    """
    bcftools mpileup -r ${regions} \
        --fasta-ref ${reference} \
        --threads ${task.cpus} \
        -Ou \
        -a FORMAT/DP,FORMAT/AD \
        ${cram} | \
    bcftools call -m --gvcf 0 --threads ${task.cpus} -Oz -o ${pop_id}.region-${region_id}.vcf.gz -
    bcftools index ${pop_id}.region-${region_id}.vcf.gz
    """

    stub:
    """
    touch ${pop_id}.region-${region_id}.vcf.gz
    touch ${pop_id}.region-${region_id}.vcf.gz.csi

    """

}