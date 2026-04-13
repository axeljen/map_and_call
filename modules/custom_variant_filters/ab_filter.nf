process ab_filter {
    tag "ab_filter"
    label 'process_small'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/ab_filtered_vcf_files", mode: 'copy'

    input:
    tuple val(region_id), val(sample), path(vcf), path(csi), val(category)

    output:
    tuple val(region_id), val(sample), path("${category}_${sample}_region-${region_id}.ab_filtered.vcf.gz"), path("${category}_${sample}_region-${region_id}.ab_filtered.vcf.gz.*"), emit: vcf

    script:
    """
    ab_filtration.py -i ${vcf} -o ${category}_${sample}_region-${region_id}.ab_filtered.vcf.gz --min-ab ${params.min_allele_balance}
    bcftools index ${category}_${sample}_region-${region_id}.ab_filtered.vcf.gz
    """

    stub:
    """
    touch ${category}_${sample}_region-${region_id}.ab_filtered.vcf.gz
    touch ${category}_${sample}_region-${region_id}.ab_filtered.vcf.gz.csi
    """
}
