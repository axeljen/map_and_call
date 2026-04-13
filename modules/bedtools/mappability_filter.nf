process mappability_filter {
    tag "mappability_filter"
    label 'process_small'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/mappability_filtered_vcf_files", mode: 'copy'

    input:
    tuple val(region_id), val(sample), path(vcf), path(csi), path(mappability_bed), val(category)

    output:
    tuple val(region_id), val(sample), path("${category}_${sample}_region-${region_id}.mappability_filtered.vcf.gz"), path("${category}_${sample}_region-${region_id}.mappability_filtered.vcf.gz.*"), val(category), emit: vcf

    script:
    """
    bcftools view -s ${sample} -Ov ${vcf} | \
    bedtools intersect -header -a /dev/stdin -b ${mappability_bed} | bgzip -c > ${category}_${sample}_region-${region_id}.mappability_filtered.vcf.gz
    bcftools index ${category}_${sample}_region-${region_id}.mappability_filtered.vcf.gz
    """

    stub:
    """
    touch ${category}_${sample}_region-${region_id}.mappability_filtered.vcf.gz
    touch ${category}_${sample}_region-${region_id}.mappability_filtered.vcf.gz.csi
    """
}
