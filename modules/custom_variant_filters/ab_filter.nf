process ab_filter {
    scratch params.use_scratch
    tag "ab_filter"
    label 'thin_medium'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/ab_filtered_vcf_files", mode: 'copy'

    input:
    tuple val(region_id), path(vcf), path(csi), val(category)

    output:
    tuple val(region_id), path("${category}_region-${region_id}.ab_filtered.vcf.gz"), path("${category}_region-${region_id}.ab_filtered.vcf.gz.*"), emit: vcf

    script:
    """
    ab_filtration.py -i ${vcf} -o ${category}_region-${region_id}.ab_filtered.vcf.gz --min-ab ${params.min_allele_balance}
    # Try to index; if it fails (e.g., empty VCF), create a valid empty index
    if ! bcftools index ${category}_region-${region_id}.ab_filtered.vcf.gz 2>/dev/null; then
        echo "Warning: Could not index ${category}_region-${region_id}.ab_filtered.vcf.gz (likely empty)"
        # Create a proper empty VCF with header and index it
        bcftools view -h ${vcf} | bgzip > ${category}_region-${region_id}.ab_filtered.vcf.gz
        bcftools index ${category}_region-${region_id}.ab_filtered.vcf.gz
    fi
    """

    stub:
    """
    touch ${category}_region-${region_id}.ab_filtered.vcf.gz
    touch ${category}_region-${region_id}.ab_filtered.vcf.gz.csi
    """
}
