process bcftools_filter_gatkindels {
    scratch params.use_scratch
    tag "filter_indels"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/filtered_indel_1_vcf_files", mode: 'copy'

    input:
    tuple val(region_id), path(vcf), path(csi)

    output:
    tuple val(region_id), path("region-${region_id}.indels.filtered1.vcf.gz"), path("region-${region_id}.indels.filtered1.vcf.gz.*"), emit: vcf

    script:
    """
    bcftools filter -e ' QUAL<30 || MQ<40 || FS>200.0 || ReadPosRankSum<-20.0' -Oz -o region-${region_id}.indels.filtered1.vcf.gz ${vcf}
    # Try to index; if it fails (e.g., empty after filtering), create a valid empty VCF with header
    if ! bcftools index region-${region_id}.indels.filtered1.vcf.gz 2>/dev/null; then
        echo "Warning: Could not index region-${region_id}.indels.filtered1.vcf.gz (likely empty after filtering)"
        bcftools view -h ${vcf} | bgzip > region-${region_id}.indels.filtered1.vcf.gz
        bcftools index region-${region_id}.indels.filtered1.vcf.gz
    fi
    """

    stub:
    """
    touch region-${region_id}.indels.filtered1.vcf.gz
    touch region-${region_id}.indels.filtered1.vcf.gz.csi
    """
}
