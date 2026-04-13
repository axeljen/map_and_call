process bcftools_filter_gatkindels {
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
    bcftools index region-${region_id}.indels.filtered1.vcf.gz
    """

    stub:
    """
    touch region-${region_id}.indels.filtered1.vcf.gz
    touch region-${region_id}.indels.filtered1.vcf.gz.csi
    """
}
