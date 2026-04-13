process select_invariant {
    tag "select_invariant"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    // publishDir "${params.outdir}/invariant_vcf_files", mode: 'copy'

    input:
    tuple val(region_id), path(vcf), path(csi)

    output:
    tuple val(region_id), path("region-${region_id}.invariant.vcf.gz"), path("region-${region_id}.invariant.vcf.gz.*"), emit: vcf

    script:
    """
    bcftools view -m 1 -M 1 -Oz -o region-${region_id}.invariant.vcf.gz ${vcf}
    bcftools index region-${region_id}.invariant.vcf.gz
    """

    stub:
    """
    touch region-${region_id}.invariant.vcf.gz
    touch region-${region_id}.invariant.vcf.gz.csi
    """
}
