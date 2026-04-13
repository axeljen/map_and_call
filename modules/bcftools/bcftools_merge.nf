
process bcftools_merge {
    tag "$region_id"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"


    input:
    tuple val(region_id), val(sample_ids), path(vcf_files), path(vcf_csi), path(reference_fasta), path(reference_fasta_index), path(reference_fasta_gzi)
    val(category)

    output:
    tuple val(region_id), path("region-${region_id}.${category}.vcf.gz"), path("region-${region_id}.${category}.vcf.gz.*"), emit: vcf

    script:
    """
    bcftools merge -g ${reference_fasta} -Ov ${vcf_files.join(' ')} |\
        # remove any invariant sites from the merged vcf
        bcftools view -m2 -Oz -o region-${region_id}.${category}.vcf.gz -
    bcftools index region-${region_id}.${category}.vcf.gz
    """

    stub:
    """
    touch region-${region_id}.${category}.vcf.gz
    touch region-${region_id}.${category}.vcf.gz.csi
    """
}