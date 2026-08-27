
process bcftools_concat {
    scratch params.use_scratch
    tag "concat_vcfs"
    label 'thin_medium'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/03_genotypes", mode: 'copy'

    input:
    tuple path(vcf_files), path(vcf_indices)
    val category

    output:
    tuple path("${category}.vcf.gz"), path("${category}.vcf.gz.*"), emit: vcf

    script:
    """
    bcftools concat -Oz -o ${category}.vcf.gz ${vcf_files.join(' ')}
    # Try to index; if it fails (e.g., empty after concat), create a valid empty VCF with header
    if ! bcftools index ${category}.vcf.gz 2>/dev/null; then
        echo "Warning: Could not index ${category}.vcf.gz (likely empty after concatenation)"
        # Get header from first input file and create proper empty VCF
        bcftools view -h ${vcf_files[0]} | bgzip > ${category}.vcf.gz
        bcftools index ${category}.vcf.gz
    fi
    """

    stub:
    """
    touch ${category}.vcf.gz
    touch ${category}.vcf.gz.csi
    """
}