
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
    # before indexing check that vcf actually has records, otherwise bcftools index will fail
    if [ \$(bcftools view -H ${category}.vcf.gz | wc -l) -gt 0 ]; then
        bcftools index ${category}.vcf.gz
    else
        echo "Warning: vcf file ${category}.vcf.gz is empty."
    fi
    """

    stub:
    """
    touch ${category}.vcf.gz
    touch ${category}.vcf.gz.csi
    """
}