
process mpileup{
    tag "$reference"
    conda "${moduleDir}/environment.yml"

    input:
    tuple path(cram), path(crai), path(reference), path(reference_indices), val(region_id), val(regions), val(pops)

    output:
    tuple val(region_id), path("bcftools.region-${region_id}.vcf.gz"), path("bcftools.region-${region_id}.vcf.gz.csi"), emit: vcf

    script:
    def regions_list = regions.join(',')
    def bamstring = cram.join(' ')
    """
    # Create popfile from the Groovy-generated string
    echo -e "${pops.join('\n')}" | sed 's/=/\t/' > bcftools_popfile.txt

    bcftools mpileup -r ${regions_list} \
        --fasta-ref ${reference} \
        --threads ${task.cpus} \
        -Ou \
        -a FORMAT/DP,FORMAT/AD \
        ${bamstring} | \
    bcftools call -m -v -G bcftools_popfile.txt --threads ${task.cpus} -Oz -o bcftools.region-${region_id}.vcf.gz -
    bcftools index bcftools.region-${region_id}.vcf.gz
    """

    stub:
    """
    touch bcftools.region-${region_id}.vcf.gz
    touch bcftools.region-${region_id}.vcf.gz.csi

    """

}