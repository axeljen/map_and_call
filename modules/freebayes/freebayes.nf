
process freebayes{
    tag "freebayes_${region_id}_${pop_id}"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/sample_vcfs", mode: 'copy'

    input:
    tuple val(pop_id), path(cram), path(crai), path(reference), path(reference_indices), val(region_id), val(regions), val(pops)

    output:
    tuple val(region_id), val(pop_id), path("${pop_id}.region-${region_id}.vcf.gz"), path("${pop_id}.region-${region_id}.vcf.gz.csi"), emit: vcf
    
    script:
    """
    # Create popfile from the Groovy-generated string
    
    echo -e "${pops.join('\n')}" | sed 's/=/\t/' > freebayes_popfile.txt

    # freebayes uses 0-based coordinates, so we need to convert the regions from 1-based to 0-based
    chrom=\$(echo ${regions} | cut -d: -f1)
    start=\$(echo ${regions} | cut -d: -f2 | cut -d- -f1)
    end=\$(echo ${regions} | cut -d: -f2 | cut -d- -f2)
    start_0_based=\$((start - 1))
    regions_0_based="\${chrom}:\${start_0_based}-\${end}"

    freebayes -f ${reference} -r \${regions_0_based} --populations freebayes_popfile.txt --ploidy ${params.ploidy} ${cram} | bgzip - > ${pop_id}.region-${region_id}.vcf.gz
    bcftools index ${pop_id}.region-${region_id}.vcf.gz
    """

    stub:
    """
    touch ${pop_id}.region-${region_id}.vcf.gz
    touch ${pop_id}.region-${region_id}.vcf.gz.csi

    """

}