
process freebayes{
    tag "freebayes_${region_id}"
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/sample_vcfs", mode: 'copy'

    input:
    tuple path(cram), path(crai), path(reference), path(reference_indices), val(region_id), val(regions), val(pops)

    output:
    tuple val(region_id), path("freebayes.region-${region_id}.vcf.gz"), path("freebayes.region-${region_id}.vcf.gz.csi"), emit: vcf
    
    script:
    def regions_list = regions.join(' ')
    def bamlist = cram.join(' ')
    """
    # Create popfile from the Groovy-generated string
    
    echo -e "${pops.join('\n')}" | sed 's/=/\t/' > freebayes_popfile.txt

    # make a target bedfile with the regions for freebayes to use
    for region in ${regions_list};
        do
        chrom=\$(echo \$region | cut -d: -f1)
        # convert start to zero based
        start=\$((\$(echo \$region | cut -d: -f2 | cut -d- -f1) - 1))
        end=\$(echo \$region | cut -d: -f2 | cut -d- -f2)
        echo -e "\${chrom}\t\${start}\t\${end}" >> freebayes_targets.bed
    done

    freebayes -f ${reference} -t freebayes_targets.bed --populations freebayes_popfile.txt --min-mapping-quality ${params.min_mapqual} --min-base-quality ${params.min_basequal} --ploidy ${params.ploidy} ${bamlist} | bgzip - > freebayes.region-${region_id}.vcf.gz
    bcftools index freebayes.region-${region_id}.vcf.gz
    """

    stub:
    """
    touch freebayes.region-${region_id}.vcf.gz
    touch freebayes.region-${region_id}.vcf.gz.csi

    """

}