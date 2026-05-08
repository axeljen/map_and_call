process finalize_masks {
    tag "finalize_masks"
    label 'medium_short'
    conda "${moduleDir}/environment.yml"


    input: 
    tuple val(sample_id), val(region_id), val(regions), path(callable_bed), path(filtered_snps), path(filtered_indels)
    
    output:
    tuple val(sample_id), path("${sample_id}_${region_id}.homref_invariants.bed.gz"), emit: homref_invariants
    tuple val(sample_id), path("${sample_id}_${region_id}.mappability_mask_total.bed.gz"), emit: mappability_mask
    tuple val(sample_id), path("${sample_id}_${region_id}.mappability_mask_snps.bed.gz"), emit: mappability_mask_snps

    script:
    def region_list = regions.join(' ')
    """
    # make dummy bedfile for filtering input bed
    for region in ${region_list};
        do
        chrom=\$(echo \$region | cut -d: -f1)
        start=\$((\$(echo \$region | cut -d: -f2 | cut -d- -f1) - 1))
        end=\$(echo \$region | cut -d: -f2 | cut -d- -f2)
        echo -e "\${chrom}\t\${start}\t\${end}" >> ${sample_id}_${region_id}_regionstring.bed
    done
    bedtools intersect -a ${callable_bed} -b ${sample_id}_${region_id}_regionstring.bed > ${sample_id}_${region_id}_callable.bed.tmp && mv ${sample_id}_${region_id}_callable.bed.tmp ${sample_id}_${region_id}_callable.bed

    # now it's just a matter of subtracting some stuff from callable bed to get the different masks
    # homref invariants should only contain regions that are not represented in the vcf files, so subtract both the filtered snps and indels from the callable bed
    bedtools subtract -a ${sample_id}_${region_id}_callable.bed -b ${filtered_snps} | \
        bedtools subtract -a - -b ${filtered_indels} > ${sample_id}_${region_id}.homref_invariants.bed
    
    # mappability mask snps should contain all regions with high enough mapping coverage to call variants, but where there are no indels, so subtract those
    bedtools subtract -a ${sample_id}_${region_id}_callable.bed -b ${filtered_indels} > ${sample_id}_${region_id}.mappability_mask_snps.bed

    # mappability mask total should contain all regions with high enough mapping coverage to call variants, so just copy the callable bed and rename it
    cp ${sample_id}_${region_id}_callable.bed ${sample_id}_${region_id}.mappability_mask_total.bed

    # gzip the output files
    gzip ${sample_id}_${region_id}.homref_invariants.bed
    gzip ${sample_id}_${region_id}.mappability_mask_total.bed
    gzip ${sample_id}_${region_id}.mappability_mask_snps.bed

    """

    stub:
    """

    """

}
