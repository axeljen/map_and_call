process finalize_masks {
    tag "finalize_masks"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"


    input: 
    tuple val(sample_id), val(region_id), val(region), path(callable_bed), path(filtered_snps), path(filtered_indels)
    
    output:
    tuple val(sample_id), path("${sample_id}_${region_id}.homref_invariants.bed"), emit: homref_invariants
    tuple val(sample_id), path("${sample_id}_${region_id}.mappability_mask_total.bed"), emit: mappability_mask
    tuple val(sample_id), path("${sample_id}_${region_id}.mappability_mask_snps.bed"), emit: mappability_mask_snps

    script:
    def chrom = region.split(':')[0]
    def start = region.split(':')[1].split('-')[0].toInteger() - 1 // convert to 0-based for bedtools
    def end = region.split(':')[1].split('-')[1]
    """
    # make dummy bedfile for filtering input bed
    echo -e ${chrom}"\t"${start}"\t"${end} > ${sample_id}_${region_id}_regionstring.bed
    bedtools intersect -a ${callable_bed} -b ${sample_id}_${region_id}_regionstring.bed > tmp && mv tmp ${sample_id}_${region_id}_callable.bed

    # now it's just a matter of subtracting some stuff from callable bed to get the different masks
    # homref invariants should only contain regions that are not represented in the vcf files, so subtract both the filtered snps and indels from the callable bed
    bedtools subtract -a ${sample_id}_${region_id}_callable.bed -b ${filtered_snps} | \
        bedtools subtract -a - -b ${filtered_indels} > ${sample_id}_${region_id}.homref_invariants.bed
    
    # mappability mask snps should contain all regions with high enough mapping coverage to call variants, but where there are no indels, so subtract those
    bedtools subtract -a ${sample_id}_${region_id}_callable.bed -b ${filtered_indels} > ${sample_id}_${region_id}.mappability_mask_snps.bed

    # mappability mask total should contain all regions with high enough mapping coverage to call variants, so just copy the callable bed and rename it
    cp ${sample_id}_${region_id}_callable.bed ${sample_id}_${region_id}.mappability_mask_total.bed

    """

    stub:
    """

    """

}
