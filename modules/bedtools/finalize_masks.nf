process finalize_masks_legacy {
    scratch params.use_scratch
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

    # delete temporary files
    rm ${sample_id}_${region_id}_callable.bed

    """

    stub:
    """

    """

}

process finalize_masks {
    scratch params.use_scratch
    tag "finalize_masks"
    label 'medium_short'
    conda "${moduleDir}/environment.yml"


    input: 
    tuple val(samples), 
    val(region_id), 
    val(regions), 
    path(bedfiles), 
    path(filtered_snps), 
    path(filtered_snps_index), 
    path(filtered_indels), 
    path(filtered_indels_index)
    
    output:
    tuple val(samples), val(region_id), path("*${region_id}.homref_invariants.bed.gz"), emit: homref_invariants
    tuple val(samples), val(region_id), path("*${region_id}.mappability_mask_total.bed.gz"), emit: mappability_mask
    tuple val(samples), val(region_id), path("*${region_id}.mappability_mask_snps.bed.gz"), emit: mappability_mask_snps

    script:
    def region_list = regions.join(' ')
    def sample_names = samples.join(' ')
    def bedfile_paths = bedfiles.join(' ')
    """
    # make dummy bedfile for filtering input bed
    for region in ${region_list};
        do
        chrom=\$(echo \$region | cut -d: -f1)
        start=\$((\$(echo \$region | cut -d: -f2 | cut -d- -f1) - 1))
        end=\$(echo \$region | cut -d: -f2 | cut -d- -f2)
        echo -e "\${chrom}\t\${start}\t\${end}" >> ${region_id}_regionstring.bed
    done

    # iterate through all the samples
    sample_list=(${sample_names})
    bedlist=(${bedfile_paths})

    for i in \${!sample_list[@]}; do
        sample_id="\${sample_list[\$i]}"
        callable_bed="\${bedlist[\$i]}"
        
        # subset the vcf file to the sample and remove missing data so we're only subtracting for this specific sample
        bcftools view -s "\${sample_id}" ${filtered_snps} | bcftools view -e 'GT="./."' -Oz -o \${sample_id}_${region_id}_filtered_snps.vcf.gz
        bcftools view -s "\${sample_id}" ${filtered_indels} | bcftools view -e 'GT="./."' -Oz -o \${sample_id}_${region_id}_filtered_indels.vcf.gz
        

        bedtools intersect -a \${callable_bed} -b ${region_id}_regionstring.bed > \${sample_id}_${region_id}_callable.bed.tmp && mv \${sample_id}_${region_id}_callable.bed.tmp \${sample_id}_${region_id}_callable.bed

        # now it's just a matter of subtracting some stuff from callable bed to get the different masks
        # homref invariants should only contain regions that are not represented in the vcf files, so subtract both the filtered snps and indels from the callable bed
        bedtools subtract -a \${sample_id}_${region_id}_callable.bed -b \${sample_id}_${region_id}_filtered_snps.vcf.gz | \
        bedtools subtract -a - -b \${sample_id}_${region_id}_filtered_indels.vcf.gz > \${sample_id}_${region_id}.homref_invariants.bed
    
        # mappability mask snps should contain all regions with high enough mapping coverage to call variants, but where there are no indels, so subtract those
        bedtools subtract -a \${sample_id}_${region_id}_callable.bed -b \${sample_id}_${region_id}_filtered_indels.vcf.gz > \${sample_id}_${region_id}.mappability_mask_snps.bed

        # mappability mask total should contain all regions with high enough mapping coverage to call variants, so just copy the callable bed and rename it
        mv \${sample_id}_${region_id}_callable.bed \${sample_id}_${region_id}.mappability_mask_total.bed

    # gzip the output files
    gzip \${sample_id}_${region_id}.homref_invariants.bed
    gzip \${sample_id}_${region_id}.mappability_mask_total.bed
    gzip \${sample_id}_${region_id}.mappability_mask_snps.bed
    rm \${sample_id}_${region_id}_filtered_snps.vcf.gz
    rm \${sample_id}_${region_id}_filtered_indels.vcf.gz

    done

    rm ${region_id}_regionstring.bed

    """

    stub:
    """

    touch \${sample_id}_${region_id}.homref_invariants.bed.gz
    touch \${sample_id}_${region_id}.mappability_mask_total.bed.gz
    touch \${sample_id}_${region_id}.mappability_mask_snps.bed.gz

    """

}
