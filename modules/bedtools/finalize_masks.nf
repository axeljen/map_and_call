process finalize_masks_legacy {
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
    source task_tmp.sh ${params.use_scratch} "\${TASK_SCRATCH_ROOT:-${params.use_scratch ? params.scratch_dir ?: '' : ''}}"
    # make dummy bedfile for filtering input bed
    for region in ${region_list};
        do
        chrom=\$(echo \$region | cut -d: -f1)
        start=\$((\$(echo \$region | cut -d: -f2 | cut -d- -f1) - 1))
        end=\$(echo \$region | cut -d: -f2 | cut -d- -f2)
        echo -e "\${chrom}\t\${start}\t\${end}" >> "\$TMPDIR/regions.bed"
    done
    bedtools intersect -a ${callable_bed} -b "\$TMPDIR/regions.bed" > "\$TMPDIR/callable.bed"

    # now it's just a matter of subtracting some stuff from callable bed to get the different masks
    # homref invariants should only contain regions that are not represented in the vcf files, so subtract both the filtered snps and indels from the callable bed
    bedtools subtract -a "\$TMPDIR/callable.bed" -b ${filtered_snps} | \
        bedtools subtract -a - -b ${filtered_indels} | gzip -c > ${sample_id}_${region_id}.homref_invariants.bed.gz
    
    # mappability mask snps should contain all regions with high enough mapping coverage to call variants, but where there are no indels, so subtract those
    bedtools subtract -a "\$TMPDIR/callable.bed" -b ${filtered_indels} | gzip -c > ${sample_id}_${region_id}.mappability_mask_snps.bed.gz

    # mappability mask total should contain all regions with high enough mapping coverage to call variants, so just copy the callable bed and rename it
    gzip -c "\$TMPDIR/callable.bed" > ${sample_id}_${region_id}.mappability_mask_total.bed.gz

    """

    stub:
    """

    """

}

process finalize_masks {
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
    source task_tmp.sh ${params.use_scratch} "\${TASK_SCRATCH_ROOT:-${params.use_scratch ? params.scratch_dir ?: '' : ''}}"
    # make dummy bedfile for filtering input bed
    for region in ${region_list};
        do
        chrom=\$(echo \$region | cut -d: -f1)
        start=\$((\$(echo \$region | cut -d: -f2 | cut -d- -f1) - 1))
        end=\$(echo \$region | cut -d: -f2 | cut -d- -f2)
        echo -e "\${chrom}\t\${start}\t\${end}" >> "\$TMPDIR/regions.bed"
    done

    # iterate through all the samples
    sample_list=(${sample_names})
    bedlist=(${bedfile_paths})

    for i in \${!sample_list[@]}; do
        sample_id="\${sample_list[\$i]}"
        callable_bed="\${bedlist[\$i]}"
        
        # subset the vcf file to the sample and remove missing data so we're only subtracting for this specific sample
        bcftools view -s "\${sample_id}" ${filtered_snps} | bcftools view -e 'GT="./."' -Oz -o "\$TMPDIR/filtered_snps.vcf.gz"
        bcftools view -s "\${sample_id}" ${filtered_indels} | bcftools view -e 'GT="./."' -Oz -o "\$TMPDIR/filtered_indels.vcf.gz"
        

        bedtools intersect -a \${callable_bed} -b "\$TMPDIR/regions.bed" > "\$TMPDIR/callable.bed"

        # now it's just a matter of subtracting some stuff from callable bed to get the different masks
        # homref invariants should only contain regions that are not represented in the vcf files, so subtract both the filtered snps and indels from the callable bed
        bedtools subtract -a "\$TMPDIR/callable.bed" -b "\$TMPDIR/filtered_snps.vcf.gz" | \
        bedtools subtract -a - -b "\$TMPDIR/filtered_indels.vcf.gz" | gzip -c > \${sample_id}_${region_id}.homref_invariants.bed.gz
    
        # mappability mask snps should contain all regions with high enough mapping coverage to call variants, but where there are no indels, so subtract those
        bedtools subtract -a "\$TMPDIR/callable.bed" -b "\$TMPDIR/filtered_indels.vcf.gz" | gzip -c > \${sample_id}_${region_id}.mappability_mask_snps.bed.gz

        # mappability mask total should contain all regions with high enough mapping coverage to call variants, so just copy the callable bed and rename it
        gzip -c "\$TMPDIR/callable.bed" > \${sample_id}_${region_id}.mappability_mask_total.bed.gz

    done

    """

    stub:
    """

    touch \${sample_id}_${region_id}.homref_invariants.bed.gz
    touch \${sample_id}_${region_id}.mappability_mask_total.bed.gz
    touch \${sample_id}_${region_id}.mappability_mask_snps.bed.gz

    """

}
