process callability_filter {
    tag "vcf_depth_filter_sample"
    conda "${moduleDir}/environment.yml"
    label "thin_medium"

    input:
    tuple val(region_id), val(regions), path(vcf), path(csi), val(sample), path(callable_regions)

    output:
    tuple val(region_id), val(sample), path("${sample}_${region_id}.dp.filtered.vcf.gz"), path("${sample}_${region_id}.dp.filtered.vcf.gz.*"), emit: vcf

    script:
    def region_list = regions.join(' ')
    """
    # extract only the focal sample
    bcftools view -s ${sample} -Oz -o ${sample}_${region_id}.vcf.gz ${vcf}
    bcftools index ${sample}_${region_id}.vcf.gz
    # grab the callable regions from the specific region
    for region in ${region_list};
        do
        chrom=\$(echo \$region | cut -d: -f1)
        startpos=\$((\$(echo \$region | cut -d: -f2 | cut -d- -f1) - 1))
        endpos=\$(echo \$region | cut -d: -f2 | cut -d- -f2)
        echo -e "\${chrom}\t\${startpos}\t\${endpos}" >> tmp.${sample}.${region_id}.bed
    done
    bedtools intersect -a ${callable_regions} -b tmp.${sample}.${region_id}.bed > tmp.${sample}.${region_id}.callable.bed

    # empty bed files will break bcftools view, so we need to add a dummy region if there are no callable regions for this sample/region
    if [ ! -s tmp.${sample}.${region_id}.callable.bed ]; then
        echo -e "\${chrom}\t0\t0" > tmp.${sample}.${region_id}.callable.bed
    fi

    # keep only genotypes that are within the callable regions
    bcftools view -R tmp.${sample}.${region_id}.callable.bed -Oz -o ${sample}_${region_id}.dp.filtered.vcf.gz ${sample}_${region_id}.vcf.gz
    bcftools index ${sample}_${region_id}.dp.filtered.vcf.gz
    
    """

    stub:
    """
    touch ${sample}_${region_id}.dp.filtered.vcf.gz
    touch ${sample}_${region_id}.dp.filtered.vcf.gz.csi
    """
}
