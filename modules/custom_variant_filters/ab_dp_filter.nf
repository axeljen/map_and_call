process ab_dp_filter {
    tag "ab_dp_filter"
    label 'medium_medium'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/ab_dp_filtered_vcf_files", mode: 'copy'

    input:
    tuple val(sample), val(region_id), val(region), path(callable_regions), path(vcf), path(csi), val(min_depth), val(max_depth), val(sex_assignment)
    val sex_linked_scaffolds
    val category

    output:
    tuple val(region_id), val(sample), path("${category}_${sample}_region-${region_id}.ab_dp_filtered.vcf.gz"), path("${category}_${sample}_region-${region_id}.ab_dp_filtered.vcf.gz.*"), emit: ab_dp_filtered_vcf

    script:
    def chrom = region.split(':')[0]
    def start = region.split(':')[1].split('-')[0].toInteger() - 1 // convert to 0-based for bedtools
    def end = region.split(':')[1].split('-')[1]
    def sexlinked_arg = sex_linked_scaffolds ? "--sex-linked-scaffolds ${sex_linked_scaffolds}" : ""
    """
    source task_tmp.sh ${params.use_scratch} "\${TASK_SCRATCH_ROOT:-${params.use_scratch ? params.scratch_dir ?: '' : ''}}"
    # first run the vcf file through bcftools to remove uncallable regions and pull out the sample of interest
    # extract callable regions for the specific region
    echo -e "${chrom}\t${start}\t${end}" > "\$TMPDIR/regionstring.bed"
    # and filter the bed file to only hold this region
    bedtools intersect -a ${callable_regions} -b "\$TMPDIR/regionstring.bed" > "\$TMPDIR/callable.bed"


    bcftools view -R "\$TMPDIR/callable.bed" -s ${sample} -Oz ${vcf} |\
    bcftools view -a -Oz -o "\$TMPDIR/input.vcf.gz"
    bcftools index "\$TMPDIR/input.vcf.gz"

    # then push it through the allel balance filtration
    ab_dp_filtration.py -i "\$TMPDIR/input.vcf.gz" -o ${category}_${sample}_region-${region_id}.ab_dp_filtered.vcf.gz \
        --min-ab ${params.min_allele_balance} --samples ${sample} \
        --min-depth ${min_depth} --max-depth ${max_depth} --sex-assignments ${sex_assignment} ${sexlinked_arg}
    # Try to index; if it fails (e.g., empty after filtering), create a valid empty VCF with header
    if ! bcftools index ${category}_${sample}_region-${region_id}.ab_dp_filtered.vcf.gz 2>/dev/null; then
        echo "Warning: Could not index ${category}_${sample}_region-${region_id}.ab_dp_filtered.vcf.gz (likely empty after filtering)"
        bcftools view -h "\$TMPDIR/input.vcf.gz" | bgzip > ${category}_${sample}_region-${region_id}.ab_dp_filtered.vcf.gz
        bcftools index ${category}_${sample}_region-${region_id}.ab_dp_filtered.vcf.gz
    fi
    """

    stub:
    """
    touch ${category}_${sample}_region-${region_id}.ab_dp_filtered.vcf.gz
    touch ${category}_${sample}_region-${region_id}.ab_dp_filtered.vcf.gz.csi
    """

}
