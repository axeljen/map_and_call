process samtools_dp {
    scratch params.use_scratch
    tag "samtools_dp"
    
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), path(cram), path(crai), val(region_id), val(regions)

    output:
    tuple val(sample_id), path("${region_id}_${sample_id}.depths.bed.gz"), emit: region_dp

    script:
    def region_list = regions.join(' ')
    
    """
    set -o pipefail
    
    for region in ${region_list};
        do
        samtools depth -r \${region} -Q ${params.min_mapqual} -q ${params.min_basequal} -a ${cram} | \
                awk -v OFS='\t' ' {print \$1 OFS \$2 - 1 OFS \$2 OFS \$3} ' | \
                bedtools groupby -i - -g 1,4 -c 2,3 -o min,max | \
                awk -v OFS='\t' ' { print \$1 OFS \$3 OFS \$4 OFS \$2} ' >> ${region_id}_${sample_id}.depths.bed
    done

    # zip it to save space
    gzip -f ${region_id}_${sample_id}.depths.bed

    """

    stub:
    """
    touch ${region_id}_${sample_id}.depths.bed
    """
}
