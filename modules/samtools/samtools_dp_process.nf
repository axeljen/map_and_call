process samtools_dp {
    tag "samtools_dp"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), path(cram), path(crai), val(region_id), val(region)

    output:
    tuple val(sample_id), path("${region_id}_${sample_id}.depths.bed"), emit: region_dp

    script:
    def chrom = region.split(':')[0]
    def start = region.split(':')[1].split('-')[0]
    def end = region.split(':')[1].split('-')[1]
    def region_length = end.toInteger() - start.toInteger() + 1
    
    """
    samtools depth -r ${region} -Q 30 -q 20 -a ${cram} | \
            awk -v OFS='\t' ' {print \$1 OFS \$2 - 1 OFS \$2 OFS \$3} ' | \
            bedtools groupby -i - -g 1,4 -c 2,3 -o min,max | \
            awk -v OFS='\t' ' { print \$1 OFS \$3 OFS \$4 OFS \$2} ' > ${region_id}_${sample_id}.depths.bed
    
    """

    stub:
    """
    touch ${region_id}_${sample_id}.depths.bed
    """
}
