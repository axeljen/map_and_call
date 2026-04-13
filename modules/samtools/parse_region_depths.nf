process parse_region_depths {
    tag "parse_depths"
    label 'process_small'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), path(bedfiles)
    path faidx

    output:
    tuple val(sample_id), path("${sample_id}.depths.bed"), emit: sample_depth_beds
    tuple val(sample_id), path("${sample_id}.depths.avg.txt"), emit: sample_depth_avg
    script:
    
    """
    # Extract chromosome order from faidx and create a sort key file
    awk '{print \$1}' ${faidx} | awk '{print \$1, NR}' > chrom_order.txt
    
    # Cat all bed files together, add sort key, sort, then remove key
    cat ${bedfiles} | \
        awk 'NR==FNR {order[\$1]=\$2; next} {print order[\$1], \$0}' chrom_order.txt - | \
        sort -k1,1n -k3,3n | \
        cut -d' ' -f2- > ${sample_id}.depths.bed

   # then do the average depths per scaffold
    cat ${sample_id}.depths.bed | \
        awk 'BEGIN {OFS="\t"} 
        {
            chrom=\$1
            len=\$3-\$2
            total_bases[chrom] += len*\$4
            total_length[chrom] += len
        } 
        END {
            for (chrom in total_length) {
                print chrom, total_length[chrom], total_bases[chrom]/total_length[chrom]
            }
        }' \
        > ${sample_id}.depths.avg.txt
    
    """

    stub:
    """
    touch ${sample_id}.depths.bed
    touch ${sample_id}.depths.avg.txt
    """
}
