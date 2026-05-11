process parse_region_depths {
    tag "parse_depths"
    label 'thin_short'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), path(bedfiles)
    path faidx

    output:
    tuple val(sample_id), path("${sample_id}.depths.bed"), emit: sample_depth_beds
    tuple val(sample_id), path("${sample_id}.depths.avg.txt"), emit: sample_depth_avg
    script:
    // Sort by numeric prefix before first underscore
    def sortedBeds = bedfiles.sort { a, b ->
        a.name.tokenize('_')[0].toInteger() <=> b.name.tokenize('_')[0].toInteger()
    }
    
    """    
    # Cat all bed files together, add sort key, sort, then remove key
    for bed in ${sortedBeds}; do
        cat \$bed >> ${sample_id}.depths.bed
    done

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
