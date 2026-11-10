process parse_region_depths {
    tag "parse_depths"
    label 'thin_short'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), path(bedfiles)
    each path(faidx)

    output:
    tuple val(sample_id), path("${sample_id}.depths.bed.gz"), emit: sample_depth_beds
    tuple val(sample_id), path("${sample_id}.depths.avg.txt"), emit: sample_depth_avg
    script:
    // Sort by numeric prefix before first underscore
    def sortedBeds = bedfiles.sort { a, b ->
        a.name.tokenize('_')[0].toInteger() <=> b.name.tokenize('_')[0].toInteger()
    }
    
    """    
    set -euo pipefail
    zcat ${sortedBeds.join(' ')} | \
        awk -v avgfile="${sample_id}.depths.avg.txt" 'BEGIN {OFS="\t"}
        {
            chrom=\$1
            len=\$3-\$2
            total_bases[chrom] += len*\$4
            total_length[chrom] += len
            print
        } 
        END {
            for (chrom in total_length) {
                print chrom, total_length[chrom], total_bases[chrom]/total_length[chrom] > avgfile
            }
        }' | gzip -c > ${sample_id}.depths.bed.gz
    touch ${sample_id}.depths.avg.txt

    """

    stub:
    """
    touch ${sample_id}.depths.bed.gz
    touch ${sample_id}.depths.avg.txt
    """
}
