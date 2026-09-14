process samtools_dp_deprecated {
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

process samtools_dp {
    scratch params.use_scratch
    tag "samtools_dp"
    
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(region_id), val(regions), val(sample_ids), path(crams), path(crais) 

    output:
    tuple val(region_id),
        val(sample_ids),
        path("*.depths.bed.gz"),
        emit: region_dp

    script:
    def region_list = regions.join(' ')
    def sample_list = sample_ids.join(' ')
    def cram_list = crams.join(' ')
    def crai_list = crais.join(' ')
    
    """
    set -o pipefail
    sample_ids=(${sample_list})
    crams=(${cram_list})
    crais=(${crai_list})
    for i in "\${!sample_ids[@]}"; do
        sample_id="\${sample_ids[\$i]}"
        cram="\${crams[\$i]}"
        crai="\${crais[\$i]}"
        for region in ${region_list};
        do
        samtools depth -r \${region} -Q ${params.min_mapqual} -q ${params.min_basequal} -a \${cram} | \
                awk -v OFS='\t' ' {print \$1 OFS \$2 - 1 OFS \$2 OFS \$3} ' | \
                bedtools groupby -i - -g 1,4 -c 2,3 -o min,max | \
                awk -v OFS='\t' ' { print \$1 OFS \$3 OFS \$4 OFS \$2} ' >> ${region_id}_\${sample_id}.depths.bed
        done
        # zip it to save space
        gzip -f ${region_id}_\${sample_id}.depths.bed
    done

    

    """

    stub:
    """
    touch ${region_id}_\${sample_id}.depths.bed
    """
}
