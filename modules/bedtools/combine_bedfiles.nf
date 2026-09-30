process combine_bedfiles {
    scratch params.use_scratch
    tag "combine_bedfiles"
    label 'medium_short'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/03_genotypes", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bedfiles), path(reference_fai), val(name)
    
    output:
    tuple val(sample_id), path("${sample_id}_${name}.bed"), emit: bedfile
    tuple val(sample_id), path("total_sites.txt"), emit: total_sites

    script:
    """
    for bedfile in ${bedfiles}; do
        cat \${bedfile} >> tmp.${sample_id}.bed
    done
    bedtools sort -faidx ${reference_fai} -i tmp.${sample_id}.bed > ${sample_id}_${name}.bed
    # print the sum of all regions to stdout
    awk '{sum += \$3 - \$2} END {print sum}' ${sample_id}_${name}.bed > total_sites.txt
    rm tmp.${sample_id}.bed
    """

    stub:
    """
    touch ${sample_id}_${name}.bed
    echo "0" > total_sites.txt
    """
}

process concatenate_masks {
    scratch params.use_scratch
    tag "concatenate_masks"
    label 'thin_short'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), path(homref_masks), path(total_masks), path(snp_masks)

    output:
    tuple val(sample_id), path("${sample_id}_homref_invariants.bed.gz"), emit: homref_mask
    tuple val(sample_id), path("${sample_id}_mappability_mask.bed.gz"), emit: mappability_mask
    tuple val(sample_id), path("${sample_id}_snp_mask.bed.gz"), emit: snp_mask

    script:
    """
    for maskfile in ${homref_masks}; do
        cat \${maskfile} >> ${sample_id}_homref_invariants.bed.gz
    done

    for maskfile in ${total_masks}; do
        cat \${maskfile} >> ${sample_id}_mappability_mask.bed.gz
    done

    for maskfile in ${snp_masks}; do
        cat \${maskfile} >> ${sample_id}_snp_mask.bed.gz
    done

    """

    stub:
    """
    touch ${sample_id}_homref_invariants.bed
    touch ${sample_id}_mappability_mask.bed
    touch ${sample_id}_snp_mask.bed
    """
}
