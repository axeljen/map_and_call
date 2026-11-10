process combine_bedfiles {
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
    source task_tmp.sh ${params.use_scratch} "\${TASK_SCRATCH_ROOT:-${params.use_scratch ? params.scratch_dir ?: '' : ''}}"
    for bedfile in ${bedfiles}; do
        cat \${bedfile} >> "\$TMPDIR/combined.bed"
    done
    bedtools sort -faidx ${reference_fai} -i "\$TMPDIR/combined.bed" > ${sample_id}_${name}.bed
    # print the sum of all regions to stdout
    awk '{sum += \$3 - \$2} END {print sum}' ${sample_id}_${name}.bed > total_sites.txt
    """

    stub:
    """
    touch ${sample_id}_${name}.bed
    echo "0" > total_sites.txt
    """
}

process concatenate_masks {
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
