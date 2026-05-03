// Check 

process damage_profiler_rescale {
    tag "$sample_id"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), val(datatype), path(bam_file), path(bam_idx)
    each path(reference_genome)

    output:
    tuple val(sample_id), path("${sample_id}_damage"), emit: damage_reports
    tuple val(sample_id), path("${sample_id}.rescaled.dedup.bam"), path("${sample_id}.rescaled.dedup.bam.bai"), emit: rescaled_bam

    script:

    """
    
    ## Mapdamage command - use fixed output directory name for reproducibility
    mapDamage -i ${bam_file} -r ${reference_genome} -d results_${sample_id} --rescale

    # move damage profiler results to a subdirectory with the sample name
    mkdir -p ${sample_id}_damage
    mv results_${sample_id}/*.pdf ${sample_id}_damage/ 2>/dev/null || true

    # find the rescaled bam, sort and index
    samtools sort -@ ${task.cpus} -o ${sample_id}.rescaled.dedup.bam results_${sample_id}/*rescaled.bam
    samtools index ${sample_id}.rescaled.dedup.bam

    """

    stub:
    """
    mkdir -p ${sample_id}_damage
    touch ${sample_id}.rescaled.dedup.bam
    touch ${sample_id}.rescaled.dedup.bam.bai
    """
}