// Check 

process damage_profiler_rescale {
    tag "$sample_id"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), val(datatype), path(bam_file), path(bam_idx)
    each path(reference_genome)

    output:
    tuple val(sample_id), path("${sample_id}_damage"), emit: damage_reports
    tuple val(sample_id), path("${sample_id}.rescaled.dedup.bam"), path("${sample_id}.rescaled.dedup.bam.bai"), emit: rescaled_bam

    script:

    """
    source task_tmp.sh ${params.use_scratch} "\${TASK_SCRATCH_ROOT:-${params.use_scratch ? params.scratch_dir ?: '' : ''}}"

    ## Mapdamage command - use fixed output directory name for reproducibility
    mapDamage -i ${bam_file} -r ${reference_genome} -d "\$TMPDIR/results_${sample_id}" --rescale

    # move damage profiler results to a subdirectory with the sample name
    mkdir -p ${sample_id}_damage
    mv "\$TMPDIR/results_${sample_id}"/*.pdf ${sample_id}_damage/ 2>/dev/null || true

    # find the rescaled bam, sort and index
    samtools sort -@ ${task.cpus} -T "\$TMPDIR/sort" -o ${sample_id}.rescaled.dedup.bam "\$TMPDIR/results_${sample_id}"/*rescaled.bam
    samtools index ${sample_id}.rescaled.dedup.bam

    """

    stub:
    """
    mkdir -p ${sample_id}_damage
    touch ${sample_id}.rescaled.dedup.bam
    touch ${sample_id}.rescaled.dedup.bam.bai
    """
}