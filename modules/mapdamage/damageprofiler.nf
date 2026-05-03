// Check 

process damage_profiler {
    tag "$sample_id"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), path(bam_file), path(bam_idx)
    each path(reference_genome)

    output:
    tuple val(sample_id), path("${sample_id}_damage"), emit: damage_reports

    script:

    """

    ## Mapdamage command - use fixed output directory name for reproducibility
    mapDamage -i ${bam_file} -r ${reference_genome} -d results_${sample_id}

    # move damage profiler results to a subdirectory with the sample name
    mkdir -p ${sample_id}_damage
    mv results_${sample_id}/*.pdf ${sample_id}_damage/ 2>/dev/null || true

    """

    stub:
    """
    mkdir -p ${sample_id}_damage
    """
}