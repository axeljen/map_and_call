// Check 

process damage_profiler {
    tag "$sample_id"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), path(bam_file), path(bam_idx)
    path reference_genome

    output:
    tuple val(sample_id), path("${sample_id}_damage"), emit: damage_reports
    tuple val(sample_id), path("${sample_id}.dedup.cram"), path("${sample_id}.dedup.cram.crai"), emit: rescaled_bam

    script:
    def rescale_arg = params.damageprofiler_rescale ? '--rescale' : ''
    def input_bam = bam_file.name.endsWith('.cram') ? "${sample_id}_input.bam" : bam_file.name
    """
    # if this is a cram file, convert to bam first
    if [[ ${bam_file} == *.cram ]]; then
        samtools view -@ ${task.cpus} -b -o ${input_bam} ${bam_file}
    else
        ln -s ${bam_file} ${input_bam}
    fi
    
    ## Mapdamage command - use fixed output directory name for reproducibility
    mapDamage -i ${input_bam} -r ${reference_genome} -d results_${sample_id} ${rescale_arg}

    # sort the rescaled bam file and convert to cram
    if [[ -f results_${sample_id}/*rescaled.bam ]]; then
        # if there is a rescaled bam, sort, convert to cram, and this will be the downstream bam
        samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam results_${sample_id}/*rescaled.bam
        samtools view -@ ${task.cpus} -C -T ${reference_genome} -o ${sample_id}.dedup.cram ${sample_id}.sorted.bam
        samtools index ${sample_id}.dedup.cram
        rm -f ${sample_id}.sorted.bam
    else
        # if there is no rescaled bam, just convert the original bam back to cram and use that as the downstream bam
        samtools view -@ ${task.cpus} -C -T ${reference_genome} -o ${sample_id}.dedup.cram ${input_bam}
        samtools index ${sample_id}.dedup.cram
    fi

    # move damage profiler results to a subdirectory with the sample name
    mkdir -p ${sample_id}_damage
    mv results_${sample_id}/*.pdf ${sample_id}_damage/ 2>/dev/null || true

    """

    stub:
    """
    touch ${sample_id}.dedup.cram
    
    """
}