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
    tuple val(sample_id), path("${sample_id}.cram"), path("${sample_id}.cram.crai"), emit: rescaled_bam

    script:
    def rescale_arg = params.damageprofiler_rescale ? '--rescale' : ''
    """
    bam_file=${bam_file}
    # if this is a cram file, convert to bam first
    if [[ ${bam_file} == *.cram ]]; then
        samtools view -@ ${task.cpus} -b -o tmp.${sample_id}.bam ${bam_file}
        bam_file=tmp.${sample_id}.bam
    fi
    
    
    ## Mapdamage command
    mapDamage -i \$bam_file -r ${reference_genome} ${rescale_arg}

    # sort the rescaled bam file and convert to cram
    rescaled_bam=\$(find . -name "*rescaled.bam" | head -n 1)
    # if there is a rescaled bam, sort, convert to cram, and this will be the downstream bam
    if [[ -f \$rescaled_bam ]]; then
    samtools sort -@ ${task.cpus} -o tmp.${sample_id}.bam \${rescaled_bam}
    mv tmp.${sample_id}.bam ${sample_id}.bam
    samtools view -@ ${task.cpus} -C -o ${sample_id}.cram ${sample_id}.bam
    samtools index ${sample_id}.cram
    # remove the intermediate bam file
    rm -f ${sample_id}.bam
    else
    # if there is no rescaled bam, just convert the original bam back to cram and use that as the downstream bam
    samtools view -@ ${task.cpus} -C -o ${sample_id}.cram \$bam_file
    samtools index ${sample_id}.cram
    fi

    # move damage profiler results to a subdirectory with the sample name
    mkdir -p ${sample_id}_damage
    mv results_*${sample_id}*/*.pdf ${sample_id}_damage/

    """

    stub:
    """
    touch ${sample_id}.bam
    
    """
}