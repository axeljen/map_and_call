// Merge bam alignments generated from different sequencing runs/lanes for the same sample
process samtools_merge {
    tag "$sample_id"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), path(bam_files)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: merged_bam

    script:
    """
    samtools merge -@ ${task.cpus} ${sample_id}.bam ${bam_files.join(' ')}
    # echo "Simulated merged BAM content for sample ${sample_id}" > ${sample_id}.bam
    """

    stub:
    """
    touch ${sample_id}.bam
    """
}

// Merge bam alignments generated from different sequencing runs/lanes for the same sample

process merge_historical_bams {
    tag "$sample_id"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), val(lane), path(bam_files)

    output:
    tuple val(sample_id), val(lane), path("${sample_id}_${lane}.bam"), path("${sample_id}_${lane}.bam.bai"), emit: merged_bam

    script:
    """
    samtools merge -@ ${task.cpus} ${sample_id}_${lane}.bam ${bam_files.join(' ')}
    samtools index ${sample_id}_${lane}.bam
    # echo "Simulated merged BAM content for sample ${sample_id}" > ${sample_id}_${lane}.bam
    """

    stub:
    """
    touch ${sample_id}.bam
    touch ${sample_id}.bam.bai
    """
}