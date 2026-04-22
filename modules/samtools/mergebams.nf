// Merge bam alignments generated from different sequencing runs/lanes for the same sample
process samtools_merge {
    tag "$sample_id"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), val(datatype), path(bam_files)

    output:
    tuple val(sample_id), val(datatype), path("${sample_id}.bam"), emit: merged_bam

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
    tuple val(sample_id), val(library), val(datatype), path(bam_files)

    output:
    tuple val(sample_id), val(library), val(datatype), path("${sample_id}_${library}.bam"), path("${sample_id}_${library}.bam.bai"), emit: merged_bam

    script:
    """
    samtools merge -@ ${task.cpus} ${sample_id}_${library}.bam ${bam_files.join(' ')}
    samtools index ${sample_id}_${library}.bam
    # echo "Simulated merged BAM content for sample ${sample_id}" > ${sample_id}_${library}.bam
    """

    stub:
    """
    touch ${sample_id}_${library}.bam
    touch ${sample_id}_${library}.bam.bai
    """
}