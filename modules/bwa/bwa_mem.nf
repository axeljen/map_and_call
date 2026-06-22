/*
 * BWA MEM - Read alignment with read group tagging
 */

process bwa_mem {
    tag "$sample_id"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), val(library), val(datatype), path(read1), path(read2), path(reference), path(reference_index)
    

    output:
    tuple val(sample_id), val(library), val(datatype), path("${sample_id}_${library}.sorted.bam"), path("${sample_id}_${library}.sorted.bam.bai"), emit: bam

    script:
    """
    ## extract the read group string from the fasta header
    # HEADER=\$(zcat ${read1} | head -1 | sed 's/^@//')
    
    ## detect platform/read group info format
    ## Detect format by structure
    #if echo "\$HEADER" | grep -qP '^\\S+:\\d+:[A-Z0-9]+:\\d+:'; then
    #    # Illumina Casava 1.8+ : instrument:run:flowcell:lane:tile:x:y
    #    PLATFORM="ILLUMINA"
    #    FLOWCELL=\$(echo "\$HEADER" | cut -d':' -f3)
    #    LANE=\$(echo "\$HEADER" | cut -d':' -f4)
    #elif echo "\$HEADER" | grep -qP 'L\\d+C\\d+R\\d+'; then
    #    # BGI/MGI: {flowcell}L{lane}C{col}R{row}{read}
    #    PLATFORM="BGI"
    #    FLOWCELL=\$(echo "\$HEADER" | sed 's/L[0-9]\\+.*//' | sed 's|/.*||')
    #    LANE=\$(echo "\$HEADER" | grep -oP '(?<=L)\\d+' | head -1)
    #else
    #    # Old Illumina or unknown: instrument:lane:tile:x:y
    #    PLATFORM="ILLUMINA"
    #    FLOWCELL=\$(echo "\$HEADER" | cut -d':' -f1)
    #    LANE=\$(echo "\$HEADER" | cut -d':' -f2)
    #fi

    FLOWCELL="fc1"
    LANE="${sample_id}_${library}"
    RGID="${sample_id}_${library}"
    PU="illumina"
    SAMPLE="${sample_id}"

    rg="@RG\\tID:\${RGID}\\tSM:\${SAMPLE}\\tLB:\${SAMPLE}\\tPL:\${PLATFORM}\\tPU:\${PU}"

    bwa mem \
         -t ${task.cpus} \
         -R "\${rg}" \
         ${reference} \
         ${read1} \
         ${read2} \
     | samtools view -@ ${task.cpus} -q ${params.min_mapqual} -b -o - - | \
     samtools sort -@ ${task.cpus} -o ${sample_id}_${library}.sorted.bam -

    samtools index ${sample_id}_${library}.sorted.bam
    """

    stub:
    """
    touch ${sample_id}_${library}.sorted.bam
    touch ${sample_id}_${library}.sorted.bam.bai
    """
}


process bwa_mem_singlereads {
    tag "$sample_id"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), val(library), val(datatype), path(reads), path(reference), path(reference_index)

    output:
    tuple val(sample_id), val(library), val(datatype), path("${sample_id}_${library}_MR.sorted.bam"), path("${sample_id}_${library}_MR.sorted.bam.bai"), emit: bam

    script:
    """
    # extract the read group string from the fasta header
    HEADER=\$(zcat ${reads} | head -1 | sed 's/^@//')
    
    ## detect platform/read group info format
    ## Detect format by structure
    #if echo "\$HEADER" | grep -qP '^\\S+:\\d+:[A-Z0-9]+:\\d+:'; then
    #    # Illumina Casava 1.8+ : instrument:run:flowcell:lane:tile:x:y
    #    PLATFORM="ILLUMINA"
    #    FLOWCELL=\$(echo "\$HEADER" | cut -d':' -f3)
    #    LANE=\$(echo "\$HEADER" | cut -d':' -f4)
    #elif echo "\$HEADER" | grep -qP 'L\\d+C\\d+R\\d+'; then
    #    # BGI/MGI: {flowcell}L{lane}C{col}R{row}{read}
    #    PLATFORM="BGI"
    #    FLOWCELL=\$(echo "\$HEADER" | sed 's/L[0-9]\\+.*//' | sed 's|/.*||')
    #    LANE=\$(echo "\$HEADER" | grep -oP '(?<=L)\\d+' | head -1)
    #else
    #    # Old Illumina or unknown: instrument:lane:tile:x:y
    #    PLATFORM="ILLUMINA"
    #    FLOWCELL=\$(echo "\$HEADER" | cut -d':' -f1)
    #    LANE=\$(echo "\$HEADER" | cut -d':' -f2)
    #fi

    FLOWCELL="fc1"
    LANE="${sample_id}_${library}"
    RGID="${sample_id}_${library}"
    PU="illumina"
    SAMPLE="${sample_id}"

    rg="@RG\\tID:\${RGID}\\tSM:\${SAMPLE}\\tLB:${library}\\tPL:\${PLATFORM}\\tPU:\${PU}"

    bwa mem \
         -t ${task.cpus} \
         -R "\${rg}" \
         ${reference} \
         ${reads} \
     | samtools view -@ ${task.cpus} -q ${params.min_mapqual} -b -o - - | \
     samtools sort -@ ${task.cpus} -o ${sample_id}_${library}_MR.sorted.bam -

    samtools index ${sample_id}_${library}_MR.sorted.bam
    """

    stub:
    """
    touch ${sample_id}_${library}_MR.sorted.bam
    touch ${sample_id}_${library}_MR.sorted.bam.bai
    """
}
