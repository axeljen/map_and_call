process extract_sample_from_bam {
    scratch params.use_scratch
    tag "$bam"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    tuple path(bam), path(bai)

    output:
    tuple path(bam), path(bai), path("${bam.baseName}.sample_id.txt"), emit: sample_bams

    script:
    """
    # Extract sample ID from BAM read group (SM tag)
    samtools view -H ${bam} | grep '^@RG' | head -n 1 | sed 's/.*SM:\\([^\\t]*\\).*/\\1/' | tr -d '\\n' > ${bam.baseName}.sample_id.txt
    """

    stub:
    """
    echo -n "SAMPLE_${bam.baseName}"
    """
}

