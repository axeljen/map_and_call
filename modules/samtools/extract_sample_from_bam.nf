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
    samtools view -H ${bam} | awk -F '\\t' '\$1 == "@RG" { for (i = 1; i <= NF; i++) if (\$i ~ /^SM:/) { split(\$i, tag, ":"); print tag[2]; exit } }' > ${bam.baseName}.sample_id.txt
    """

    stub:
    """
    echo -n "SAMPLE_${bam.baseName}"
    """
}

