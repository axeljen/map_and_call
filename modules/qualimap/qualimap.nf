// Check 

process qualimap {
    tag "$sample_id"
    label 'medium_medium'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), path(bam_file), path(index_file)

    output:
    tuple val(sample_id), path("${sample_id}_qualimap.pdf"), emit: qualimap_report

    script:
    def mem = (task.memory.toMega() * 0.8) as long
    """
    bamfile=${bam_file}
    # convert from cram to bam if in cram format
    if [[ ${bam_file} == *.cram ]]; then
        source task_tmp.sh ${params.use_scratch} "\${TASK_SCRATCH_ROOT:-${params.use_scratch ? params.scratch_dir ?: '' : ''}}"
        samtools view -@ ${task.cpus} -b -o "\$TMPDIR/${sample_id}.bam" ${bam_file}
        bamfile="\$TMPDIR/${sample_id}.bam"
    fi
    qualimap bamqc --java-mem-size=${mem}M -bam "\$bamfile" -outfile ${sample_id}_qualimap.pdf -outdir . -outformat PDF
    # move the final pdf here
    pdf_file=\$(find . -name "*qualimap.pdf" -type f -print -quit)
    if [[ "\$pdf_file" != "./${sample_id}_qualimap.pdf" ]]; then
        mv "\$pdf_file" ${sample_id}_qualimap.pdf
    fi


    """

    stub:
    """
    touch ${sample_id}_qualimap.pdf
    """
}