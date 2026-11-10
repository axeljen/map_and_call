process sratools_fastq_dump {
    tag "sra_fastqdump"
    label 'medium_long'
    // conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), val(lane), val(datatype), val(library), val(sra_id)
    

    output:
    tuple val(sample_id), val(lane), val(datatype), val(library), path("${sra_id}_1.f*q.gz"), path("${sra_id}_2.f*q.gz"), emit: reads


    script:
    """
    source task_tmp.sh ${params.use_scratch} "\${TASK_SCRATCH_ROOT:-${params.use_scratch ? params.scratch_dir ?: '' : ''}}"
    df -hT "\$TMPDIR"

    prefetch $sra_id --max-size u -O "\$TMPDIR"
    echo "Starting fasterq-dump"
    fasterq-dump "\$TMPDIR/$sra_id" -t "\$TMPDIR" -O "\$TMPDIR" -e $task.cpus
    echo "Starting compression"
    pigz -p $task.cpus -c "\$TMPDIR/${sra_id}_1.fastq" > ${sra_id}_1.fastq.gz
    rm "\$TMPDIR/${sra_id}_1.fastq"
    pigz -p $task.cpus -c "\$TMPDIR/${sra_id}_2.fastq" > ${sra_id}_2.fastq.gz

    """

    stub:
    """
    touch ${sra_id}_1.fastq.gz ${sra_id}_2.fastq.gz
    """
}