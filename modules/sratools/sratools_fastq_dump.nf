process sratools_fastq_dump {
    scratch (params.use_scratch ? params.scratch_dir : false)
    tag "sra_fastqdump"
    label 'medium_long'
    // conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), val(lane), val(datatype), val(library), val(sra_id)
    

    output:
    tuple val(sample_id), val(lane), val(datatype), val(library), path("${sra_id}_1.f*q.gz"), path("${sra_id}_2.f*q.gz"), emit: reads


    script:
    """
    echo "Task directory:"
    pwd
    df -hT .

    prefetch $sra_id --max-size u -O .
    echo "Starting fasterq-dump"
    fasterq-dump ./$sra_id -t . -O . -e $task.cpus
    echo "Starting compression"
    pigz -p $task.cpus ${sra_id}_1.fastq ${sra_id}_2.fastq

    """

    stub:
    """
    touch ${sra_id}_1.fastq.gz ${sra_id}_2.fastq.gz
    """
}