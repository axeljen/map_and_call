process adapterremoval {
    tag "adapterremoval"
    label "wide"
    conda "${moduleDir}/environment.yml"

    //publishDir "results/adapterremoval", mode: 'copy'

    input:
    tuple val(sample_id), val(lane), val(datatype), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), val(lane), path("${sample_id}_${lane}_trimmed_R1.fastq.gz"), path("${sample_id}_${lane}_trimmed_R2.fastq.gz"), path("${sample_id}_${lane}_collapsed.fastq.gz"), path("${sample_id}_${lane}_singletons.fastq.gz"), emit: trimmed_reads

    script:
    """
    # Adapter trimming and read merging with adapterremoval
    AdapterRemoval --file1 ${reads_1} --file2 ${reads_2} \
        --threads 8 \
        --basename ${sample_id}_${lane} \
        --trimns --trimqualities --collapse \
        --output1 ${sample_id}_${lane}_trimmed_R1.fastq.gz \
        --output2 ${sample_id}_${lane}_trimmed_R2.fastq.gz \
        --outputcollapsed ${sample_id}_${lane}_collapsed.fastq.gz \
        --singleton ${sample_id}_${lane}_singletons.fastq.gz \
        --gzip


    """
    stub:
    """
    touch ${sample_id}_${lane}_trimmed_R1.fastq.gz
    touch ${sample_id}_${lane}_trimmed_R2.fastq.gz
    touch ${sample_id}_${lane}_collapsed.fastq.gz
    touch ${sample_id}_${lane}_singletons.fastq.gz
    """

}
