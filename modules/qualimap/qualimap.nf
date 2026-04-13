// Check 

process qualimap {
    tag "$sample_id"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), path(bam_file), path(index_file)

    output:
    tuple val(sample_id), path("${sample_id}_qualimap.pdf"), emit: qualimap_report

    script:
    """
    bamfile=${bam_file}
    # convert from cram to bam if in cram format
    if [[ ${bam_file} == *.cram ]]; then
        samtools view -@ ${task.cpus} -b -o tmp.${sample_id}.bam ${bam_file}
        bamfile=tmp.${sample_id}.bam
    fi
    qualimap bamqc -bam \$bamfile -outfile ${sample_id}_qualimap.pdf -outdir . -outformat PDF
    # move the final pdf here
    pdf_file=\$(find . -name "*qualimap.pdf" | head -n 1)
    mv \$pdf_file ${sample_id}_qualimap.pdf

    # remove the intermediate bam file if it was created
    if [[ ${bam_file} == *.cram ]]; then
        rm -f tmp.${sample_id}.bam
    fi


    """

    stub:
    """
    touch ${sample_id}_qualimap.pdf
    """
}