
process datasets_download {
    tag "datasets_download"
    label 'medium_long'
    // conda "${moduleDir}/environment.yml"

    input:
    val(accession_id)
    

    output:
    path("${accession_id}_genomic.fna.gz"), emit: reference_genome


    script:
    """
    source task_tmp.sh ${params.use_scratch} "\${TASK_SCRATCH_ROOT:-${params.use_scratch ? params.scratch_dir ?: '' : ''}}"
    # Set SSL certificate file for ncbi-datasets-cli
    export SSL_CERT_FILE=/opt/conda/ssl/cacert.pem
    datasets download genome accession "${accession_id}" --filename "\$TMPDIR/${accession_id}.zip"

    unzip "\$TMPDIR/${accession_id}.zip" -d "\$TMPDIR/${accession_id}"

    # find the genomic FASTA file and move here
    genome=\$(find "\$TMPDIR/${accession_id}" -type f -name "*_genomic.f*a")
    if [[ -z "\$genome" ]]; then
        echo "Genomic FASTA file not found"
        exit 1
    fi
    # if more than one file
    if [[ \$(echo "\$genome" | wc -l) -gt 1 ]]; then
        echo "Multiple genomic FASTA files found"
        exit 1
    fi
    # bgzip the reference
    bgzip -c "\$genome" > "${accession_id}_genomic.fna.gz"

    """

    stub:
    """
    touch "${accession_id}_genomic.fna.gz"
    """
}