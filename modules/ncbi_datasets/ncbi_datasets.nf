
process datasets_download {
    scratch (params.use_scratch ? params.scratch_dir : false)
    tag "datasets_download"
    label 'medium_long'
    // conda "${moduleDir}/environment.yml"

    input:
    val(accession_id)
    

    output:
    path("${accession_id}_genomic.fna.gz"), emit: reference_genome


    script:
    """
    # Set SSL certificate file for ncbi-datasets-cli
    export SSL_CERT_FILE=/opt/conda/ssl/cacert.pem
    datasets download genome accession "${accession_id}" --filename "${accession_id}.zip"

    unzip "${accession_id}.zip" -d "${accession_id}"

    # find the genomic FASTA file and move here
    genome=\$(find "${accession_id}" -type f -name "*_genomic.f*a")
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

    rm -rf "${accession_id}" "${accession_id}.zip"

    """

    stub:
    """
    touch "${accession_id}_genomic.fna.gz"
    """
}