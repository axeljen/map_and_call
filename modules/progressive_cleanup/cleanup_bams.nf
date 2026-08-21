// Check 

process cleanup_bams {
    tag "$sample_id"
    label 'thin_short'

    input:
    tuple val(sample_id), val(datatype), val(bam), val(bai)

    output:
    tuple val(sample_id), val(datatype),  path("${sample_id}_cleaned_bams.txt"), emit: deleted_files

    script:

    """
    # get work root for safety check, not to remove anything outside of the work dir
    work_root=\$(realpath ${workflow.workDir})

    # delete intermediate BAM files once the final, clean BAMs have been generated, to save storage space
    target=\$(readlink -f ${bam})
    # check that the target exists, otherwise error
    if [[ -z "\${target}" || ! -f "\${target}" ]];
    then
        echo "Error: target BAM file does not exist: \${target}"
        exit 1
    fi
    # check if the target is within the work root, otherwise error
    if [[ "\${target}" != \${work_root}/* ]];
    then
        echo "Error: target BAM file is outside of the work root: \${target}"
        exit 1
    fi
    # otherwise delete target, but for now simply echo the target to be deleted, for safety
    echo "Deleting intermediate BAM file: \${target}" >> ${sample_id}_cleaned_bams.txt
    rm -f "\${target}"

    """

    stub:
    """
    touch ${sample_id}_cleaned_bams.txt
    """
}
