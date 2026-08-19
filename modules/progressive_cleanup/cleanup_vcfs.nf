// Check 

process cleanup_vcfs {
    tag "$vcf"
    label 'thin_short'

    input:
    val vcf

    output:
    path "cleaned_vcfs.txt", emit: deleted_files

    script:

    """
    # get work root for safety check, not to remove anything outside of the work dir
    work_root=\$(realpath ${workflow.workDir})

    # delete intermediate BAM files once the final, clean BAMs have been generated, to save storage space
    target=\$(readlink -f ${vcf})
    # check that the target exists, otherwise error
    if [[ -z "\${target}" || ! -f "\${target}" ]];
    then
        echo "Error: target VCF file does not exist: \${target}"
        exit 1
    fi
    # check if the target is within the work root, otherwise error
    if [[ "\${target}" != \${work_root}/* ]];
    then
        echo "Error: target VCF file is outside of the work root: \${target}"
        exit 1
    fi
    # otherwise delete target, but for now simply echo the target to be deleted, for safety
    echo "Deleting intermediate VCF file: \${target}" >> cleaned_vcfs.txt
    rm -f "\${target}"

    """

    stub:
    """
    touch cleaned_vcfs.txt
    """
}
