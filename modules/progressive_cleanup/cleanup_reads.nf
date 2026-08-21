// Check 

process cleanup_reads {
    tag "$sample_id"
    label 'thin_short'

    input:
    tuple val(sample_id), val(library), val(datatype), val(reads)

    output:
    tuple val(sample_id), val(library), val(datatype), path("${sample_id}_${library}_cleaned_reads.txt"), emit: deleted_files

    script:
    def readlist = reads.collect { read -> "'${read.toString()}'" }.join(' ')

    """
    # get work root for safety check, not to remove anything outside of the work dir
    work_root=\$(realpath ${workflow.workDir})

    # delete intermediate read files once the final, clean reads have been generated, to save storage space
    for read in ${readlist};
        do
        target=\$(readlink -f \${read})
        # check that the target exists, otherwise error
        if [[ -z "\${target}" || ! -f "\${target}" ]];
        then
            echo "Error: target read file does not exist: \${target}"
            exit 1
        fi
        # check if the target is within the work root, otherwise error
        if [[ "\${target}" != \${work_root}/* ]];
        then
            echo "Error: target read file is outside of the work root: \${target}"
            exit 1
        fi
        # otherwise delete target, but for now simply echo the target to be deleted, for safety
        echo "Deleting intermediate read file: \${target}" >> ${sample_id}_${library}_cleaned_reads.txt
        rm -f "\${target}"
    done

    """

    stub:
    """
    touch ${sample_id}_${library}_cleaned_reads.txt
    """
}
