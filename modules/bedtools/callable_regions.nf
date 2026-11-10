process callable_regions {
    tag "callable_regions"
    label 'medium_short'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/03_sample_callability", mode: 'copy'

    input: 
    tuple val(sample_id), val(min_dp), val(max_dp), val(sex_assignment), path(bedfile), val(sex_limited_scaffolds), val(non_sex_limited_scaffolds), path(faidx)
    
    output:
    tuple val(sample_id), path("${sample_id}.callable_regions.bed.gz"), emit: callable

    script:
    def sex_limited_scaffolds_list = sex_limited_scaffolds.join(' ')
    def non_sex_limited_scaffolds_list = non_sex_limited_scaffolds.join(' ')
    // if min dp is an integer, this should be used also on sex-linked contigs, otherwise half it if the sex_assignment is hemizygous
    def min_dp_sexlinked = sex_assignment == "hemizygous" ? (min_dp instanceof Integer ? min_dp / 2 : min_dp) : min_dp
    def max_dp_sexlinked = sex_assignment == "hemizygous" ? (max_dp instanceof Integer ? max_dp / 2 : max_dp) : max_dp

    """
    source task_tmp.sh ${params.use_scratch} "\${TASK_SCRATCH_ROOT:-${params.use_scratch ? params.scratch_dir ?: '' : ''}}"
    # make a temporary bedfile while filtering
    zcat ${bedfile} > "\$TMPDIR/callable.bed"

    # make sex-limited and non-sexlimited bedstrings based on the input lists
    for scaffold in $sex_limited_scaffolds_list; do
        start=0
        end=\$(awk -v scaffold="\$scaffold" '\$1 == scaffold {print \$2}' ${faidx})
        echo -e "\${scaffold}\\t\${start}\\t\${end}"
        # remove this scaffold from the temporary bedfile
        awk -v scaffold="\$scaffold" '\$1 != scaffold' "\$TMPDIR/callable.bed" > "\$TMPDIR/callable.tmp.bed"
        mv "\$TMPDIR/callable.tmp.bed" "\$TMPDIR/callable.bed"
    done > "\$TMPDIR/sex_limited_scaffolds.bed"

    for scaffold in $non_sex_limited_scaffolds_list; do
        start=0
        end=\$(awk -v scaffold="\$scaffold" '\$1 == scaffold {print \$2}' ${faidx})
        echo -e "\${scaffold}\\t\${start}\\t\${end}"
        # remove this scaffold from the temporary bedfile
        awk -v scaffold="\$scaffold" '\$1 != scaffold' "\$TMPDIR/callable.bed" > "\$TMPDIR/callable.tmp.bed"
        mv "\$TMPDIR/callable.tmp.bed" "\$TMPDIR/callable.bed"
    done > "\$TMPDIR/non_sex_limited_scaffolds.bed"

    # filter the temporary bedfile based on min and max depth
    awk -v min_dp="${min_dp}" -v max_dp="${max_dp}" -v OFS='\\t' ' \$4 >= min_dp && \$4 <= max_dp {print \$1, \$2, \$3}' "\$TMPDIR/callable.bed" > "\$TMPDIR/callable.filtered.bed"
    mv "\$TMPDIR/callable.filtered.bed" "\$TMPDIR/callable.bed"

    # if sex assignment is hemizygous, half the min and max dp for both the sex-limited and non-sex-limited scaffolds
    if [[ "${sex_assignment}" == "hemizygous" ]]; then
        # and filter both these files into the final callable bedfile
        cat "\$TMPDIR/sex_limited_scaffolds.bed" "\$TMPDIR/non_sex_limited_scaffolds.bed" | awk -v min_dp=${min_dp_sexlinked} -v max_dp=${max_dp_sexlinked} -v OFS='\\t' ' \$4 >= min_dp && \$4 <= max_dp {print \$1, \$2, \$3}' >> "\$TMPDIR/callable.bed"
    # otherwise if homozygous, keep same min and max dp and simply skip the sex-limited scaffolds (if any)
    else
        cat "\$TMPDIR/non_sex_limited_scaffolds.bed" | awk -v min_dp="${min_dp}" -v max_dp="${max_dp}" -v OFS='\\t' ' \$4 >= min_dp && \$4 <= max_dp {print \$1, \$2, \$3}' >> "\$TMPDIR/callable.bed"
    fi

    # finally, mask out any regions that are in the reference mask bedfile
    # if a reference mask was provided as a parameter, use that, otherwise just create an empty bedfile to use as the reference mask
    if [[ -f "${params.reference_mask}" ]]; then
        reference_mask="${params.reference_mask}"
    else
        touch "\$TMPDIR/reference_mask.bed"
        reference_mask="\$TMPDIR/reference_mask.bed"
    fi
    bedtools subtract -a "\$TMPDIR/callable.bed" -b "\${reference_mask}" > "\$TMPDIR/callable_regions.bed"

    # Sort and merge overlapping/adjacent regions
    # First merge per-scaffold using GNU sort (disk-based) to avoid OOM with large files,
    # then re-sort by faidx scaffold order on the much smaller merged output
    bedtools merge -i "\$TMPDIR/callable_regions.bed" > "\$TMPDIR/merged.bed"
    # and sort the final output so that the sex-limited scaffolds show up where they should
    bedtools sort -faidx ${faidx} -i "\$TMPDIR/merged.bed" | gzip -c > ${sample_id}.callable_regions.bed.gz

    """

    stub:
    """
    touch ${sample_id}.callable_regions.bed.gz
    """

}
