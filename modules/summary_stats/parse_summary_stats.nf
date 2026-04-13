process parse_summary_stats {
    tag "parse_summary_stats"
    label 'process_small'
    // conda "${moduleDir}/environment.yml"


    input: 
    tuple val(sample_id), val(pop), path(invariants), path(callable), path(snpmask), val(raw_bamdepth_autosomes), val(raw_bamdepth_non_sex_limited), val(raw_bamdepth_sex_limited), val(ratio_non_sex_limited_to_autosomes), val(sex_assignment), val(ds_bamdepth_autosomes), val(ds_bamdepth_non_sex_limited), val(ds_bamdepth_sex_limited), path(combined_stats_snps), path(combined_stats_indels)
    val(sex_chromosome_system)
    output:
    path("${sample_id}_summary_statistics.txt"), emit: summary_statistics

    script:
    """
    
    # total number of callable sites
    callable_sites=\$(awk '{sum += \$3 - \$2} END {print sum}' ${callable})
    # fetch stats from snpstats
    read -r num_records num_homref num_het num_homalt num_missing < <(
        awk -v sample="${sample_id}" '\$1==sample {print \$2, \$3, \$4, \$5, \$6; exit}' *_snps_sample_stats.tsv
        )
    # number of homozygous reference calls outside of the vcf file
    homref_invariants=\$(awk '{sum += \$3 - \$2} END {print sum}' ${invariants})

    # total number of homozygous calls
    total_hom=\$((num_homref + homref_invariants + num_homalt))

    # heterozygosity as snps / snps + total_hom
    heterozygosity=\$(echo "scale=6; \$num_het / (\$num_het + \$total_hom)" | bc)

    # pending the sex chromosome system and sex assignment, tidy the labels to male/female
    if [[ "${sex_assignment}" == "hemizygous" && ${sex_chromosome_system} == "XY" ]]; then
        parsed_sex_assignment="male"
    elif [[ "${sex_assignment}" == "homozygous" && ${sex_chromosome_system} == "XY" ]]; then
        parsed_sex_assignment="female"
    elif [[ "${sex_assignment}" == "hemizygous" && ${sex_chromosome_system} == "ZW" ]]; then
        parsed_sex_assignment="female"
    elif [[ "${sex_assignment}" == "homozygous" && ${sex_chromosome_system} == "ZW" ]]; then
        parsed_sex_assignment="male"
    else
        parsed_sex_assignment="${sex_assignment}"
    fi

    # and better names for sex_linked vs non_sex_linked
    if [[ "${sex_chromosome_system}" == "XY" ]]; then
        parsed_sex_limited_label="Y"
        parsed_non_sex_limited_label="X"
    elif [[ "${sex_chromosome_system}" == "ZW" ]]; then
        parsed_sex_limited_label="W"
        parsed_non_sex_limited_label="Z"
    else
        parsed_sex_limited_label="sex_limited"
        parsed_non_sex_limited_label="non_sex_limited"
    fi

    # print a table with the following columns: sample_id, population, sex_assignment, bam_depth_autosomes, bam_depth_non_sex_limited, bam_depth_sex_limited, ratio_non_sex_limited_to_autosomes, downsampled_bam_depth_autosomes, downsampled_bam_depth_non_sex_limited, downsampled_bam_depth_sex_limited, callable_sites, num_snps, num_missing_snps, heterozygosity
    echo -e "sample_id\tpopulation\tsex_assignment\tbam_depth_autosomes\tbam_depth_\${parsed_non_sex_limited_label}\tbam_depth_\${parsed_sex_limited_label}\tratio_\${parsed_non_sex_limited_label}_to_autosomes\tdownsampled_bam_depth_autosomes\tdownsampled_bam_depth_\${parsed_non_sex_limited_label}\tdownsampled_bam_depth_\${parsed_sex_limited_label}\tcallable_sites\tnum_snps\tnum_missing_snps\theterozygosity" > ${sample_id}_summary_statistics.txt
    echo -e "${sample_id}\t${pop}\t\${parsed_sex_assignment}\t${raw_bamdepth_autosomes}\t${raw_bamdepth_non_sex_limited}\t${raw_bamdepth_sex_limited}\t${ratio_non_sex_limited_to_autosomes}\t${ds_bamdepth_autosomes}\t${ds_bamdepth_non_sex_limited}\t${ds_bamdepth_sex_limited}\t\${callable_sites}\t\${num_records}\t\${num_missing}\t\${heterozygosity}" >> ${sample_id}_summary_statistics.txt

    """

    stub:
    """
    echo -e "" > ${sample_id}_summary_statistics.txt
    """

}
