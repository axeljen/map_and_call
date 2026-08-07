// ═══════════════════════════════════════════════════════════════════════════════
//                   SHARED WORKFLOW UTILITY FUNCTIONS
//
// Common helper functions used across multiple workflow files
// ═══════════════════════════════════════════════════════════════════════════════

class WorkflowUtils {
    
    /*
     * Setup sex chromosome system and identify sex-linked contigs
     * Returns map with: sex_chrom_system, sex_linked_list, sex_limited_list, non_sex_limited_list
     */
    static def setupSexChromosomeSystem(params) {
        def result = [
            sex_chrom_system: 'unknown',
            sex_linked_list: [],
            sex_limited_list: [],
            non_sex_limited_list: []
        ]
        
        if (params.x_scaffolds || params.y_scaffolds) {
            if (params.z_scaffolds || params.w_scaffolds) {
                throw new Exception("Please only specify either X and/or Y, OR Z and/or W scaffolds, not both.")
            }
            result.sex_chrom_system = 'XY'
            result.sex_linked_list = [params.x_scaffolds, params.y_scaffolds].flatten()
            result.sex_limited_list = [params.y_scaffolds].flatten()
            result.non_sex_limited_list = [params.x_scaffolds].flatten()
        }
        else if (params.z_scaffolds || params.w_scaffolds) {
            result.sex_chrom_system = 'ZW'
            result.sex_linked_list = [params.z_scaffolds, params.w_scaffolds].flatten()
            result.sex_limited_list = [params.w_scaffolds].flatten()
            result.non_sex_limited_list = [params.z_scaffolds].flatten()
        }
        
        return result
    }
    
    /*
     * Calculate sample depth statistics and sex assignments
     * Inputs:
     *   - depth_avg_ch: channel from parse_region_depths.sample_depth_avg
     *   - sex_limited_list: list of sex-limited scaffolds (Y/W)
     *   - non_sex_limited_list: list of non-sex-limited sex chromosomes (X/Z)
     *   - lower_threshold: coverage ratio lower threshold
     *   - upper_threshold: coverage ratio upper threshold
     * Output:
     *   channel with tuple(sample, autosome_depth, non_sex_limited_depth, sex_limited_depth, ratio, sex)
     */
    static def calculateDepthAndSex(depth_avg_ch, sex_limited_list, non_sex_limited_list, lower_threshold, upper_threshold) {
        return depth_avg_ch
            .flatMap { sample, file ->
                file.splitCsv(header: false, sep: '\t').collect { row ->
                    def chrom = row[0]
                    def chrom_length = row[1].toDouble()
                    def avg_depth = row[2].toDouble()
                    return tuple(sample, chrom, chrom_length, avg_depth)
                }
            }
            // Categorize chromosomes
            .map { sample, chrom, chrom_length, avg_depth ->
                def category = 'autosomes'
                if (sex_limited_list.contains(chrom)) {
                    category = 'sex_limited'
                }
                if (non_sex_limited_list.contains(chrom)) {
                    category = 'non_sex_limited'
                }
                return tuple(sample, chrom, chrom_length, avg_depth, category)
            }
            // Group by sample and category, calculate weighted average depth per category
            .groupTuple(by: [0, 4])
            .map { sample, chroms, lengths, values, category ->
                // Weighted average: sum(depth * length) / sum(length)
                def total_bases = [lengths, values].transpose().collect { len, depth -> len * depth }.sum()
                def total_length = lengths.sum()
                def avg_depth = total_length > 0 ? total_bases / total_length : 0.0
                return tuple(sample, [(category): avg_depth])
            }
            .groupTuple(by: 0)
            .map { sample, category_maps ->
                // Merge all category maps into a single map per sample
                def depth_map = [:]
                category_maps.each { map ->
                    depth_map.putAll(map)
                }
                tuple(sample, depth_map)
            }
            // Calculate ratio and assign sex
            .map { sample, depth_map ->
                def autosome_depth = depth_map.get('autosomes', 0.0)
                def non_sex_limited_depth = depth_map.get('non_sex_limited', 0.0)
                def sex_limited_depth = depth_map.get('sex_limited', 0.0)
                
                def ratio = (autosome_depth > 0 && non_sex_limited_depth > 0) ? 
                            non_sex_limited_depth / autosome_depth : 0.0
                
                def sex = 'unknown'
                if (ratio > 0) {
                    if (ratio < lower_threshold) {
                        sex = 'unknown'
                    } else if (ratio >= 1.5) {
                        sex = 'unknown'
                    } else if (lower_threshold <= ratio && ratio <= upper_threshold) {
                        sex = 'hemizygous'
                    } else {
                        sex = 'homozygous'
                    }
                }
                
                tuple(sample, autosome_depth, non_sex_limited_depth, sex_limited_depth, ratio, sex)
            }
    }
}
