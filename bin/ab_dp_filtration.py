#!/usr/bin/env python3

import argparse
import pysam

args = argparse.ArgumentParser(description="Calculate allele balance for each variant in a VCF file and filter variants based on allele balance.")

args.add_argument("-i", "--input", help="Input VCF file")
args.add_argument("-o", "--output", help="Output VCF file")
args.add_argument("--min-ab", type=float, help="Minimum allele balance threshold")
args.add_argument('-s', '--samples', help='comma-separated list of samples to include (default: all)', default=None)
args.add_argument('--min-depth', help='Minimum depth to keep a genotype', type=float, default=0)
args.add_argument('--max-depth', help='Maximum depth to keep a genotype', type=float, default=200)
args.add_argument('--sex-assignments', help='Comma separated list of sex assignments as hemizygous or homozygous of included samples. IMPORTANT: If multiple samples are included, the order must correspond exactly to the order of samples specified in --samples or, if this argument is omitted, the order of samples in the VCF header. Unless any sex-linked scaffolds are provided, this argument has no effect.')
args.add_argument('--sex-linked-scaffolds', help='Comma separated list of sex-linked scaffolds. If this is provided in combination with sex assignments, the coverage thresholds will be changed to half min-depth and max-depth for samples assigned as hemizygous on these scaffolds.')

args = args.parse_args()

# Open the input VCF file
vcf_in = pysam.VariantFile(args.input)
# Create the output VCF file
vcf_out = pysam.VariantFile(args.output, 'w', header=vcf_in.header)

# parse depth cutoffs, sex assignments and sex-linked scaffolds
if args.samples:
    samples = args.samples.split(',')
else:
    samples = list(vcf_in.header.samples)
sex_assignments = {sample: 'homozygous' for sample in samples}  # default to homozygous if not provided
if args.sex_assignments:
    sex_assignments_list = args.sex_assignments.split(',')
    for sample, assignment in zip(samples, sex_assignments_list):
        sex_assignments[sample] = assignment
sex_linked_scaffolds = set(args.sex_linked_scaffolds.split(',')) if args.sex_linked_scaffolds else set()

# print warnings/errors if there are mismatches between sex assignments/number of samples/sex-linked contigs
if args.sex_assignments and len(sex_assignments) != len(samples):
    raise ValueError(f"Number of sex assignments provided does not match number of samples. Please provide sex assignments for all samples or omit the argument to default to homozygous for all samples.")
if len(sex_linked_scaffolds) > 0 and args.sex_linked_scaffolds is not None and args.sex_assignments is None:
    print(f"Warning: Providing sex-linked scaffolds without sex assignments will not have any affect at all, all scaffolds will be filtered with the provided min and max depth thresholds without any adjustments.")
if len(sex_linked_scaffolds) == 0 and args.sex_assignments and args.sex_assignments is not None:
    print(f"Warning: No sex-linked scaffolds provided, sex assignment argument will have no effect.")

# Iterate through each variant in the input VCF
for record in vcf_in:
    # Calculate allele balance for each sample
    for sample in record.samples:
        # if sample is hemizygous and we're on a sex_linked contig, adjust depth thresholds to half
        if record.chrom in sex_linked_scaffolds and sex_assignments[sample] == 'hemizygous':
            min_depth = args.min_depth / 2
            max_depth = args.max_depth / 2
        else:
            min_depth = args.min_depth
            max_depth = args.max_depth
        gt = record.samples[sample]['GT']
        if None in gt:
            continue  # Skip if genotype is missing
        ad = record.samples[sample]['AD']
        if not None in ad:
            dp = sum(ad)
        else:            
            try:
                dp = record.samples[sample]['DP']
            except KeyError:
                dp = 0
        if dp < min_depth or dp > max_depth:
            # set genotype to missing if depth is outside thresholds
            record.samples[sample]['GT'] = (None, None)
        # if heterozygous, check allele balance
        elif gt[0] != gt[1]: # Check for heterozygous genotype
            ad = record.samples[sample]['AD']  # Allele depth
            #try:
            if ad is not None and sum(ad) > 0:  # Avoid division by zero
                ab = min(ad) / sum(ad)  # Calculate allele balance
                if ab < args.min_ab:  # Filter based on allele balance threshold
                    # set genotype to missing if allele balance is below threshold
                    record.samples[sample]['GT'] = (None, None)
            # except:
            #     print(record)
    # Write the (potentially modified) record to the output VCF
    vcf_out.write(record)