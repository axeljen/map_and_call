#!/usr/bin/env python3

import argparse
import pysam

args = argparse.ArgumentParser(description="Calculate allele balance for each variant in a VCF file and filter variants based on allele balance.")

args.add_argument("-i", "--input", help="Input VCF file")
args.add_argument("-o", "--output", help="Output VCF file")
args.add_argument("--min-ab", type=float, help="Minimum allele balance threshold")
args.add_argument('-s', '--samples', help='comma-separated list of samples to include (default: all)', default=None)

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

# Iterate through each variant in the input VCF
for record in vcf_in:
    # Calculate allele balance for each sample
    for sample in record.samples:
        gt = record.samples[sample]['GT']
        if gt[0] == gt[1]: # Check for homozygous genotype
            continue  # Skip if homozygous
        if None in gt:
            continue  # Skip if genotype is missing
        ad = record.samples[sample]['AD']  # Allele depth
        if ad is not None and sum(ad) > 0:  # Avoid division by zero
            ab = min(ad) / sum(ad)  # Calculate allele balance
            if ab < args.min_ab:  # Filter based on allele balance threshold
                # set genotype to missing if allele balance is below threshold
                record.samples[sample]['GT'] = (None, None)
    # Write the (potentially modified) record to the output VCF
    vcf_out.write(record)