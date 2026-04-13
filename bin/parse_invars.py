#!/usr/bin/env python3

import argparse
import pysam

args = argparse.ArgumentParser(description='Parse and filter invariant genotypes')

args.add_argument('-i', '--input', help='input vcf file', required=True)
args.add_argument('-o', '--output', help='output prefix for bed file(s).', required=True)
args.add_argument('-s', '--samples', help='comma-separated list of samples to include (default: all)', default=None)
args.add_argument('--min-depth', help='Minimum depth to keep a genotype', type=float, default=0)
args.add_argument('--max-depth', help='Maximum depth to keep a genotype', type=float, default=200)
args.add_argument('--sex-assignments', help='Comma separated list of sex assignments as hemizygous or homozygous of included samples. IMPORTANT: If multiple samples are included, the order must correspond exactly to the order of samples specified in --samples or, if this argument is omitted, the order of samples in the VCF header. Unless any sex-linked scaffolds are provided, this argument has no effect.')
args.add_argument('--sex-linked-scaffolds', help='Comma separated list of sex-linked scaffolds. If this is provided in combination with sex assignments, the coverage thresholds will be changed to half min-depth and max-depth for samples assigned as hemizygous on these scaffolds.')

args = args.parse_args()

# open the input VCF file using pysam
vcf = pysam.VariantFile(args.input)

# initialize a dict for storing invariant intervals
invars = {sample: [] for sample in args.samples.split(',')} if args.samples else {sample: [] for sample in vcf.header.samples}

# parse sex assignments and sex-linked scaffolds if provided
sex_assignments = {sample: 'homozygous' for sample in invars.keys()}  # default to homozygous if not provided
if args.sex_assignments:
    sex_assignments_list = args.sex_assignments.split(',')
    for sample, assignment in zip(invars.keys(), sex_assignments_list):
        sex_assignments[sample] = assignment
sex_linked_scaffolds = set(args.sex_linked_scaffolds.split(',')) if args.sex_linked_scaffolds else set()

if len(sex_linked_scaffolds) > 0 and args.sex_linked_scaffolds is not None:
    print(f"Warning: No sex-linked scaffolds provided, sex assignment argument will have no effect.")
if args.sex_assignments and len(sex_assignments) != len(invars):
    raise ValueError(f"Number of sex assignments provided does not match number of samples. Please provide sex assignments for all samples or omit the argument to default to homozygous for all samples.")
if args.sex_linked_scaffolds is not None and args.sex_assignments is None:
    print(f"Warning: Providing sex-linked scaffolds without sex assignments will not have any affect at all, all scaffolds will be filtered with the provided min and max depth thresholds without any adjustments.")

# initialize variables for tracking the current interval
current_intervals = {sample: (None, None, None) for sample in invars.keys()}  # (chrom, start, end)

def new_interval(current_intervals, sample, chrom, start):
    # append the last interval to the list if it exists
    if current_intervals[sample][0] is not None:
        invars[sample].append(current_intervals[sample])
    # start a new interval
    current_intervals[sample] = (chrom, start, start)
    return current_intervals

def extend_interval(current_intervals, sample, chrom, pos):
    # extend the current interval
    current_intervals[sample] = (chrom, current_intervals[sample][1], pos)
    return current_intervals

def is_consecutive(current_intervals, sample, chrom, pos):
    # check if the current position is consecutive to the last interval
    return (current_intervals[sample][0] == chrom and
            current_intervals[sample][2] is not None and
            pos == current_intervals[sample][2] + 1)

def break_interval(current_intervals, sample):
    # append the last interval to the list if it exists
    if current_intervals[sample][0] is not None:
        invars[sample].append(current_intervals[sample])
    # reset the current interval
    current_intervals[sample] = (None, None, None)
    return current_intervals

def add_interval(current_intervals, sample, chrom, pos):
    # wrapper function to check whether we should extend or start a new interval
    if is_consecutive(current_intervals, sample, chrom, pos):
        return extend_interval(current_intervals, sample, chrom, pos)
    else:        
        return new_interval(current_intervals, sample, chrom, pos)

# initialize a list of lists for storing intervals
for record in vcf.fetch():
    for sample in invars.keys():
        # if we're on a sex linked scaffold and the sample is assigned as hemizygous, adjust the depth thresholds to half
        if record.chrom in sex_linked_scaffolds and sex_assignments[sample] == 'hemizygous':
            min_depth = args.min_depth / 2
            max_depth = args.max_depth / 2
        else:
            min_depth = args.min_depth
            max_depth = args.max_depth
        if sample not in record.samples:
            continue
            #print(f"Warning: Sample {sample} not found in VCF header, skipping")
            #raise ValueError(f'Sample {sample} not found in VCF header')
        gt = record.samples[sample]['GT']
        if gt[0] == gt[1] == 0:
            # this would be a homozygous reference genotype, which is what we are looking for
            ad = record.samples[sample].get('AD', None)
            if ad is not None and None not in ad:
                # AD can be a single value or a tuple/list
                depth = sum(ad)
                if depth < min_depth or depth > max_depth:
                    current_intervals = break_interval(current_intervals, sample)
                    continue
                else:
                    current_intervals = add_interval(current_intervals, sample, record.chrom, record.pos)
            else:
                current_intervals = break_interval(current_intervals, sample)

# close any remaining open intervals
for sample in invars.keys():
    if current_intervals[sample][0] is not None:
        invars[sample].append(current_intervals[sample])

# write one bedfile per sample
for sample in invars.keys():
    with open(args.output + "_" + sample + ".bed", 'w') as f:
        for chrom, start, end in invars[sample]:
            f.write('\t'.join([chrom, str(start - 1), str(end)]) + '\n')