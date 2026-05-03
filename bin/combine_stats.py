#!/usr/bin/env python3
import argparse
import random

args = argparse.ArgumentParser(description="Combine multiple outputs from vcfstats.py.")

args.add_argument("--ab-dp", help="Comma separated list of allele balance and genotype depth files.", required=True)
args.add_argument("--gt-stats", help="Comma separated list of genotype stats files.", required=True)
args.add_argument("--sample-stats", help="Comma separated list of sample summary stats files.", required=True)
args.add_argument("--rec-counts", help="Comma separated list of record type counts files.", required=True)
args.add_argument("-o", "--output", help="Output files prefix.", required=True)
args.add_argument("--max-records", help="Maximum number of records to write for the ab_dp and qual_fmiss_maf files. If input exceeds this number, we'll randomly sample this many records. Defaults to 100K.", type=int, default=10000)

args = args.parse_args()

# check that the number of files match
ab_dp_files = args.ab_dp.split(",")
gt_stats_files = args.gt_stats.split(",")
sample_stats_files = args.sample_stats.split(",")
rec_counts_files = args.rec_counts.split(",")
if not (len(ab_dp_files) == len(gt_stats_files) == len(sample_stats_files) == len(rec_counts_files)):
    raise ValueError("The number of files provided for each category must match.")

# prep data structures for combined stats
combined_ab_dp = {}
combined_gt_stats = []
combined_sample_stats = {}
combined_rec_counts = {}

samples = []

# ab stats, start with reading the first file to get sample names
with open(ab_dp_files[0], "r") as f:
    for line in f:
        items = line.strip().split("\t")
        stat = items[0]
        sample = items[1]
        values = [float(x) for x in items[2:]]
        if sample not in combined_ab_dp:
            combined_ab_dp[sample] = {}
        combined_ab_dp[sample][stat] = values
        if not sample in samples:
            samples.append(sample)

# read the rest of the ab files
for file in ab_dp_files[1:]:
    _samples = []
    with open(file, "r") as f:
        for line in f:
            items = line.strip().split("\t")
            stat = items[0]
            sample = items[1]
            values = [float(x) for x in items[2:]]
            if sample not in combined_ab_dp:
                raise ValueError(f"Sample {sample} found in file {file} not found in the first ab-dp file.")
            combined_ab_dp[sample][stat] = values
            _samples.append(sample)
    if not set(_samples) == set(samples):
        raise ValueError(f"Samples in file {file} do not match samples in the first ab-dp file.")
# gt stats
for file in gt_stats_files:
    with open(file, "r") as f:
        for i,line in enumerate(f):
            if i == 0:
                continue
            items = line.strip().split("\t")
            qual = items[0]
            fmiss = items[1]
            maf = items[2]
            dp = items[3]
            combined_gt_stats.append((qual,fmiss,maf,dp))
            
# sample stats
for file in sample_stats_files:
    _samples = []
    with open(file, "r") as f:
        for i,line in enumerate(f):
            if i == 0:
                continue
            items = line.strip().split("\t")
            sample = items[0]
            num_records = int(items[1])
            num_hom_ref = int(items[2])
            num_het = int(items[3])
            num_hom_alt = int(items[4])
            num_missing = int(items[5])
            if sample not in combined_sample_stats:
                combined_sample_stats[sample] = {"num_records": num_records, "num_hom_ref": num_hom_ref, "num_het": num_het, "num_hom_alt": num_hom_alt, "num_missing": num_missing}
            else:
                combined_sample_stats[sample]["num_records"] += num_records
                combined_sample_stats[sample]["num_hom_ref"] += num_hom_ref
                combined_sample_stats[sample]["num_het"] += num_het
                combined_sample_stats[sample]["num_hom_alt"] += num_hom_alt
                combined_sample_stats[sample]["num_missing"] += num_missing
            _samples.append(sample)
    if not set(_samples) == set(samples):
        raise ValueError(f"Samples in file {file} do not match samples in the first ab-dp file.")
# record counts
for file in rec_counts_files:
    with open(file, "r") as f:
        for i,line in enumerate(f):
            if i == 0:
                continue
            items = line.strip().split("\t")
            rec_type = items[0]
            count = items[1]
            if rec_type in combined_rec_counts:
                combined_rec_counts[rec_type] += int(count)
            else:
                combined_rec_counts[rec_type] = int(count)

# subsample ab_dp and gt stats if needed
total_records = len(combined_gt_stats)
if total_records > args.max_records:
    sampled_items = random.sample(combined_gt_stats, args.max_records)
    combined_gt_stats = sampled_items

for sample in samples:
    total_ab = len(combined_ab_dp[sample]["allele_balance"])
    total_dp = len(combined_ab_dp[sample]["genotype_depth"])
    if total_ab > args.max_records:
        sampled_ab = random.sample(combined_ab_dp[sample]["allele_balance"], args.max_records)
        combined_ab_dp[sample]["allele_balance"] = sampled_ab
    if total_dp > args.max_records:
        sampled_dp = random.sample(combined_ab_dp[sample]["genotype_depth"], args.max_records)
        combined_ab_dp[sample]["genotype_depth"] = sampled_dp


# write combined stats to files
with open(f"{args.output}_ab.tsv", "w") as f:
    print(f)
    f.write("sample\tminor_allele_support\n")
    for sample in samples:
        for i in combined_ab_dp[sample]["allele_balance"]:
            f.write(sample + "\t" + str(i) + "\n")
with open(f"{args.output}_dp.tsv", "w") as f:
    f.write("sample\tgenotype_depth\n")
    for sample in samples:
       for i in combined_ab_dp[sample]["genotype_depth"]:
           f.write(sample + "\t" + str(i) + "\n")
with open(f"{args.output}_qual_fmiss_maf_dp.tsv", "w") as f:
    f.write("qual\tfmiss\tmaf\tdp\n")
    for qual, fmiss, maf, dp in combined_gt_stats:
        f.write(qual + "\t" + fmiss + "\t" + maf + "\t" + dp + "\n")
with open(f"{args.output}_sample_stats.tsv", "w") as f:
    f.write("sample\tnum_records\tnum_hom_ref\tnum_het\tnum_hom_alt\tnum_missing\n")
    for sample in samples:
        f.write(sample + "\t" + str(combined_sample_stats[sample]["num_records"]) + "\t" + str(combined_sample_stats[sample]["num_hom_ref"]) + "\t" + str(combined_sample_stats[sample]["num_het"]) + "\t" + str(combined_sample_stats[sample]["num_hom_alt"]) + "\t" + str(combined_sample_stats[sample]["num_missing"]) + "\n")
with open(f"{args.output}_record_counts.tsv", "w") as f:
    f.write("record_type\tcount\n")
    for rec_type in combined_rec_counts:
        f.write(rec_type + "\t" + str(combined_rec_counts[rec_type]) + "\n")
