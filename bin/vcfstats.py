#!/usr/bin/env python3

import pysam
import argparse
import sys

def record_type(record):
    if record.alts is None:
        return 'no_alt'
    if len(record.ref) == 1 and all(len(alt) == 1 for alt in record.alts):
        return 'snp_biallelic' if len(record.alts) == 1 else 'snp_multiallelic'
    elif len(record.ref) != 1 or any(len(alt) != 1 for alt in record.alts):
        return 'indel_biallelic' if len(record.alts) == 1 else 'indel_multiallelic'
    elif len(record.alts) == 0:
        return 'no_alt'
    else:
        return 'other'
    
def biallelic_snpstats(record, samples, stats, allele_balance, qual, fmiss, maf, dp):
    called_genotypes = 0
    freqs = (0,0)
    qual.append(round(record.qual, 3) if record.qual is not None else 'nan')
    dp.append(record.info['DP'] if 'DP' in record.info else 'nan')
    for sample in samples:
        if sample not in record.samples:
            raise ValueError(f"Sample {sample} not found in VCF file.")
        gt = record.samples[sample]["GT"]
        stats[sample]["num_records"] += 1
        if None in gt:
            stats[sample]["num_missing"] += 1
        elif gt[0] == gt[1] == 0:
            called_genotypes += 1
            freqs = (freqs[0] + 2, freqs[1])
            stats[sample]["num_hom_ref"] += 1
            allele_balance[sample]['genotype_depth'].append(sum([i for i in record.samples[sample]["AD"] if i is not None]))
        elif gt[0] == gt[1] and gt[0] != 0:
            called_genotypes += 1
            freqs = (freqs[0], freqs[1] + 2)
            stats[sample]["num_hom_alt"] += 1
            allele_balance[sample]['genotype_depth'].append(sum([i for i in record.samples[sample]["AD"] if i is not None]))
        elif gt[0] != gt[1]:
            called_genotypes += 1
            freqs = (freqs[0] + 1, freqs[1] + 1)
            stats[sample]["num_het"] += 1
            allele_balance[sample]['genotype_depth'].append(sum([i for i in record.samples[sample]["AD"] if i is not None]))
            try:
                allele_balance[sample]['allele_balance'].append(round(min(record.samples[sample]["AD"][0], record.samples[sample]["AD"][1]) / sum(record.samples[sample]["AD"]), 3))
            except:
                allele_balance[sample]['allele_balance'].append(0)
    if called_genotypes > 0:
        fmiss.append(round(1 - called_genotypes / len(samples), 3))
        if sum(freqs) > 0:
            if min(freqs) == 0:
                maf.append(0)
            else:
                maf.append(round(min(freqs) / sum(freqs), 3))
        else:
            maf.append('NA')
    else:
        fmiss.append(1)
        maf.append('NA')
    return stats, allele_balance, qual, fmiss, maf, dp

def other_stats(record, stats):
    for sample in samples:
        if sample not in record.samples:
            raise ValueError(f"Sample {sample} not found in VCF file.")
        gt = record.samples[sample]["GT"]
        stats[sample]["num_records"] += 1
        if None in gt:
            stats[sample]["num_missing"] += 1
        elif gt[0] == gt[1] == 0:
            stats[sample]["num_hom_ref"] += 1
        elif gt[0] == gt[1] and gt[0] != 0:
            stats[sample]["num_hom_alt"] += 1
        else:            
            stats[sample]["num_het"] += 1
    return stats

args = argparse.ArgumentParser(description="Fetch some sample-based statistics from vcf file.")
args.add_argument("-i", "--input", help="Input vcf file.", required=True)
args.add_argument("-o", "--output", help="Output files prefix.", required=True)
args = args.parse_args()

vcf = pysam.VariantFile(args.input)
samples = list(vcf.header.samples)

# set up a dictionary for storing the stats
stats = {sample: {"num_records": 0, "num_hom_ref": 0, "num_het": 0, "num_hom_alt": 0, "num_missing": 0} for sample in samples}

# allele balance and genotype depth for each sample and genotype
allele_balance = {sample: {'allele_balance': [], 'genotype_depth': []} for sample in samples}
# genotype quality, missingness and minor allele frequency for each genotype
qual = []
fmiss = []
maf = []
dp = []
# record types
rec_types = {'snp_biallelic': 0, 'snp_multiallelic': 0, 
'indel_biallelic': 0, 'indel_multiallelic': 0, 'other': 0, 'no_alt': 0,}

# count records
nrecs = 0

for record in vcf.fetch():
    nrecs += 1
    rec_types[record_type(record)] += 1
    if record_type(record) == 'snp_biallelic' or record_type(record) == 'indel_biallelic':
         stats, allele_balance, qual, fmiss, maf, dp = biallelic_snpstats(record, samples, stats, allele_balance, qual, fmiss, maf, dp)
    else:
        stats = other_stats(record, stats)    

with open(args.output + "_ab_dp.txt", "w") as f:
    for sample in samples:
        f.write(f"allele_balance\t{sample}\t" + "\t".join(map(str, allele_balance[sample]['allele_balance'])) + "\n")
        f.write(f"genotype_depth\t{sample}\t" + "\t".join(map(str, allele_balance[sample]['genotype_depth'])) + "\n")

with open(args.output + "_sample_sumstats.txt", "w") as f:
    f.write("Sample\tNum_Records\tNum_Hom_Ref\tNum_Het\tNum_Hom_Alt\tNum_Missing\n")
    for sample in samples:
        f.write(f"{sample}\t{stats[sample]['num_records']}\t{stats[sample]['num_hom_ref']}\t{stats[sample]['num_het']}\t{stats[sample]['num_hom_alt']}\t{stats[sample]['num_missing']}\n")

# write quality, missingness and minor allele frequency for each genotype
with open(args.output + "_qual_fmiss_maf_dp.txt", "w") as f:
    f.write("QUAL\tF_MISS\tMAF\tDP\n")
    for q, f_miss, maf_val, dp_val in zip(qual, fmiss, maf, dp):
        f.write(f"{q}\t{f_miss}\t{maf_val}\t{dp_val}\n")

# write record types
with open(args.output + "_rec_counts.txt", "w") as f:
    f.write("Record_Type\tCount\n")
    for rec_type, count in rec_types.items():
        f.write(f"{rec_type}\t{count}\n")