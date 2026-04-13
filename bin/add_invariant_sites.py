#!/usr/bin/env python3

# this scripts takes a vcf file with snps and a bedfile with homozygous reference calls, and adds invariant sites to the vcf file, so that psmc can be run on the resulting file vcf file
# optionally provide a region upon which to work, in the format "chr:start-end"

import argparse
from pysam import FastaFile, VariantFile
from pybedtools import BedTool, Interval

def format_homref_interval(interval, ref_fasta):
    chrom = interval.chrom
    start = interval.start
    end = interval.end
    ref = ref_fasta.fetch(chrom, start, end)
    
    # Create tab-separated string with all fields
    fields = [chrom, str(start), str(end), ref, '.', 'PASS', '.','GT', '0/0']
    line = '\t'.join(fields)
    
    return line


def filter_bed(bed, chrom, start, end):
    if not chrom:
        return bed
    if start and end:
        return bed.filter(lambda x: x.chrom == chrom and int(x.start) >= start and int(x.end) <= end)
    else:
        return bed.filter(lambda x: x.chrom == chrom)

def fetch_refallele(chrom, pos, fastafile):
    # fetch reference allele from fasta file
    return fastafile.fetch(chrom, pos, pos+1)

args = argparse.ArgumentParser(description="Add invariant sites to a vcf file, given a bed file with homozygous reference calls")
args.add_argument("--vcf", help="Input vcf file with snps", required=True)
args.add_argument("--homref", help="Input bed file with homozygous reference calls", required=True)
args.add_argument("--reference", help="Reference fasta file", required=True)
args.add_argument("--output", help="Output vcf file with invariant sites added", required=True)
args.add_argument("--region", help="Region to work on, in the format 'chr:start-end'", required=False)
args = args.parse_args()


homrefs = "/cfs/klemming/projects/supr/naiss2025-23-567/dev/nf-varcall/output/03_genotypes/Eurostopodus_archboldi_B10K-D-KJ-64.homozygous_reference.bed"
homrefs = args.homref
refgenome_file = "/cfs/klemming/projects/supr/naiss2025-23-567/dev/nf-varcall/output/01_reference_genome/Eurostopodus_archboldi_B10K-D-KJ-64.fasta.gz"
refgenome_file = args.reference
vcf_in = "/cfs/klemming/projects/supr/naiss2025-23-567/dev/nf-varcall/output/03_genotypes/snps.vcf.gz"
vcf_in = args.vcf
region = "scaffold_37:1-1000000"
region = args.region
chrom = None
start = None
end = None
if region:
    if len(region.split(":")) == 2:
        chrom, positions = region.split(":")
        if len(positions.split("-")) == 2:
            start, end = positions.split("-")
            region = (chrom, int(start) -1 , int(end))
        else:
            start = int(positions) - 1
            end = start
    else:
        chrom = region
        start = 0
        # fetch the length of the chromosome from the fasta file
        with FastaFile(refgenome_file) as f:
            end = f.get_reference_length(chrom)
    # make a bedlike region string for filtering the bed file
    region_bed = BedTool(f"{chrom}\t{start}\t{end}", from_string=True)


# fetch the vcf file
vcf_bed = BedTool(vcf_in)

# open bed file
bed = BedTool(homrefs)

# open fasta reference
ref = FastaFile(refgenome_file)

# filter if needed
if region:
    filtered_homref_bed = bed.intersect(region_bed)
    filtered_vcf_bed = vcf_bed.intersect(region_bed)

# Add fields to homref entries
records = [format_homref_interval(x, ref) for x in filtered_homref_bed]

for rec in filtered_vcf_bed:
    records.append('\t'.join([rec[0], str(int(rec[1])-1), rec[1], rec[3], rec[4], rec[5], rec[6], rec[7], rec[8]]))

combined_bed = BedTool('\n'.join(records) + '\n', from_string=True).sort()
for i in combined_bed:
    print(i)

for i in combined_bed:
    chrom = i.chrom
    start = i.start
    end = i.end
    for j in range(int(start), int(end) + 1):
        print(j)
        ref_allele = fetch_refallele(chrom, j, ref)
        print(f"{chrom}\t{str(j)}\t.\t{ref_allele}\t{ref_allele}\t.\tPASS\t.")

