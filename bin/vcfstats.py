#!/usr/bin/env python3
"""Compute sample- and record-level statistics for one VCF (or VCF chunk).

This is designed to run once per chunk in a larger pipeline, where many
chunks' outputs are later merged genome-wide by combine_stats.py. To keep
that merge step's memory bounded regardless of how many chunks or records
there are, distribution-shaped values (allele balance, genotype depth, QUAL,
F_MISS, MAF) are NOT written out in full here. Instead, each is kept as:

  * a bounded reservoir sample (see ReservoirSampler) -- a fixed-size, uniform
    random sample of the values seen, used for plotting distributions/violins.
  * exact streaming moments (see RunningMoments) -- n/sum/sum-of-squares/
    min/max, which combine_stats.py can merge across chunks in O(1) per
    chunk to get exact genome-wide mean/SD/min/max without ever holding all
    the raw values in memory.

Per-sample genotype tallies and record type counts are cheap to combine
exactly (just sums), so those are still written out in full.
"""

import argparse
import random
from dataclasses import dataclass
from typing import Optional

import pysam


# ---------------------------------------------------------------------------
# Pure classification helpers
# ---------------------------------------------------------------------------

def classify_record(record):
    """Classify a VCF record by variant type and allele count.

    Returns one of: 'no_alt', 'snp_biallelic', 'snp_multiallelic',
    'indel_biallelic', 'indel_multiallelic'.
    """
    if not record.alts:
        return 'no_alt'
    is_snp = len(record.ref) == 1 and all(len(alt) == 1 for alt in record.alts)
    biallelic = len(record.alts) == 1
    if is_snp:
        return 'snp_biallelic' if biallelic else 'snp_multiallelic'
    return 'indel_biallelic' if biallelic else 'indel_multiallelic'


def classify_genotype(gt):
    """Classify a diploid genotype tuple.

    Returns one of: 'missing', 'hom_ref', 'hom_alt', 'het'.
    """
    if gt is None or None in gt:
        return 'missing'
    if gt[0] == gt[1]:
        return 'hom_ref' if gt[0] == 0 else 'hom_alt'
    return 'het'


def safe_ad_sum(ad):
    """Sum an AD (allelic depth) tuple, skipping any None entries.

    Returns None if ad is None or every entry is missing.
    """
    if ad is None:
        return None
    values = [v for v in ad if v is not None]
    if not values:
        return None
    return sum(values)


def compute_allele_balance(ad):
    """Compute minor-allele support (min(AD[0], AD[1]) / sum(AD)) for a het call.

    Returns None (rather than raising) if AD is missing, doesn't cover both
    alleles, or sums to zero -- i.e. whenever allele balance genuinely can't
    be computed from the data.
    """
    if ad is None or len(ad) < 2 or ad[0] is None or ad[1] is None:
        return None
    total = safe_ad_sum(ad)
    if not total:
        return None
    return round(min(ad[0], ad[1]) / total, 3)


# ---------------------------------------------------------------------------
# Bounded stream summarization
# ---------------------------------------------------------------------------

class ReservoirSampler:
    """Fixed-capacity uniform random sample of a stream (Algorithm R).

    After any number of calls to add(), `values` is a uniform random sample
    of at most `capacity` items from everything ever added -- memory never
    grows past `capacity`, regardless of stream length.
    """

    def __init__(self, capacity, rng=None):
        if capacity < 0:
            raise ValueError("capacity must be >= 0")
        self.capacity = capacity
        self.rng = rng if rng is not None else random.Random()
        self._values = []
        self.n_seen = 0

    def add(self, value):
        self.n_seen += 1
        if len(self._values) < self.capacity:
            self._values.append(value)
        else:
            j = self.rng.randrange(self.n_seen)
            if j < self.capacity:
                self._values[j] = value

    @property
    def values(self):
        return list(self._values)


@dataclass
class RunningMoments:
    """Streaming n/sum/sum-of-squares/min/max.

    Mergeable across chunks in O(1) (see merge()), which is what makes exact
    genome-wide mean/SD/min/max possible without ever concatenating raw
    per-chunk values.
    """

    n: int = 0
    sum: float = 0.0
    sumsq: float = 0.0
    min: Optional[float] = None
    max: Optional[float] = None

    def add(self, value):
        self.n += 1
        self.sum += value
        self.sumsq += value * value
        self.min = value if self.min is None else min(self.min, value)
        self.max = value if self.max is None else max(self.max, value)

    @property
    def mean(self):
        return self.sum / self.n if self.n else None

    @property
    def sd(self):
        if self.n < 2:
            return None
        variance = self.sumsq / self.n - self.mean ** 2
        return variance ** 0.5 if variance > 0 else 0.0

    def merge(self, other):
        """Return a new RunningMoments combining self and other."""
        if self.n == 0:
            return RunningMoments(other.n, other.sum, other.sumsq, other.min, other.max)
        if other.n == 0:
            return RunningMoments(self.n, self.sum, self.sumsq, self.min, self.max)
        return RunningMoments(
            n=self.n + other.n,
            sum=self.sum + other.sum,
            sumsq=self.sumsq + other.sumsq,
            min=min(self.min, other.min),
            max=max(self.max, other.max),
        )


# ---------------------------------------------------------------------------
# Per-sample genotype tallies (exact, cheap -- unaffected by the above)
# ---------------------------------------------------------------------------

@dataclass
class SampleGenotypeCounts:
    num_records: int = 0
    num_hom_ref: int = 0
    num_het: int = 0
    num_hom_alt: int = 0
    num_missing: int = 0

    _ATTR_BY_CLASS = {
        'hom_ref': 'num_hom_ref',
        'het': 'num_het',
        'hom_alt': 'num_hom_alt',
        'missing': 'num_missing',
    }

    def update(self, genotype_class):
        self.num_records += 1
        attr = self._ATTR_BY_CLASS[genotype_class]
        setattr(self, attr, getattr(self, attr) + 1)


# ---------------------------------------------------------------------------
# Accumulator: single stateful object driving the per-record processing
# ---------------------------------------------------------------------------

class StatsAccumulator:
    """Accumulates per-sample and per-record VCF statistics for one chunk."""

    REC_TYPES = (
        'snp_biallelic', 'snp_multiallelic',
        'indel_biallelic', 'indel_multiallelic', 'no_alt',
    )
    SAMPLE_METRICS = ('allele_balance', 'genotype_depth')
    SITE_METRICS = ('qual', 'fmiss', 'maf', 'dp')

    def __init__(self, samples, max_sampled_values=5000, rng=None):
        self.samples = list(samples)
        rng = rng if rng is not None else random.Random()

        self.rec_type_counts = {rt: 0 for rt in self.REC_TYPES}
        self.sample_counts = {s: SampleGenotypeCounts() for s in self.samples}
        self.sample_moments = {
            s: {m: RunningMoments() for m in self.SAMPLE_METRICS} for s in self.samples
        }
        self.sample_reservoirs = {
            s: {m: ReservoirSampler(max_sampled_values, rng) for m in self.SAMPLE_METRICS}
            for s in self.samples
        }
        self.site_moments = {m: RunningMoments() for m in self.SITE_METRICS}
        self.site_reservoir = ReservoirSampler(max_sampled_values, rng)

    def process_record(self, record):
        rec_type = classify_record(record)
        self.rec_type_counts[rec_type] += 1
        if rec_type in ('snp_biallelic', 'indel_biallelic'):
            self._process_biallelic(record)
        else:
            self._process_other(record)

    def _process_other(self, record):
        for sample in self.samples:
            self._require_sample(record, sample)
            gt = record.samples[sample]["GT"]
            self.sample_counts[sample].update(classify_genotype(gt))

    def _process_biallelic(self, record):
        called_genotypes = 0
        allele_counts = [0, 0]  # [ref alleles, alt alleles] across called genotypes

        for sample in self.samples:
            self._require_sample(record, sample)
            gt = record.samples[sample]["GT"]
            gt_class = classify_genotype(gt)
            self.sample_counts[sample].update(gt_class)
            if gt_class == 'missing':
                continue

            called_genotypes += 1
            ad = record.samples[sample]["AD"]
            depth = safe_ad_sum(ad)
            if depth is not None:
                self._add_sample_value(sample, 'genotype_depth', depth)

            if gt_class == 'hom_ref':
                allele_counts[0] += 2
            elif gt_class == 'hom_alt':
                allele_counts[1] += 2
            else:  # het
                allele_counts[0] += 1
                allele_counts[1] += 1
                ab = compute_allele_balance(ad)
                # Preserve historical behavior: an unmeasurable allele balance
                # on a het call is recorded as 0 rather than dropped.
                self._add_sample_value(sample, 'allele_balance', ab if ab is not None else 0)

        qual = round(record.qual, 3) if record.qual is not None else None
        dp = record.info['DP'] if 'DP' in record.info else None

        if called_genotypes > 0:
            fmiss = round(1 - called_genotypes / len(self.samples), 3)
            total_alleles = sum(allele_counts)
            maf = round(min(allele_counts) / total_alleles, 3) if total_alleles else None
        else:
            fmiss = 1.0
            maf = None

        self._add_site_value('qual', qual)
        self._add_site_value('fmiss', fmiss)
        self._add_site_value('maf', maf)
        self._add_site_value('dp', dp)
        self.site_reservoir.add((qual, fmiss, maf, dp))

    def _add_sample_value(self, sample, metric, value):
        self.sample_moments[sample][metric].add(value)
        self.sample_reservoirs[sample][metric].add(value)

    def _add_site_value(self, metric, value):
        if value is not None:
            self.site_moments[metric].add(value)

    @staticmethod
    def _require_sample(record, sample):
        # pysam's VariantRecordSamples raises KeyError from `in` itself for an
        # unknown sample name rather than returning False, so this can't be a
        # plain `if sample not in record.samples:` check.
        try:
            found = sample in record.samples
        except KeyError:
            found = False
        if not found:
            raise ValueError(f"Sample {sample} not found in VCF file.")


# ---------------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------------

def _format_value(value):
    return 'NA' if value is None else str(value)


def write_outputs(accumulator, output_prefix):
    _write_sample_sumstats(accumulator, output_prefix)
    _write_rec_counts(accumulator, output_prefix)
    _write_ab_dp(accumulator, output_prefix)
    _write_qual_fmiss_maf_dp(accumulator, output_prefix)


def _write_sample_sumstats(acc, prefix):
    with open(f"{prefix}_sample_sumstats.txt", "w") as f:
        f.write("Sample\tNum_Records\tNum_Hom_Ref\tNum_Het\tNum_Hom_Alt\tNum_Missing\n")
        for sample in acc.samples:
            c = acc.sample_counts[sample]
            f.write(f"{sample}\t{c.num_records}\t{c.num_hom_ref}\t{c.num_het}\t"
                    f"{c.num_hom_alt}\t{c.num_missing}\n")


def _write_rec_counts(acc, prefix):
    with open(f"{prefix}_rec_counts.txt", "w") as f:
        f.write("Record_Type\tCount\n")
        for rec_type in StatsAccumulator.REC_TYPES:
            f.write(f"{rec_type}\t{acc.rec_type_counts[rec_type]}\n")


def _write_ab_dp(acc, prefix):
    # Same `stat<TAB>sample<TAB>v1<TAB>v2...` shape as before, now bounded to
    # a reservoir sample, plus a "<metric>_moments" row per (sample, metric)
    # carrying the exact n/sum/sumsq/min/max for that metric.
    with open(f"{prefix}_ab_dp.txt", "w") as f:
        for sample in acc.samples:
            for metric in StatsAccumulator.SAMPLE_METRICS:
                values = acc.sample_reservoirs[sample][metric].values
                f.write(metric + "\t" + sample + "\t" + "\t".join(map(str, values)) + "\n")
                m = acc.sample_moments[sample][metric]
                f.write(f"{metric}_moments\t{sample}\t{m.n}\t{m.sum}\t{m.sumsq}\t"
                        f"{_format_value(m.min)}\t{_format_value(m.max)}\n")


def _write_qual_fmiss_maf_dp(acc, prefix):
    # Same header + per-record row shape as before, now bounded to a
    # reservoir sample, followed by a `#MOMENTS` block with the exact
    # n/sum/sumsq/min/max per metric.
    with open(f"{prefix}_qual_fmiss_maf_dp.txt", "w") as f:
        f.write("QUAL\tF_MISS\tMAF\tDP\n")
        for qual, fmiss, maf, dp in acc.site_reservoir.values:
            f.write("\t".join(_format_value(v) for v in (qual, fmiss, maf, dp)) + "\n")
        f.write("#MOMENTS\n")
        f.write("metric\tn\tsum\tsumsq\tmin\tmax\n")
        for metric in StatsAccumulator.SITE_METRICS:
            m = acc.site_moments[metric]
            f.write(f"{metric}\t{m.n}\t{m.sum}\t{m.sumsq}\t"
                    f"{_format_value(m.min)}\t{_format_value(m.max)}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Fetch some sample-based statistics from vcf file.")
    parser.add_argument("-i", "--input", help="Input vcf file.", required=True)
    parser.add_argument("-o", "--output", help="Output files prefix.", required=True)
    parser.add_argument(
        "--max-sampled-values", type=int, default=1000,
        help="Maximum number of values to keep, via reservoir sampling, per "
             "distribution-shaped metric (allele balance, genotype depth, QUAL, "
             "F_MISS, MAF) -- bounds this chunk's output size regardless of how "
             "many records it has. Exact genome-wide moments (n/sum/sumsq/min/max) "
             "are always kept in full precision, unaffected by this cap. "
             "Defaults to 1000.")
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Random seed for reservoir sampling, for reproducibility.")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    vcf = pysam.VariantFile(args.input)
    samples = list(vcf.header.samples)
    rng = random.Random(args.seed)

    accumulator = StatsAccumulator(samples, max_sampled_values=args.max_sampled_values, rng=rng)
    for record in vcf.fetch():
        accumulator.process_record(record)

    write_outputs(accumulator, args.output)


if __name__ == "__main__":
    main()
