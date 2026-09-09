#!/usr/bin/env python3
"""Combine per-chunk vcfstats.py outputs into genome-wide summary statistics.

Each vcfstats.py chunk output already caps its distribution-shaped values
(allele balance, genotype depth, QUAL, F_MISS, MAF) to a bounded reservoir
sample plus exact streaming moments (n/sum/sumsq/min/max) -- see
vcfstats.py's StatsAccumulator. This script merges those bounded per-chunk
pieces genome-wide: moments merge in O(1) per chunk, and reservoir samples
are concatenated (O(cap) per chunk) before one final random subsample -- so
memory here scales with num_chunks * cap, not with the total number of
genotypes/records genome-wide.
"""

import argparse
import random


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_optional_float(value):
    return None if value == "NA" else float(value)


def parse_ab_dp_file(path):
    """Parse one vcfstats.py `_ab_dp.txt` chunk file.

    Returns (samples, values, moments):
      samples: sample names in the order first encountered in the file.
      values: {(sample, metric): [bounded reservoir sample from this chunk]}
      moments: {(sample, metric): {n, sum, sumsq, min, max}}
    """
    samples = []
    seen = set()
    values = {}
    moments = {}
    with open(path) as f:
        for line in f:
            items = line.rstrip("\n").split("\t")
            tag, sample = items[0], items[1]
            rest = items[2:]
            if sample not in seen:
                seen.add(sample)
                samples.append(sample)
            if tag.endswith("_moments"):
                metric = tag[: -len("_moments")]
                n, s, sumsq, vmin, vmax = rest
                moments[(sample, metric)] = {
                    "n": int(n),
                    "sum": float(s),
                    "sumsq": float(sumsq),
                    "min": _parse_optional_float(vmin),
                    "max": _parse_optional_float(vmax),
                }
            else:
                metric = tag
                values.setdefault((sample, metric), []).extend(
                    float(v) for v in rest if v != "")
    return samples, values, moments


def parse_qual_fmiss_maf_dp_file(path):
    """Parse one vcfstats.py `_qual_fmiss_maf_dp.txt` chunk file.

    Returns (rows, moments):
      rows: [(qual, fmiss, maf, dp), ...] -- this chunk's bounded reservoir
        sample of records, values or None (for 'NA') per field.
      moments: {metric: {n, sum, sumsq, min, max}}
    """
    rows = []
    moments = {}
    with open(path) as f:
        lines = iter(f)
        next(lines, None)  # header: QUAL  F_MISS  MAF  DP
        for line in lines:
            line = line.rstrip("\n")
            if line == "#MOMENTS":
                break
            items = line.split("\t")
            rows.append(tuple(_parse_optional_float(v) for v in items))
        next(lines, None)  # moments header: metric  n  sum  sumsq  min  max
        for line in lines:
            metric, n, s, sumsq, vmin, vmax = line.rstrip("\n").split("\t")
            moments[metric] = {
                "n": int(n),
                "sum": float(s),
                "sumsq": float(sumsq),
                "min": _parse_optional_float(vmin),
                "max": _parse_optional_float(vmax),
            }
    return rows, moments


def parse_sample_stats_file(path):
    """Parse one vcfstats.py `_sample_sumstats.txt` chunk file.

    Returns {sample: {num_records, num_hom_ref, num_het, num_hom_alt, num_missing}}.
    """
    keys = ("num_records", "num_hom_ref", "num_het", "num_hom_alt", "num_missing")
    stats = {}
    with open(path) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            items = line.rstrip("\n").split("\t")
            sample = items[0]
            stats[sample] = dict(zip(keys, (int(x) for x in items[1:6])))
    return stats


def parse_rec_counts_file(path):
    """Parse one vcfstats.py `_rec_counts.txt` chunk file. Returns {rec_type: count}."""
    counts = {}
    with open(path) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            rec_type, count = line.rstrip("\n").split("\t")
            counts[rec_type] = int(count)
    return counts


# ---------------------------------------------------------------------------
# Merging
# ---------------------------------------------------------------------------

def merge_moments(moments_list):
    """Merge a list of {n, sum, sumsq, min, max} dicts into one, in O(len(list))."""
    total = {"n": 0, "sum": 0.0, "sumsq": 0.0, "min": None, "max": None}
    for m in moments_list:
        if m["n"] == 0:
            continue
        total["n"] += m["n"]
        total["sum"] += m["sum"]
        total["sumsq"] += m["sumsq"]
        total["min"] = m["min"] if total["min"] is None else min(total["min"], m["min"])
        total["max"] = m["max"] if total["max"] is None else max(total["max"], m["max"])
    return total


def moments_mean(m):
    return m["sum"] / m["n"] if m["n"] else None


def moments_sd(m):
    if m["n"] < 2:
        return None
    variance = m["sumsq"] / m["n"] - moments_mean(m) ** 2
    return variance ** 0.5 if variance > 0 else 0.0


def merge_reservoir_values(value_lists, max_records, rng=None):
    """Concatenate bounded per-chunk samples and subsample down to max_records."""
    rng = rng if rng is not None else random.Random()
    combined = [v for values in value_lists for v in values]
    if len(combined) > max_records:
        combined = rng.sample(combined, max_records)
    return combined


# ---------------------------------------------------------------------------
# Per-category combination across chunk files
# ---------------------------------------------------------------------------

def combine_ab_dp_files(paths):
    """Parse and merge multiple `_ab_dp.txt` chunk files.

    Returns (samples, combined_values, combined_moments) -- see
    parse_ab_dp_file for the shape of values/moments; combined_values is the
    concatenation (not yet subsampled) across all chunks.
    """
    samples = None
    all_values = {}
    all_moments = {}
    for path in paths:
        file_samples, values, moments = parse_ab_dp_file(path)
        if samples is None:
            samples = file_samples
        elif set(file_samples) != set(samples):
            raise ValueError(
                f"Samples in file {path} do not match samples in the first ab-dp file.")
        for key, vals in values.items():
            all_values.setdefault(key, []).extend(vals)
        for key, m in moments.items():
            all_moments.setdefault(key, []).append(m)
    merged_moments = {key: merge_moments(m_list) for key, m_list in all_moments.items()}
    return samples or [], all_values, merged_moments


def combine_gt_stats_files(paths):
    """Parse and merge multiple `_qual_fmiss_maf_dp.txt` chunk files.

    Returns (rows, merged_moments); rows is the concatenation (not yet
    subsampled) of every chunk's bounded reservoir sample.
    """
    all_rows = []
    all_moments = {}
    for path in paths:
        rows, moments = parse_qual_fmiss_maf_dp_file(path)
        all_rows.extend(rows)
        for metric, m in moments.items():
            all_moments.setdefault(metric, []).append(m)
    merged_moments = {metric: merge_moments(m_list) for metric, m_list in all_moments.items()}
    return all_rows, merged_moments


def combine_sample_stats_files(paths, expected_samples=None):
    """Parse and exactly sum multiple `_sample_sumstats.txt` chunk files."""
    keys = ("num_records", "num_hom_ref", "num_het", "num_hom_alt", "num_missing")
    combined = {}
    for path in paths:
        file_stats = parse_sample_stats_file(path)
        if expected_samples and set(file_stats) != set(expected_samples):
            raise ValueError(
                f"Samples in file {path} do not match samples in the first ab-dp file.")
        for sample, stats in file_stats.items():
            target = combined.setdefault(sample, {k: 0 for k in keys})
            for key in keys:
                target[key] += stats[key]
    return combined


def combine_rec_counts_files(paths):
    """Parse and exactly sum multiple `_rec_counts.txt` chunk files."""
    combined = {}
    for path in paths:
        for rec_type, count in parse_rec_counts_file(path).items():
            combined[rec_type] = combined.get(rec_type, 0) + count
    return combined


# ---------------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------------

def _format_value(value):
    return "NA" if value is None else str(value)


def write_ab(prefix, samples, values):
    with open(f"{prefix}_ab.tsv", "w") as f:
        f.write("sample\tminor_allele_support\n")
        for sample in samples:
            for v in values.get((sample, "allele_balance"), []):
                f.write(f"{sample}\t{v}\n")


def write_dp(prefix, samples, values):
    with open(f"{prefix}_dp.tsv", "w") as f:
        f.write("sample\tgenotype_depth\n")
        for sample in samples:
            for v in values.get((sample, "genotype_depth"), []):
                f.write(f"{sample}\t{v}\n")


def write_qual_fmiss_maf_dp(prefix, rows):
    with open(f"{prefix}_qual_fmiss_maf_dp.tsv", "w") as f:
        f.write("qual\tfmiss\tmaf\tdp\n")
        for qual, fmiss, maf, dp in rows:
            f.write("\t".join(_format_value(v) for v in (qual, fmiss, maf, dp)) + "\n")


def write_sample_stats(prefix, samples, stats):
    with open(f"{prefix}_sample_stats.tsv", "w") as f:
        f.write("sample\tnum_records\tnum_hom_ref\tnum_het\tnum_hom_alt\tnum_missing\n")
        for sample in samples:
            s = stats[sample]
            f.write(f"{sample}\t{s['num_records']}\t{s['num_hom_ref']}\t{s['num_het']}\t"
                    f"{s['num_hom_alt']}\t{s['num_missing']}\n")


def write_record_counts(prefix, counts):
    with open(f"{prefix}_record_counts.tsv", "w") as f:
        f.write("record_type\tcount\n")
        for rec_type, count in counts.items():
            f.write(f"{rec_type}\t{count}\n")


def write_dist_moments(prefix, samples, ab_dp_moments, site_moments):
    """Bonus output: exact genome-wide n/mean/sd/min/max, not sample-derived.

    Not currently consumed by plot_variantstats.R, but picked up automatically
    by combine_stats.nf's `path("${category}*")` glob for future use.
    """
    with open(f"{prefix}_dist_moments.tsv", "w") as f:
        f.write("scope\tkey\tmetric\tn\tmean\tsd\tmin\tmax\n")
        for metric in ("qual", "fmiss", "maf", "dp"):
            m = site_moments.get(metric)
            if m:
                f.write(f"site\tALL\t{metric}\t{m['n']}\t{_format_value(moments_mean(m))}\t"
                        f"{_format_value(moments_sd(m))}\t{_format_value(m['min'])}\t"
                        f"{_format_value(m['max'])}\n")
        for sample in samples:
            for metric in ("allele_balance", "genotype_depth"):
                m = ab_dp_moments.get((sample, metric))
                if m:
                    f.write(f"sample\t{sample}\t{metric}\t{m['n']}\t"
                            f"{_format_value(moments_mean(m))}\t{_format_value(moments_sd(m))}\t"
                            f"{_format_value(m['min'])}\t{_format_value(m['max'])}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Combine multiple outputs from vcfstats.py.")
    parser.add_argument("--ab-dp", required=True,
                         help="Comma separated list of allele balance and genotype depth files.")
    parser.add_argument("--gt-stats", required=True,
                         help="Comma separated list of genotype/quality stats files.")
    parser.add_argument("--sample-stats", required=True,
                         help="Comma separated list of sample summary stats files.")
    parser.add_argument("--rec-counts", required=True,
                         help="Comma separated list of record type counts files.")
    parser.add_argument("-o", "--output", required=True, help="Output files prefix.")
    parser.add_argument(
        "--max-records", type=int, default=10000,
        help="Maximum number of records to write for the ab_dp and qual_fmiss_maf "
             "files. If the combined per-chunk reservoir samples exceed this "
             "number, randomly subsample down to it. Defaults to 10000.")
    parser.add_argument("--seed", type=int, default=None,
                         help="Random seed for the final subsample.")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    ab_dp_files = args.ab_dp.split(",")
    gt_stats_files = args.gt_stats.split(",")
    sample_stats_files = args.sample_stats.split(",")
    rec_counts_files = args.rec_counts.split(",")
    if not (len(ab_dp_files) == len(gt_stats_files)
            == len(sample_stats_files) == len(rec_counts_files)):
        raise ValueError("The number of files provided for each category must match.")

    rng = random.Random(args.seed)

    samples, ab_dp_values, ab_dp_moments = combine_ab_dp_files(ab_dp_files)
    site_rows, site_moments = combine_gt_stats_files(gt_stats_files)
    sample_stats = combine_sample_stats_files(sample_stats_files, samples)
    rec_counts = combine_rec_counts_files(rec_counts_files)

    ab_dp_values = {
        key: merge_reservoir_values([vals], args.max_records, rng)
        for key, vals in ab_dp_values.items()
    }
    site_rows = merge_reservoir_values([site_rows], args.max_records, rng)

    write_ab(args.output, samples, ab_dp_values)
    write_dp(args.output, samples, ab_dp_values)
    write_qual_fmiss_maf_dp(args.output, site_rows)
    write_sample_stats(args.output, samples, sample_stats)
    write_record_counts(args.output, rec_counts)
    write_dist_moments(args.output, samples, ab_dp_moments, site_moments)


if __name__ == "__main__":
    main()
