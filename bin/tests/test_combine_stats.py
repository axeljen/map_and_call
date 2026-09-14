import random
import statistics

import pytest

from combine_stats import (
    combine_ab_dp_files,
    combine_gt_stats_files,
    merge_moments,
    merge_reservoir_values,
    moments_mean,
    moments_sd,
    parse_ab_dp_file,
    parse_qual_fmiss_maf_dp_file,
)


def _raw_moments(values):
    return {
        "n": len(values),
        "sum": sum(values),
        "sumsq": sum(v * v for v in values),
        "min": min(values),
        "max": max(values),
    }


# ---------------------------------------------------------------------------
# merge_moments / moments_mean / moments_sd
# ---------------------------------------------------------------------------

def test_merge_moments_matches_concatenation():
    a = [1.0, 2.0, 3.0]
    b = [4.0, 5.0]
    merged = merge_moments([_raw_moments(a), _raw_moments(b)])
    assert merged == _raw_moments(a + b)


def test_merge_moments_skips_empty_chunks():
    empty = {"n": 0, "sum": 0.0, "sumsq": 0.0, "min": None, "max": None}
    merged = merge_moments([empty, _raw_moments([1.0, 2.0])])
    assert merged == _raw_moments([1.0, 2.0])


def test_merge_moments_of_nothing():
    assert merge_moments([]) == {"n": 0, "sum": 0.0, "sumsq": 0.0, "min": None, "max": None}


def test_moments_mean_and_sd_match_stdlib():
    values = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]
    m = _raw_moments(values)
    assert moments_mean(m) == pytest.approx(statistics.mean(values))
    assert moments_sd(m) == pytest.approx(statistics.pstdev(values))


# ---------------------------------------------------------------------------
# merge_reservoir_values
# ---------------------------------------------------------------------------

def test_merge_reservoir_values_under_cap_keeps_everything():
    result = merge_reservoir_values([[1, 2], [3, 4]], max_records=10)
    assert result == [1, 2, 3, 4]


def test_merge_reservoir_values_over_cap_subsamples():
    population = list(range(1000))
    result = merge_reservoir_values([population], max_records=50, rng=random.Random(0))
    assert len(result) == 50
    assert set(result) <= set(population)


def test_merge_reservoir_values_handles_empty_chunks():
    result = merge_reservoir_values([[], [1, 2], []], max_records=10)
    assert result == [1, 2]


# ---------------------------------------------------------------------------
# Round-trip parsing of vcfstats.py's chunk output formats
# ---------------------------------------------------------------------------

def test_ab_dp_round_trip(tmp_path):
    path = tmp_path / "chunk_ab_dp.txt"
    path.write_text(
        "allele_balance\ts1\t0.4\t0.6\n"
        "allele_balance_moments\ts1\t2\t1.0\t0.52\t0.4\t0.6\n"
        "genotype_depth\ts1\t10\t12\n"
        "genotype_depth_moments\ts1\t2\t22\t244\t10\t12\n"
    )
    samples, values, moments = parse_ab_dp_file(str(path))

    assert samples == ["s1"]
    assert values[("s1", "allele_balance")] == [0.4, 0.6]
    assert values[("s1", "genotype_depth")] == [10.0, 12.0]
    assert moments[("s1", "allele_balance")] == {
        "n": 2, "sum": 1.0, "sumsq": 0.52, "min": 0.4, "max": 0.6,
    }
    assert moments[("s1", "genotype_depth")] == {
        "n": 2, "sum": 22.0, "sumsq": 244.0, "min": 10.0, "max": 12.0,
    }


def test_ab_dp_round_trip_handles_empty_reservoir(tmp_path):
    path = tmp_path / "chunk_ab_dp.txt"
    path.write_text(
        "allele_balance\ts1\t\n"
        "allele_balance_moments\ts1\t0\t0.0\t0.0\tNA\tNA\n"
    )
    _, values, moments = parse_ab_dp_file(str(path))
    assert values[("s1", "allele_balance")] == []
    assert moments[("s1", "allele_balance")]["n"] == 0
    assert moments[("s1", "allele_balance")]["min"] is None


def test_qual_fmiss_maf_dp_round_trip(tmp_path):
    path = tmp_path / "chunk_site.txt"
    path.write_text(
        "QUAL\tF_MISS\tMAF\tDP\n"
        "30.0\t0.0\t0.25\t15\n"
        "NA\t0.5\t0.0\tNA\n"
        "#MOMENTS\n"
        "metric\tn\tsum\tsumsq\tmin\tmax\n"
        "qual\t1\t30.0\t900.0\t30.0\t30.0\n"
        "fmiss\t2\t0.5\t0.25\t0.0\t0.5\n"
        "maf\t2\t0.25\t0.0625\t0.0\t0.25\n"
        "dp\t1\t15\t225\t15\t15\n"
    )
    rows, moments = parse_qual_fmiss_maf_dp_file(str(path))

    assert rows == [(30.0, 0.0, 0.25, 15.0), (None, 0.5, 0.0, None)]
    assert moments["qual"] == {"n": 1, "sum": 30.0, "sumsq": 900.0, "min": 30.0, "max": 30.0}
    assert moments["dp"] == {"n": 1, "sum": 15.0, "sumsq": 225.0, "min": 15.0, "max": 15.0}


# ---------------------------------------------------------------------------
# Combining multiple chunk files
# ---------------------------------------------------------------------------

def test_combine_ab_dp_files_sums_moments_and_concatenates_values(tmp_path):
    chunk1 = tmp_path / "chunk1_ab_dp.txt"
    chunk1.write_text(
        "allele_balance\ts1\t0.4\n"
        "allele_balance_moments\ts1\t1\t0.4\t0.16\t0.4\t0.4\n"
        "genotype_depth\ts1\t10\n"
        "genotype_depth_moments\ts1\t1\t10\t100\t10\t10\n"
    )
    chunk2 = tmp_path / "chunk2_ab_dp.txt"
    chunk2.write_text(
        "allele_balance\ts1\t0.6\n"
        "allele_balance_moments\ts1\t1\t0.6\t0.36\t0.6\t0.6\n"
        "genotype_depth\ts1\t12\n"
        "genotype_depth_moments\ts1\t1\t12\t144\t12\t12\n"
    )

    samples, values, moments = combine_ab_dp_files([str(chunk1), str(chunk2)])

    assert samples == ["s1"]
    assert sorted(values[("s1", "allele_balance")]) == [0.4, 0.6]
    assert moments[("s1", "allele_balance")]["n"] == 2
    assert moments[("s1", "allele_balance")]["sum"] == pytest.approx(1.0)
    assert moments[("s1", "genotype_depth")]["min"] == 10.0
    assert moments[("s1", "genotype_depth")]["max"] == 12.0


def test_combine_ab_dp_files_rejects_mismatched_samples(tmp_path):
    chunk1 = tmp_path / "chunk1_ab_dp.txt"
    chunk1.write_text("allele_balance\ts1\t0.4\nallele_balance_moments\ts1\t1\t0.4\t0.16\t0.4\t0.4\n")
    chunk2 = tmp_path / "chunk2_ab_dp.txt"
    chunk2.write_text("allele_balance\ts2\t0.6\nallele_balance_moments\ts2\t1\t0.6\t0.36\t0.6\t0.6\n")

    with pytest.raises(ValueError):
        combine_ab_dp_files([str(chunk1), str(chunk2)])


def test_combine_gt_stats_files_concatenates_rows_and_merges_moments(tmp_path):
    chunk1 = tmp_path / "chunk1_site.txt"
    chunk1.write_text(
        "QUAL\tF_MISS\tMAF\tDP\n"
        "30.0\t0.0\t0.25\t15\n"
        "#MOMENTS\n"
        "metric\tn\tsum\tsumsq\tmin\tmax\n"
        "qual\t1\t30.0\t900.0\t30.0\t30.0\n"
        "fmiss\t1\t0.0\t0.0\t0.0\t0.0\n"
        "maf\t1\t0.25\t0.0625\t0.25\t0.25\n"
        "dp\t1\t15\t225\t15\t15\n"
    )
    chunk2 = tmp_path / "chunk2_site.txt"
    chunk2.write_text(
        "QUAL\tF_MISS\tMAF\tDP\n"
        "40.0\t0.5\t0.1\t20\n"
        "#MOMENTS\n"
        "metric\tn\tsum\tsumsq\tmin\tmax\n"
        "qual\t1\t40.0\t1600.0\t40.0\t40.0\n"
        "fmiss\t1\t0.5\t0.25\t0.5\t0.5\n"
        "maf\t1\t0.1\t0.01\t0.1\t0.1\n"
        "dp\t1\t20\t400\t20\t20\n"
    )

    rows, moments = combine_gt_stats_files([str(chunk1), str(chunk2)])

    assert rows == [(30.0, 0.0, 0.25, 15.0), (40.0, 0.5, 0.1, 20.0)]
    assert moments["qual"]["n"] == 2
    assert moments["qual"]["sum"] == pytest.approx(70.0)
    assert moments["dp"]["max"] == 20.0
