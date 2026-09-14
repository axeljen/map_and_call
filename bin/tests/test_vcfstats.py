import random
import statistics
from dataclasses import dataclass

import pysam
import pytest

from vcfstats import (
    ReservoirSampler,
    RunningMoments,
    StatsAccumulator,
    classify_genotype,
    classify_record,
    compute_allele_balance,
    safe_ad_sum,
)


# ---------------------------------------------------------------------------
# classify_record
# ---------------------------------------------------------------------------

@dataclass
class FakeRecord:
    """Minimal stand-in for a pysam VariantRecord -- classify_record only
    ever looks at .ref/.alts, so a real pysam object isn't needed here."""
    ref: str
    alts: tuple


@pytest.mark.parametrize("ref, alts, expected", [
    ("A", ("T",), "snp_biallelic"),
    ("A", ("T", "C"), "snp_multiallelic"),
    ("AT", ("A",), "indel_biallelic"),
    ("A", ("AT", "ATG"), "indel_multiallelic"),
    ("A", None, "no_alt"),
    ("A", (), "no_alt"),  # empty (not None) alts -- previously misclassified
])
def test_classify_record(ref, alts, expected):
    assert classify_record(FakeRecord(ref, alts)) == expected


# ---------------------------------------------------------------------------
# classify_genotype
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("gt, expected", [
    ((0, 0), "hom_ref"),
    ((1, 1), "hom_alt"),
    ((2, 2), "hom_alt"),
    ((0, 1), "het"),
    ((1, 0), "het"),
    ((None, None), "missing"),
    ((0, None), "missing"),
    (None, "missing"),
])
def test_classify_genotype(gt, expected):
    assert classify_genotype(gt) == expected


# ---------------------------------------------------------------------------
# safe_ad_sum / compute_allele_balance
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("ad, expected", [
    ((5, 10), 15),
    ((5, None), 5),
    ((None, None), None),
    (None, None),
    ((), None),
])
def test_safe_ad_sum(ad, expected):
    assert safe_ad_sum(ad) == expected


@pytest.mark.parametrize("ad, expected", [
    ((5, 10), round(5 / 15, 3)),
    ((0, 0), None),  # zero total depth
    ((5, None), None),  # can't compute without both alleles
    (None, None),
    ((5,), None),  # fewer than two alleles
])
def test_compute_allele_balance(ad, expected):
    assert compute_allele_balance(ad) == expected


# ---------------------------------------------------------------------------
# ReservoirSampler
# ---------------------------------------------------------------------------

def test_reservoir_sampler_respects_capacity():
    sampler = ReservoirSampler(capacity=5, rng=random.Random(0))
    for i in range(1000):
        sampler.add(i)
    assert len(sampler.values) == 5
    assert sampler.n_seen == 1000


def test_reservoir_sampler_deterministic_with_seed():
    a = ReservoirSampler(capacity=10, rng=random.Random(42))
    b = ReservoirSampler(capacity=10, rng=random.Random(42))
    for i in range(500):
        a.add(i)
        b.add(i)
    assert a.values == b.values


def test_reservoir_sampler_under_capacity_keeps_everything():
    sampler = ReservoirSampler(capacity=100, rng=random.Random(0))
    for i in range(10):
        sampler.add(i)
    assert sorted(sampler.values) == list(range(10))


def test_reservoir_sampler_sample_mean_converges():
    rng = random.Random(1234)
    population = list(range(100000))
    sampler = ReservoirSampler(capacity=5000, rng=rng)
    for v in population:
        sampler.add(v)
    sample_mean = statistics.mean(sampler.values)
    population_mean = statistics.mean(population)
    assert abs(sample_mean - population_mean) < 0.05 * population_mean


# ---------------------------------------------------------------------------
# RunningMoments
# ---------------------------------------------------------------------------

def test_running_moments_matches_direct_computation():
    values = [1.0, 5.0, 2.5, 9.0, -3.0, 4.25]
    moments = RunningMoments()
    for v in values:
        moments.add(v)
    assert moments.n == len(values)
    assert moments.mean == pytest.approx(statistics.mean(values))
    assert moments.sd == pytest.approx(statistics.pstdev(values))
    assert moments.min == min(values)
    assert moments.max == max(values)


def test_running_moments_merge_matches_concatenation():
    values = [1.0, 5.0, 2.5, 9.0, -3.0, 4.25, 7.0]
    first_half, second_half = values[:3], values[3:]

    a, b, combined = RunningMoments(), RunningMoments(), RunningMoments()
    for v in first_half:
        a.add(v)
    for v in second_half:
        b.add(v)
    for v in values:
        combined.add(v)

    merged = a.merge(b)
    assert merged.n == combined.n
    assert merged.sum == pytest.approx(combined.sum)
    assert merged.sumsq == pytest.approx(combined.sumsq)
    assert merged.min == combined.min
    assert merged.max == combined.max
    assert merged.mean == pytest.approx(combined.mean)
    assert merged.sd == pytest.approx(combined.sd)


def test_running_moments_merge_with_empty():
    a = RunningMoments()
    b = RunningMoments()
    b.add(3.0)
    assert a.merge(b) == RunningMoments(n=1, sum=3.0, sumsq=9.0, min=3.0, max=3.0)


# ---------------------------------------------------------------------------
# StatsAccumulator.process_record, end-to-end over a tiny in-memory VCF
# ---------------------------------------------------------------------------

def _make_header():
    header = pysam.VariantHeader()
    header.add_line("##contig=<ID=chr1,length=1000>")
    header.add_line('##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">')
    for sample in ("s1", "s2"):
        header.add_sample(sample)
    return header


@pytest.fixture
def three_record_vcf(tmp_path):
    """A tiny VCF: one biallelic SNP, one biallelic indel, one multiallelic
    (non-biallelic) SNP -- exercising both the biallelic and "other" code
    paths of StatsAccumulator.process_record."""
    header = _make_header()
    path = tmp_path / "test.vcf"

    with pysam.VariantFile(str(path), "w", header=header) as vcf_out:
        rec1 = vcf_out.new_record(contig="chr1", start=0, alleles=("A", "T"), qual=30)
        rec1.info["DP"] = 15
        rec1.samples["s1"]["GT"] = (0, 0)
        rec1.samples["s1"]["AD"] = (10, 0)
        rec1.samples["s2"]["GT"] = (0, 1)
        rec1.samples["s2"]["AD"] = (5, 5)
        vcf_out.write(rec1)

        rec2 = vcf_out.new_record(contig="chr1", start=1, alleles=("AT", "A"))
        rec2.samples["s1"]["GT"] = (1, 1)
        rec2.samples["s1"]["AD"] = (0, 8)
        rec2.samples["s2"]["GT"] = (None, None)
        vcf_out.write(rec2)

        rec3 = vcf_out.new_record(contig="chr1", start=2, alleles=("A", "C", "G"))
        rec3.samples["s1"]["GT"] = (0, 1)
        rec3.samples["s2"]["GT"] = (2, 2)
        vcf_out.write(rec3)

    with pysam.VariantFile(str(path)) as vcf_in:
        yield list(vcf_in.fetch())


def test_process_record_end_to_end(three_record_vcf):
    acc = StatsAccumulator(["s1", "s2"], max_sampled_values=100, rng=random.Random(0))
    for record in three_record_vcf:
        acc.process_record(record)

    assert acc.rec_type_counts == {
        "snp_biallelic": 1,
        "snp_multiallelic": 1,
        "indel_biallelic": 1,
        "indel_multiallelic": 0,
        "no_alt": 0,
    }

    s1 = acc.sample_counts["s1"]
    assert (s1.num_records, s1.num_hom_ref, s1.num_het, s1.num_hom_alt, s1.num_missing) \
        == (3, 1, 1, 1, 0)

    s2 = acc.sample_counts["s2"]
    assert (s2.num_records, s2.num_hom_ref, s2.num_het, s2.num_hom_alt, s2.num_missing) \
        == (3, 0, 1, 1, 1)

    # genotype depth only accrues for biallelic records with resolvable AD
    assert acc.sample_moments["s1"]["genotype_depth"].n == 2
    assert acc.sample_moments["s1"]["genotype_depth"].sum == pytest.approx(18)
    assert acc.sample_moments["s2"]["genotype_depth"].n == 1  # rec2 missing, rec3 not biallelic

    # allele balance only accrues for het calls on biallelic records
    assert acc.sample_moments["s2"]["allele_balance"].n == 1
    assert acc.sample_moments["s2"]["allele_balance"].sum == pytest.approx(0.5)
    assert acc.sample_moments["s1"]["allele_balance"].n == 0  # s1 never het on a biallelic site

    # site-level: only the 2 biallelic records contribute
    assert acc.site_moments["qual"].n == 1  # rec2 has no QUAL
    assert acc.site_moments["qual"].sum == pytest.approx(30)
    assert acc.site_moments["dp"].n == 1  # rec2 has no INFO/DP
    assert len(acc.site_reservoir.values) == 2

    rows = acc.site_reservoir.values
    rec1_row = next(r for r in rows if r[0] == 30.0)
    assert rec1_row == (30.0, 0.0, 0.25, 15)  # fmiss=0 (both called), maf=min(3,1)/4

    rec2_row = next(r for r in rows if r[0] is None)
    assert rec2_row == (None, 0.5, 0.0, None)  # 1/2 samples called, maf=min(0,2)/2


def test_process_record_missing_sample_raises(three_record_vcf):
    acc = StatsAccumulator(["s1", "s2", "not_in_vcf"], rng=random.Random(0))
    with pytest.raises(ValueError, match="not_in_vcf"):
        acc.process_record(three_record_vcf[0])
