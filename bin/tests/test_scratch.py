import gzip
import os
from pathlib import Path
import re
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[2]
HELPER = ROOT / "bin" / "task_tmp.sh"


def run_shell(work, script, enabled=True, scratch=None):
    return subprocess.run(
        ["bash", "-uec", 'source "$1" "$2" "$3"\n' + script, "test",
         str(HELPER), str(enabled).lower(), str(scratch or "")],
        cwd=work, text=True, capture_output=True,
    )


@pytest.mark.parametrize("enabled", [False, True])
@pytest.mark.parametrize("status", [0, 17])
def test_temporary_workspace_is_cleaned_without_moving_outputs(tmp_path, enabled, status):
    work = tmp_path / "shared work"
    scratch = tmp_path / "local scratch"
    work.mkdir()
    scratch.mkdir()
    sentinel = scratch / "unrelated"
    sentinel.write_text("keep")
    result = run_shell(
        work,
        'printf "%s\\n" "$PWD" "$TMPDIR" > locations.txt\n'
        'printf temporary > "$TMPDIR/private.bam"\n'
        f"printf output > final.bam\nexit {status}",
        enabled, scratch,
    )
    assert result.returncode == status, result.stderr
    cwd, temporary = (work / "locations.txt").read_text().splitlines()
    assert Path(cwd) == work
    assert Path(temporary).parent == (scratch if enabled else work)
    assert not Path(temporary).exists()
    assert (work / "final.bam").read_text() == "output"
    assert sentinel.read_text() == "keep"
    assert {p.name for p in work.iterdir()} == {"locations.txt", "final.bam"}


@pytest.mark.parametrize("scratch", ["", "missing"])
def test_missing_scratch_fails_without_fallback(tmp_path, scratch):
    result = run_shell(
        tmp_path, "touch should-not-exist", scratch=tmp_path / scratch if scratch else None,
    )
    assert result.returncode != 0
    assert "scratch" in result.stderr or "temporary directory root" in result.stderr
    assert not (tmp_path / "should-not-exist").exists()
    assert not list(tmp_path.glob("map-and-call.*"))


def test_pipeline_failure_is_not_hidden_by_compression_or_cleanup(tmp_path):
    result = run_shell(tmp_path, "false | gzip -c > output.gz", enabled=False)
    assert result.returncode != 0
    assert not list(tmp_path.glob("map-and-call.*"))


def test_termination_cleans_temporary_workspace(tmp_path):
    result = run_shell(
        tmp_path, 'printf "%s" "$TMPDIR" > location.txt\nkill -TERM $$',
        enabled=False,
    )
    assert result.returncode == 143
    assert not Path((tmp_path / "location.txt").read_text()).exists()


def test_modules_never_enable_whole_task_scratch():
    for module in (ROOT / "modules").glob("*/*.nf"):
        assert not any(line.strip().startswith("scratch ") for line in module.read_text().splitlines()), module


def cached_image(pattern):
    matches = sorted((ROOT / ".singularity").glob(pattern))
    if not matches or not shutil.which("singularity"):
        pytest.skip("Integration test requires cached Singularity images")
    return matches[0]


def execute(command, cwd):
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def test_all_modules_compile(tmp_path):
    if not shutil.which("nextflow"):
        pytest.skip("Compilation test requires Nextflow on PATH")
    imports = []
    for module in sorted((ROOT / "modules").glob("*/*.nf")):
        names = re.findall(r"^process\s+(\w+)", module.read_text(), re.MULTILINE)
        for name in names:
            imports.append(f"include {{ {name} as checked_{len(imports)} }} from '{module}'")
    (tmp_path / "compile.nf").write_text("\n".join(imports) + "\nworkflow {}\n")
    (tmp_path / "empty.config").write_text("")
    execute(["nextflow", "-C", "empty.config", "run", "compile.nf", "-ansi-log", "false"], tmp_path)


@pytest.mark.parametrize("enabled", [False, True])
def test_download_modules_keep_only_compressed_outputs(tmp_path, enabled):
    if not all(shutil.which(tool) for tool in ["nextflow", "pigz", "unzip"]):
        pytest.skip("Download routing test requires Nextflow, pigz and unzip")
    image = cached_image("*library-bcftools_bedtools_pysam_python-*.img")
    scratch = tmp_path / "scratch space"
    scratch.mkdir()
    mock_bin = tmp_path / "mock_bin"
    mock_bin.mkdir()
    provider = mock_bin / "provider"
    provider.write_text("""#!/usr/bin/env python3
import os
from pathlib import Path
import sys
import zipfile

name = Path(sys.argv[0]).name
args = sys.argv[1:]
temporary = Path(os.environ["TMPDIR"])
if name == "prefetch":
    output = Path(args[args.index("-O") + 1])
    assert output == temporary
    (output / args[0]).mkdir()
    (output / args[0] / (args[0] + ".sra")).write_text("archive")
elif name == "fasterq-dump":
    output = Path(args[args.index("-O") + 1])
    assert output == temporary
    assert Path(args[args.index("-t") + 1]) == temporary
    assert Path(args[0]).parent == temporary
    accession = Path(args[0]).name
    for mate in [1, 2]:
        (output / f"{accession}_{mate}.fastq").write_text("@read\\nACGT\\n+\\nIIII\\n")
elif name == "datasets":
    output = Path(args[args.index("--filename") + 1])
    assert output.parent == temporary
    with zipfile.ZipFile(output, "w") as archive:
        archive.writestr("data/genome_genomic.fna", ">chr1\\nACGT\\n")
else:
    raise ValueError(name)
""")
    provider.chmod(0o755)
    for tool in ["prefetch", "fasterq-dump", "datasets"]:
        (mock_bin / tool).symlink_to(provider)
    bgzip = mock_bin / "bgzip"
    bgzip.write_text(
        f'#!/bin/bash\nexec singularity exec --bind "{tmp_path}" "{image}" bgzip "$@"\n'
    )
    bgzip.chmod(0o755)
    (tmp_path / "downloads.nf").write_text(f"""
include {{ sratools_fastq_dump }} from '{ROOT}/modules/sratools/sratools_fastq_dump'
include {{ datasets_download }} from '{ROOT}/modules/ncbi_datasets/ncbi_datasets'
workflow {{
    sratools_fastq_dump(tuple('sample', 'lane1', 1, 'lib', 'TEST123'))
    datasets_download('TEST_ASSEMBLY')
}}
""")
    (tmp_path / "test.config").write_text(f"""
includeConfig '{ROOT}/nextflow.config'
plugins = []
conda.enabled = false
trace.enabled = false
report.enabled = false
timeline.enabled = false
params.use_scratch = {str(enabled).lower()}
params.scratch_dir = '$SCRATCH_TEST_ROOT'
process.executor = 'local'
process.cpus = 1
env.PATH = '{mock_bin}:{ROOT}/bin:' + System.getenv('PATH')
""")
    env = dict(os.environ, SCRATCH_TEST_ROOT=str(scratch), NXF_OFFLINE="true")
    result = subprocess.run(
        ["nextflow", "-C", "test.config", "run", "downloads.nf", "-ansi-log", "false"],
        cwd=tmp_path, env=env, text=True, capture_output=True, timeout=180,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert not list(scratch.iterdir())
    outputs = {}
    for command in (tmp_path / "work").glob("*/*/.command.sh"):
        assert not list(command.parent.glob("map-and-call.*"))
        for output in command.parent.glob("*"):
            if not output.name.startswith("."):
                assert output.is_file() and output.name.endswith(".gz"), output
                with gzip.open(output, "rt") as compressed:
                    outputs[output.name] = compressed.read()
    assert outputs == {
        "TEST123_1.fastq.gz": "@read\nACGT\n+\nIIII\n",
        "TEST123_2.fastq.gz": "@read\nACGT\n+\nIIII\n",
        "TEST_ASSEMBLY_genomic.fna.gz": ">chr1\nACGT\n",
    }


@pytest.mark.parametrize("enabled", [False, True])
def test_nextflow_temporary_data_and_streamed_outputs(tmp_path, enabled):
    if not shutil.which("nextflow"):
        pytest.skip("Integration test requires Nextflow on PATH")
    samtools_image = cached_image("*library-bedtools_samtools_gzip-*.img")
    filter_image = cached_image("*library-bcftools_bedtools_pysam_python-*.img")
    minibwa_image = cached_image("*library-minibwa_samtools-*.img")
    bbmap_image = cached_image("*library-bbmap-*.img")
    fastp_image = cached_image("*library-fastp-*.img")
    scratch = tmp_path / "scratch space"
    scratch.mkdir()
    (scratch / "unrelated").write_text("keep")
    (tmp_path / "reference.fa").write_text(">chr1\n" + "A" * 100 + "\n")
    (tmp_path / "reads.sam").write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n"
        "@RG\tID:rg1\tSM:s1\tLB:lib1\n"
        "read1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rg1\n"
    )
    runtime = ["singularity", "exec", "--bind", str(tmp_path), str(samtools_image)]
    execute(runtime + ["samtools", "faidx", "reference.fa"], tmp_path)
    execute(runtime + ["samtools", "view", "-b", "-o", "reads.bam", "reads.sam"], tmp_path)
    execute(runtime + ["samtools", "index", "reads.bam"], tmp_path)
    vcf_header = (
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
    )
    for name, row in [
        ("snps", "chr1\t5\t.\tA\tC\t60\tPASS\t.\tGT:AD\t0/1:5,5\t0/0:10,0\n"),
        ("indels", "chr1\t8\t.\tAA\tA\t60\tPASS\t.\tGT:AD\t0/1:5,5\t0/0:10,0\n"),
    ]:
        (tmp_path / f"{name}.vcf").write_text(vcf_header + row)
        runtime = ["singularity", "exec", "--bind", str(tmp_path), str(filter_image)]
        execute(runtime + ["bcftools", "view", "-Oz", "-o", f"{name}.vcf.gz", f"{name}.vcf"], tmp_path)
        execute(runtime + ["bcftools", "index", f"{name}.vcf.gz"], tmp_path)
    for sample in ["s1", "s2"]:
        with gzip.open(tmp_path / f"{sample}.bed.gz", "wt") as bed:
            bed.write("chr1\t0\t20\n")
    for region, content in [(2, "chr1\t10\t20\t4\n"), (1, "chr1\t0\t10\t2\n")]:
        with gzip.open(tmp_path / f"{region}_s1.depths.bed.gz", "wt") as bed:
            bed.write(content)
    (tmp_path / "first.bed").write_text("chr1\t0\t10\n")
    (tmp_path / "second.bed").write_text("chr1\t10\t20\n")
    sequence = "ACGTTGCAAGTCGATCGTAGCTAGTCAGTACGATCGTACGATGCTAGCTAGTCGATCGTA"
    for mate in [1, 2]:
        with gzip.open(tmp_path / f"reads{mate}.fastq.gz", "wt") as reads:
            reads.write(f"@read1/{mate}\n{sequence}\n+\n" + "I" * len(sequence) + "\n")

    workflow = r"""
include { samtools_dp; samtools_dp_deprecated } from '__ROOT__/modules/samtools/samtools_dp_process'
include { parse_region_depths } from '__ROOT__/modules/samtools/parse_region_depths'
include { samtools_markdups } from '__ROOT__/modules/samtools/samtools_markdups'
include { callability_filter; callability_filter_legacy } from '__ROOT__/modules/bedtools/callability_filter'
include { finalize_masks; finalize_masks_legacy } from '__ROOT__/modules/bedtools/finalize_masks'
include { callable_regions } from '__ROOT__/modules/bedtools/callable_regions'
include { combine_bedfiles } from '__ROOT__/modules/bedtools/combine_bedfiles'
include { ab_dp_filter } from '__ROOT__/modules/custom_variant_filters/ab_dp_filter'
include { minibwa_index } from '__ROOT__/modules/minibwa/minibwa_index'
include { minibwa_map } from '__ROOT__/modules/minibwa/minibwa_map'
include { clumpify_single; clumpify_paired } from '__ROOT__/modules/bbmap/bbmap'
include { fastp } from '__ROOT__/modules/fastp/trimming'

workflow {
    def bam = file('reads.bam')
    def bai = file('reads.bam.bai')
    def fai = file('reference.fa.fai')
    def snps = file('snps.vcf.gz')
    def snpidx = file('snps.vcf.gz.csi')
    def indels = file('indels.vcf.gz')
    def indelidx = file('indels.vcf.gz.csi')
    def beds = [file('s1.bed.gz'), file('s2.bed.gz')]
    def depths = [file('2_s1.depths.bed.gz'), file('1_s1.depths.bed.gz')]
    samtools_dp(tuple(1, ['chr1:1-20'], ['s1'], [bam], [bai]))
    samtools_dp_deprecated(tuple('old', bam, bai, 1, ['chr1:1-20']))
    parse_region_depths(tuple('s1', depths), fai)
    samtools_markdups(tuple('s1', 1, bam, bai))
    callability_filter(tuple(1, ['chr1:1-20'], snps, snpidx, ['s1', 's2'], beds, 'snps'))
    callability_filter_legacy(tuple(1, ['chr1:1-20'], snps, snpidx, 's1', beds[0]))
    finalize_masks(tuple(['s1', 's2'], 1, ['chr1:1-20'], beds, snps, snpidx, indels, indelidx))
    finalize_masks_legacy(tuple('old', 1, ['chr1:1-20'], beds[0], snps, indels))
    callable_regions(tuple('s1', 1, 10, 'homozygous', depths[1], [], [], fai))
    combine_bedfiles(tuple('s1', [file('first.bed'), file('second.bed')], fai, 'combined'))
    ab_dp_filter(tuple('s1', 1, 'chr1:1-20', beds[0], snps, snpidx, 1, 20, 'homozygous'), '', 'snps')
    def r1 = file('reads1.fastq.gz')
    def r2 = file('reads2.fastq.gz')
    minibwa_index(file('reference.fa'))
    minibwa_map(minibwa_index.out.reference.map { ref, idx -> tuple('s1', 'lib', 1, r1, r2, ref, idx) })
    clumpify_single(tuple('s1', 'single', 1, r1))
    clumpify_paired(tuple('s1', 'paired', 1, r1, r2))
    fastp(tuple('s1', 'lane1', 1, 'lib', r1, r2))
}
"""
    (tmp_path / "smoke.nf").write_text(workflow.replace("__ROOT__", str(ROOT)))
    config = f"""
includeConfig '{ROOT}/nextflow.config'
plugins = []
conda.enabled = false
singularity.enabled = true
singularity.autoMounts = true
singularity.cacheDir = '{ROOT}/.singularity'
singularity.runOptions = '--bind {ROOT}'
trace.enabled = false
report.enabled = false
timeline.enabled = false
params.with_singularity = true
params.min_mapqual = 0
params.min_basequal = 0
params.clumpify_threads = 1
params.use_scratch = {str(enabled).lower()}
params.scratch_dir = '$SCRATCH_TEST_ROOT'
params.outdir = '{tmp_path}/reports'
process.executor = 'local'
process.cpus = 1
process.memory = '1 GB'
process.container = '{filter_image}'
process {{
    withName: 'minibwa_index|minibwa_map' {{
        container = '{minibwa_image}'
        cpus = 1
        memory = '1 GB'
    }}
    withName: 'clumpify_single|clumpify_paired' {{
        container = '{bbmap_image}'
        cpus = 1
        memory = '4 GB'
    }}
    withName: 'fastp' {{
        container = '{fastp_image}'
        cpus = 1
        memory = '1 GB'
    }}
    withName: 'callability_filter.*|finalize_masks.*|callable_regions|combine_bedfiles|ab_dp_filter' {{
        container = '{filter_image}'
    }}
    withName: 'samtools_.*|parse_region_depths' {{
        container = '{samtools_image}'
        cpus = 1
        memory = '1 GB'
    }}
}}
env.PATH = '{ROOT}/bin:/opt/conda/bin:/usr/bin:/bin'
"""
    (tmp_path / "test.config").write_text(config)
    env = dict(os.environ, SCRATCH_TEST_ROOT=str(scratch), NXF_OFFLINE="true")
    result = subprocess.run(
        ["nextflow", "-C", "test.config", "run", "smoke.nf", "-ansi-log", "false"],
        cwd=tmp_path, env=env, text=True, capture_output=True, timeout=180,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert {p.name for p in scratch.iterdir()} == {"unrelated"}

    tasks = {}
    for command in (tmp_path / "work").glob("*/*/.command.sh"):
        script = command.read_text()
        wrapper = command.with_name(".command.run").read_text()
        assert "NXF_SCRATCH=''" in wrapper
        assert not list(command.parent.glob("map-and-call.*"))
        if "source task_tmp.sh" in script:
            log = command.with_name(".command.err").read_text()
            temporary = next(line.removeprefix("Temporary workspace: ")
                             for line in log.splitlines() if line.startswith("Temporary workspace: "))
            assert Path(temporary).parent == (scratch if enabled else command.parent)
            assert not Path(temporary).exists()
        for line in wrapper.splitlines():
            if line.startswith("### name:"):
                name = line.split("'")[1].split(" ")[0]
                tasks[name] = command.parent
        if enabled:
            assert '--bind "$SCRATCH_TEST_ROOT"' in wrapper
            assert '--env "TASK_SCRATCH_ROOT=$SCRATCH_TEST_ROOT"' in wrapper
    assert len(tasks) == 16

    def read_gzip(task, name):
        with gzip.open(tasks[task] / name, "rt") as output:
            return output.read()

    assert read_gzip("parse_region_depths", "s1.depths.bed.gz") == (
        "chr1\t0\t10\t2\nchr1\t10\t20\t4\n"
    )
    assert (tasks["parse_region_depths"] / "s1.depths.avg.txt").read_text() == "chr1\t20\t3\n"
    assert not list(tasks["parse_region_depths"].glob("*.depths.bed"))
    assert read_gzip("samtools_dp", "1_s1.depths.bed.gz") == (
        "chr1\t0\t10\t1\nchr1\t10\t20\t0\n"
    )
    assert not list(tasks["samtools_dp"].glob("*.depths.bed"))
    assert read_gzip("callable_regions", "s1.callable_regions.bed.gz") == "chr1\t0\t10\n"
    for sample in ["s1", "s2"]:
        assert read_gzip("finalize_masks", f"{sample}_1.mappability_mask_total.bed.gz") == "chr1\t0\t20\n"
    assert read_gzip("finalize_masks", "s1_1.homref_invariants.bed.gz") == (
        "chr1\t0\t4\nchr1\t5\t7\nchr1\t9\t20\n"
    )
    assert (tasks["samtools_markdups"] / "s1.dedup.bam.bai").is_file()
    assert (tasks["combine_bedfiles"] / "s1_combined.bed").read_text() == (
        "chr1\t0\t10\nchr1\t10\t20\n"
    )
    for name in ["callability_filter", "callability_filter_legacy", "ab_dp_filter"]:
        assert len([p for p in tasks[name].glob("*.vcf.gz") if not p.is_symlink()]) == 1
        assert len([p for p in tasks[name].glob("*.vcf.gz.csi") if not p.is_symlink()]) == 1
    assert (tasks["minibwa_map"] / "s1_lib.sorted.bam.bai").is_file()
    assert read_gzip("clumpify_single", "s1_single_dedup_fastq.gz").splitlines()[1] == sequence
    assert read_gzip("clumpify_paired", "s1_paired_dedup_R1.fastq.gz").splitlines()[1] == sequence
    assert read_gzip("fastp", "s1_lane1_trimmed_R1.fastq.gz").splitlines()[1] == sequence
    assert "Temporary workspace:" not in (tasks["fastp"] / ".command.err").read_text()
