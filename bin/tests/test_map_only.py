import re
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize(
    "mode,skip,variant_calling,config_mode,downsample",
    [
        ("map_only", "false", False, False, False),
        ("map_only", "false", False, True, True),
        ("map_and_call", "true", False, False, False),
        ("map_and_call", "false", True, False, False),
    ],
)
def test_mapping_mode_workflow_graph(
    tmp_path, mode, skip, variant_calling, config_mode, downsample
):
    if not shutil.which("nextflow"):
        pytest.skip("Workflow test requires Nextflow on PATH")

    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGT\n")
    sheet = tmp_path / "samples.csv"
    sheet.write_text("sample_id,data_type,library,read_1,read_2\n")
    config = tmp_path / "test.config"
    config.write_text(
        "conda.enabled = false\n"
        + (f"params.mode = '{mode}'\n" if config_mode else "")
    )
    mode_args = [] if config_mode else ["--mode", mode]
    downsample_args = ["--downsample_bams", "true"] if downsample else []
    result = subprocess.run(
        [
            "nextflow", "-c", str(config), "run", str(ROOT / "main.nf"),
            "-stub-run", "-ansi-log", "false", "-with-dag", str(tmp_path / "dag.dot"),
            "--outdir", str(tmp_path / "output"),
            *mode_args, "--skip_variant_calling", skip, *downsample_args,
            "--input", str(sheet),
            "--reference", str(reference),
        ],
        cwd=tmp_path, text=True, capture_output=True, timeout=180,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    log = (tmp_path / ".nextflow.log").read_text()
    processes = set(re.findall(r"Creating process '([^']+)'", log))
    assert "PREPROCESS_MODERN:fastp" in processes
    assert "PROCESS_BAMS:samtools_merge" in processes
    assert "INDEX_REFERENCE:bwa_index" in processes
    assert any(p.startswith("MAP_MODERN:") for p in processes)
    assert any(p.startswith("VARIANT_CALLING:") for p in processes) == variant_calling
    assert any(p.startswith("VARIANT_FILTERS:") for p in processes) == variant_calling
    assert ("PROCESS_BAMS:samtools_downsample" in processes) == downsample
    if not variant_calling:
        assert not (tmp_path / "output" / "03_genotypes").exists()


@pytest.mark.parametrize(
    "mode,mapper,historical,expected_indexes",
    [
        ("map_only", "bwa_mem", False, {"bwa_index"}),
        ("map_only", "minibwa", False, {"minibwa_index"}),
        ("map_only", "minibwa", True, {"bwa_index", "minibwa_index"}),
        ("map_and_call", "minibwa", True, {"bwa_index", "minibwa_index"}),
        ("call_variants", "bwa_mem", False, set()),
    ],
)
def test_reference_indexes_for_mapping_modes(tmp_path, mode, mapper, historical, expected_indexes):
    if not shutil.which("nextflow"):
        pytest.skip("Workflow test requires Nextflow on PATH")

    script = tmp_path / "index.nf"
    script.write_text(
        f"include {{ INDEX_REFERENCE }} from '{ROOT / 'subworkflows/index_reference'}'\n"
        "workflow {\n"
        f"    INDEX_REFERENCE(channel.empty(), null, 20, [], [], [], [], {str(historical).lower()})\n"
        "}\n"
    )
    config = tmp_path / "test.config"
    config.write_text(
        f"params.mode = '{mode}'\n"
        f"params.modern_mapper = '{mapper}'\n"
        "params.historical_mapper = 'bwa_mem'\n"
        "conda.enabled = false\n"
    )
    result = subprocess.run(
        ["nextflow", "-C", str(config), "run", str(script), "-ansi-log", "false"],
        cwd=tmp_path, text=True, capture_output=True, timeout=90,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    log = (tmp_path / ".nextflow.log").read_text()
    processes = set(re.findall(r"Creating process '([^']+)'", log))
    for index in ("bwa_index", "minibwa_index"):
        assert (f"INDEX_REFERENCE:{index}" in processes) == (index in expected_indexes)
