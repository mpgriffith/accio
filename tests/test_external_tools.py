"""Tests for external tool wrappers in accio.wrappers.external_tools."""

import pytest
import subprocess
from pathlib import Path
from unittest.mock import patch, MagicMock

from accio.wrappers.external_tools import BLASTRunner, ExternalToolError, check_tool_availability


def test_makeblastdb_error_handling(monkeypatch):
    """Test that BLASTRunner.run_makeblastdb raises ExternalToolError on failure."""
    runner = BLASTRunner()

    def mock_run(*args, **kwargs):
        raise subprocess.CalledProcessError(
            returncode=1,
            cmd=['makeblastdb'],
            stderr=b"makeblastdb error: invalid fasta file"
        )

    monkeypatch.setattr("subprocess.run", mock_run)

    with pytest.raises(ExternalToolError) as exc_info:
        runner.run_makeblastdb("dummy.fasta", "dummy_out")
    assert "makeblastdb error: invalid fasta file" in str(exc_info.value)


def test_check_tool_availability_bwa_and_new_tools(monkeypatch):
    """Test that check_tool_availability correctly detects bwa and new database tools."""
    def mock_which(cmd):
        # Pretend all tools are present on system PATH
        return f"/usr/bin/{cmd}"

    def mock_run(cmd, *args, **kwargs):
        binary = cmd[0]
        # bwa prints usage to stderr and exits with 1
        if binary == 'bwa':
            return subprocess.CompletedProcess(
                cmd, returncode=1, stdout=b"", stderr=b"Program: bwa (alignment via Burrows-Wheeler transformation)\n"
            )
        return subprocess.CompletedProcess(cmd, returncode=0, stdout=b"v1.0", stderr=b"")

    monkeypatch.setattr("shutil.which", mock_which)
    monkeypatch.setattr("subprocess.run", mock_run)

    avail = check_tool_availability()

    # bwa should be available even though it exits with 1
    assert avail['bwa']['available'] is True

    # Check newly added tools from issue 12
    for tool in ['skani', 'amrfinder', 'pling', 'plasnet', 'mob_typer', 'mob_cluster']:
        assert tool in avail, f"{tool} should be in checked tools"
        assert avail[tool]['available'] is True


def test_validate_databases(tmp_path):
    """Test that validate_databases checks the actual built DB files and not obsolete ResFinder."""
    from accio.wrappers.external_tools import validate_databases

    # Empty dir should fail all
    status_empty = validate_databases(str(tmp_path))
    assert not any(status_empty.values())
    assert 'resfinder_db' not in status_empty, "ResFinder should not be checked"

    # Create all required DB files
    (tmp_path / 'plasmidDB_info.csv').write_text("info")
    (tmp_path / 'plasmidDB.fasta').write_text(">p\nATGC\n")
    (tmp_path / 'plasmidDB.nhr').write_text("nhr")
    (tmp_path / 'plasmidDB.fasta.msh').write_text("msh")
    (tmp_path / 'PlasmidFinder.nhr').write_text("pf_nhr")
    (tmp_path / 'repetitive.dna.fas.nhr').write_text("rep_nhr")

    status_full = validate_databases(str(tmp_path))
    assert all(status_full.values()), f"All database files should exist, got: {status_full}"
    assert status_full['plasmidfinder_blast'] is True
    assert status_full['rep_blast'] is True


def test_mob_runner_tool_resolution_path(monkeypatch):
    """Test MOBRunner resolves directly if tool is on PATH."""
    from accio.wrappers.external_tools import MOBRunner
    runner = MOBRunner()

    monkeypatch.setattr("shutil.which", lambda cmd: f"/bin/{cmd}" if cmd == "mob_recon" else None)
    resolved = runner._resolve_tool("mob_recon")
    assert resolved == ["mob_recon"]


def test_mob_runner_tool_resolution_conda(monkeypatch):
    """Test MOBRunner resolves via conda/mamba if binary is not on PATH."""
    from accio.wrappers.external_tools import MOBRunner
    from accio.config import AnalysisConfig

    cfg = AnalysisConfig()
    cfg.MOB_CONDA_ENV = "my_mob_env"
    runner = MOBRunner(config=cfg)

    def mock_which(cmd):
        if cmd == "conda":
            return "/usr/bin/conda"
        return None

    monkeypatch.setattr("shutil.which", mock_which)
    resolved = runner._resolve_tool("mob_recon")
    assert resolved == ["/usr/bin/conda", "run", "--no-capture-output", "-n", "my_mob_env", "mob_recon"]


def test_mob_runner_tool_resolution_custom_dir(tmp_path):
    """Test MOBRunner resolves via MOB_BIN_DIR."""
    from accio.wrappers.external_tools import MOBRunner
    from accio.config import AnalysisConfig

    fake_bin = tmp_path / "mob_recon"
    fake_bin.write_text("#!/bin/sh\n")

    cfg = AnalysisConfig()
    cfg.MOB_BIN_DIR = str(tmp_path)
    runner = MOBRunner(config=cfg)

    resolved = runner._resolve_tool("mob_recon")
    assert resolved == [str(fake_bin)]


def test_mob_runner_error_format():
    """Test MOBRunner produces clear diagnostic message mentioning conda environment."""
    from accio.wrappers.external_tools import MOBRunner
    runner = MOBRunner()

    msg = runner._format_mob_error("MOB-recon", Exception("Command not found"))
    assert "conda create -n mob_suite -c bioconda mob_suite" in msg
    assert "pandas <= 1.5.3" in msg


def test_blastn_handles_paths_with_spaces_and_threads(monkeypatch):
    """Test that BLASTRunner.run_blastn passes paths with spaces intact and includes -num_threads."""
    from accio.config import AnalysisConfig
    cfg = AnalysisConfig()
    cfg.THREADS = 12
    runner = BLASTRunner(config=cfg)

    captured_cmd = []
    def mock_check_output(cmd, **kwargs):
        captured_cmd.extend(cmd)
        return b""

    monkeypatch.setattr("subprocess.check_output", mock_check_output)
    runner.run_blastn("/path/with spaces/query.fasta", "/path/with spaces/db")

    assert "/path/with spaces/query.fasta" in captured_cmd
    assert "/path/with spaces/db" in captured_cmd
    assert "-num_threads" in captured_cmd
    assert captured_cmd[captured_cmd.index("-num_threads") + 1] == "12"


def test_nucmer_handles_paths_with_spaces(monkeypatch, tmp_path):
    """Test that NucmerRunner.run_nucmer passes paths with spaces intact without splitting them."""
    from accio.wrappers.external_tools import NucmerRunner
    runner = NucmerRunner()

    prefix = tmp_path / "space dir" / "out_pref"
    prefix.parent.mkdir(parents=True)
    (Path(f"{prefix}.coords")).write_text("")
    (Path(f"{prefix}.snps")).write_text("")

    captured_cmds = []
    def mock_check_call(cmd, **kwargs):
        captured_cmds.append(cmd)

    monkeypatch.setattr("subprocess.check_call", mock_check_call)
    runner.run_nucmer("/ref with spaces/ref.fa", "/query with spaces/q.fa", str(prefix))

    assert captured_cmds[0][0] == "nucmer"
    assert "/ref with spaces/ref.fa" in captured_cmds[0]
    assert "/query with spaces/q.fa" in captured_cmds[0]


def test_read_mapping_runner_uses_mapping_preset(monkeypatch, tmp_path):
    """Test that ReadMappingRunner.map_reads uses config.MAPPING_PRESET when use_minimap=True."""
    from accio.wrappers.external_tools import ReadMappingRunner
    from accio.config import AnalysisConfig

    cfg = AnalysisConfig()
    cfg.MAPPING_PRESET = "sr"
    cfg.THREADS = 6
    runner = ReadMappingRunner(config=cfg)

    popen_cmds = []
    class MockProc:
        def __init__(self, cmd, **kwargs):
            popen_cmds.append(cmd)
            self.stdout = MagicMock()
            self.returncode = 0
        def wait(self):
            return 0

    monkeypatch.setattr("subprocess.Popen", MockProc)
    monkeypatch.setattr("subprocess.check_call", lambda cmd, **kw: None)

    out_bam = runner.map_reads("ref.fa", ["r1.fq", "r2.fq"], str(tmp_path / "test"), use_minimap=True)
    assert out_bam.endswith(".bam")
    minimap_cmd = popen_cmds[0]
    assert minimap_cmd[:5] == ["minimap2", "-a", "-P", "-x", "sr"]
    assert "-t" in minimap_cmd and minimap_cmd[minimap_cmd.index("-t") + 1] == "6"


def test_mash_screen_and_skani_dist_pass_threads(monkeypatch):
    """Test that MashRunner.run_screen and SkaniRunner.run_ska_dist pass configured threads."""
    from accio.wrappers.external_tools import MashRunner, SkaniRunner
    from accio.config import AnalysisConfig

    cfg = AnalysisConfig()
    cfg.THREADS = 10
    mash_runner = MashRunner(config=cfg)
    skani_runner = SkaniRunner(config=cfg)

    recorded_cmds = []
    def mock_run(cmd, **kwargs):
        recorded_cmds.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout=b"", stderr=b"")

    monkeypatch.setattr("subprocess.run", mock_run)
    mash_runner.run_screen("db.msh", "query.fa")
    skani_runner.run_ska_dist("plasmids.fa", "out.tsv")

    assert recorded_cmds[0] == ["mash", "screen", "-p", "10", "db.msh", "query.fa"]
    assert "-t" in recorded_cmds[1] and recorded_cmds[1][recorded_cmds[1].index("-t") + 1] == "10"



