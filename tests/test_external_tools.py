"""Tests for external tool wrappers in accio.wrappers.external_tools."""

import pytest
import subprocess
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

