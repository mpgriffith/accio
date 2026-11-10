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

