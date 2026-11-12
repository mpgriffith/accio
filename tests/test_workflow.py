"""Tests for PlasmidAnalysisWorkflow in accio.core.workflow."""

import pytest
import pandas as pd
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from unittest.mock import MagicMock

from accio.config import AnalysisConfig
from accio.core.workflow import PlasmidAnalysisWorkflow
from accio.core.contig import Contig


@pytest.fixture
def workflow(base_config):
    return PlasmidAnalysisWorkflow(base_config)


def test_load_contigs(workflow, tmp_path):
    fasta_file = tmp_path / "assembly.fasta"
    seq1 = SeqRecord(Seq("ATGC" * 100), id="contig_circ", description="circ contig circular=true")
    seq2 = SeqRecord(Seq("CGTA" * 100), id="contig_linear", description="linear contig")
    # Sequence larger than 1.5 * MAX_PLASMID_LENGTH (450,000 bp)
    seq3 = SeqRecord(Seq("AAAA" * 150000), id="chromosome_contig", description="very large contig")

    from Bio import SeqIO
    with open(fasta_file, "w") as f:
        SeqIO.write([seq1, seq2, seq3], f, "fasta")

    contigs = workflow._load_contigs(str(fasta_file))

    assert "contig_circ" in contigs
    assert contigs["contig_circ"].circular is True
    assert "contig_linear" in contigs
    assert contigs["contig_linear"].circular is False
    assert "chromosome_contig" not in contigs  # filtered out by length


def test_process_novel_plasmids_skipped_when_flag_false(base_config):
    base_config.FIND_NOVEL_PLASMIDS = False
    wf = PlasmidAnalysisWorkflow(base_config)
    wf._run_and_parse_mob_recon = MagicMock()

    contig = Contig("c1", SeqRecord(Seq("ATGC"), id="c1"))
    contig.type = "plasmid"
    contig.assigned = False

    result = wf._process_novel_plasmids({"c1": contig}, "out_dir", "sample1")
    assert result['contigs'] == []
    assert result['mob_results'].empty
    assert result['contig_results'].empty
    wf._run_and_parse_mob_recon.assert_not_called()


def test_generate_outputs_skip_novel_plasmids(base_config, tmp_path):
    base_config.FIND_NOVEL_PLASMIDS = False
    wf = PlasmidAnalysisWorkflow(base_config)

    contig = Contig("c1", SeqRecord(Seq("ATGC" * 100), id="c1"))
    contig.type = "plasmid"
    contig.assigned = False

    assignment_results = {'plasmids_assigned': [], 'contigs_assigned': {}}
    novel_results = {'contigs': [], 'mob_results': pd.DataFrame(), 'contig_results': pd.DataFrame()}

    output_files = wf._generate_outputs(
        assignment_results, novel_results, {"c1": contig},
        str(tmp_path), "sample1"
    )

    # Plasmids file should be generated
    assert 'plasmids' in output_files
    plasmids_df = pd.read_csv(output_files['plasmids'])
    # Should not contain any novel_plasmid rows
    assert not (plasmids_df['ACCIO_Plasmids'] == 'novel_plasmid').any()

    # Contigs summary should have novel_plasmid_contig = False
    assert 'contigs' in output_files
    contigs_df = pd.read_csv(output_files['contigs'])
    assert (contigs_df['novel_plasmid_contig'] == False).all()


def test_sequence_analysis_skips_classifiers_when_flag_false(base_config, tmp_path):
    base_config.FIND_NOVEL_PLASMIDS = False
    wf = PlasmidAnalysisWorkflow(base_config)

    # Mock tool runners
    wf.blast_runner.run_blastn = MagicMock(return_value=pd.DataFrame())
    wf.mash_runner.run_screen = MagicMock(return_value=pd.DataFrame())
    wf.nucmer_runner.run_nucmer = MagicMock(return_value=(pd.DataFrame(), pd.DataFrame(), {}))
    wf._run_amrfinder = MagicMock(return_value=pd.DataFrame())
    wf.plasme_runner.run_plasme = MagicMock()
    wf.genomad_runner.run_genomad = MagicMock()

    fasta = tmp_path / "test.fasta"
    fasta.write_text(">c1\nATGC\n")

    results = wf._run_sequence_analysis(
        str(fasta), str(tmp_path), str(tmp_path), "sample1", keep_intermediate=False
    )

    # PLASMe and geNomad runners should NOT be called
    wf.plasme_runner.run_plasme.assert_not_called()
    wf.genomad_runner.run_genomad.assert_not_called()
    assert results['plasme'].empty
    assert results['genomad'].empty


def test_process_contigs_amr_data(base_config):
    import warnings
    wf = PlasmidAnalysisWorkflow(base_config)

    contig = Contig("c1", SeqRecord(Seq("ATGC" * 100), id="c1"))
    contigs = {"c1": contig}

    amr_genes_df = pd.DataFrame([
        {
            'Contig id': 'c1',
            'Element symbol': 'blaCTX-M-15',
            '% Identity to reference': 100.0,
            '% Coverage of reference': 100.0
        },
        {
            'Contig id': 'c1',
            'Element symbol': 'low_match_gene',
            '% Identity to reference': 60.0,
            '% Coverage of reference': 60.0  # sim_score = 36 <= 80
        }
    ])

    blast_results = {
        'plasme': pd.DataFrame(),
        'genomad': pd.DataFrame(),
        'plasmidfinder': pd.DataFrame(),
        'repetitive': pd.DataFrame(),
        'amr_genes': amr_genes_df,
        'nucmer_coords': pd.DataFrame(columns=['qseqid']),
        'nucmer_snps': pd.DataFrame(columns=['query']),
    }

    with warnings.catch_warnings():
        warnings.simplefilter("error", pd.errors.SettingWithCopyWarning)
        wf._process_contigs(contigs, blast_results, coverage_data=None)

    assert contig.amr_genes == ['blaCTX-M-15']


def test_classifier_failure_records_warnings(base_config, tmp_path):
    """Test that classifier failures are tracked in workflow warnings and status."""
    from accio.wrappers.external_tools import ExternalToolError
    base_config.CLASSIFIER = 'both'
    base_config.FIND_NOVEL_PLASMIDS = True
    wf = PlasmidAnalysisWorkflow(base_config)

    wf.blast_runner.run_blastn = MagicMock(return_value=pd.DataFrame())
    wf.mash_runner.run_screen = MagicMock(return_value=pd.DataFrame())
    wf.nucmer_runner.run_nucmer = MagicMock(return_value=(pd.DataFrame(), pd.DataFrame(), {}))
    wf._run_amrfinder = MagicMock(return_value=pd.DataFrame())

    # Simulate PLASMe and geNomad throwing ExternalToolError
    wf.plasme_runner.run_plasme = MagicMock(side_effect=ExternalToolError("PLASMe out of memory"))
    wf.genomad_runner.run_genomad = MagicMock(side_effect=ExternalToolError("geNomad out of memory"))

    fasta = tmp_path / "test.fasta"
    fasta.write_text(">c1\nATGC\n")

    results = wf._run_sequence_analysis(
        str(fasta), str(tmp_path), str(tmp_path), "sample1", keep_intermediate=False
    )

    assert len(wf.warnings) == 2
    assert any("PLASMe" in w for w in wf.warnings)
    assert any("geNomad" in w for w in wf.warnings)

    summary = wf._create_summary({'plasmids_assigned': []}, {'contigs': []})
    assert summary['status'] == 'completed_with_warnings'
    assert len(summary['warnings']) == 2


def test_genomad_splits_passed_to_cmd(base_config):
    """Test that GENOMAD_SPLITS config passes --splits to genomad command."""
    from unittest.mock import patch
    from accio.wrappers.external_tools import GenomadRunner

    base_config.GENOMAD_SPLITS = 8
    runner = GenomadRunner(base_config)

    with patch("subprocess.check_call") as mock_call:
        runner.run_genomad("in.fa", "out_dir")
        cmd = mock_call.call_args[0][0]
        assert '--splits' in cmd
        splits_idx = cmd.index('--splits')
        assert cmd[splits_idx + 1] == '8'


