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
        if hasattr(pd.errors, "SettingWithCopyWarning"):
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


def test_clean_up_files_handles_nested_directories_and_preserves_novel_fasta(base_config, tmp_path):
    """Test that _clean_up_files removes nested directories without error and keeps novel_plasmid_contigs.fasta."""
    wf = PlasmidAnalysisWorkflow(base_config)

    # Create nested directory like genomad output
    nested_dir = tmp_path / "genomad" / "sample1_summary"
    nested_dir.mkdir(parents=True)
    (nested_dir / "sample1_plasmid_summary.tsv").write_text("seq_name\tplasmid_score\n")

    # Create novel plasmid fasta
    novel_fasta = tmp_path / "sample1_novel_plasmid_contigs.fasta"
    novel_fasta.write_text(">c_novel\nATGC\n")

    contig = Contig("c_novel", SeqRecord(Seq("ATGC" * 100), id="c_novel"), circular=True)
    assignment_results = {'plasmids_assigned': [], 'contigs_assigned': {}}
    novel_results = {
        'contigs': [contig],
        'mob_results': pd.DataFrame(),
        'contig_results': pd.DataFrame([{'contig_id': 'c_novel', 'molecule_type': 'chromosome'}])
    }

    output_files = wf._generate_outputs(
        assignment_results, novel_results, {"c_novel": contig},
        str(tmp_path), "sample1"
    )
    assert 'novel_plasmids_fasta' in output_files

    wf._clean_up_files(str(tmp_path), list(output_files.values()))
    assert not (tmp_path / "genomad").exists()
    assert novel_fasta.exists()


def test_load_contigs_circular_only(base_config, tmp_path):
    """Test that CIRCULAR_ONLY filters out linear contigs during _load_contigs."""
    base_config.CIRCULAR_ONLY = True
    wf = PlasmidAnalysisWorkflow(base_config)

    fasta_file = tmp_path / "asm.fasta"
    fasta_file.write_text(">c_circ circular=true\nATGC\n>c_lin\nATGC\n")

    contigs = wf._load_contigs(str(fasta_file))
    assert "c_circ" in contigs
    assert "c_lin" not in contigs


def test_generate_outputs_partial_mob_and_per_plasmid_fasta(base_config, tmp_path):
    """Test _generate_outputs when mob_results is missing a cluster ID and verify per-plasmid FASTA creation."""
    from accio.core.plasmid import Plasmid
    wf = PlasmidAnalysisWorkflow(base_config)

    c_assigned = Contig("c_assigned", SeqRecord(Seq("ATGC" * 100), id="c_assigned"))
    c_novel = Contig("c_novel", SeqRecord(Seq("GGCC" * 100), id="c_novel"), circular=True)
    c_novel.add_rep_type(pd.DataFrame([{'sseqid': 'IncFIB_1'}]))
    c_novel.add_amr_genes(['blaTEM-1'])

    p_assigned = Plasmid("pRef1", "cl1", "subcomm_1", 400, "db.msh")
    p_assigned.contigs = {"c_assigned": c_assigned}

    assignment_results = {'plasmids_assigned': [p_assigned], 'contigs_assigned': {'c_assigned': 'pRef1'}}
    # mob_contigs has cluster 'novel_cl_1', while mob_results has a different cluster (testing fallback without IndexError)
    mob_results = pd.DataFrame([{'primary_cluster_id': 'other_cluster', 'size': 400}])
    mob_contigs = pd.DataFrame([{'contig_id': 'c_novel', 'molecule_type': 'plasmid', 'primary_cluster_id': 'novel_cl_1', 'size': 400}])
    novel_results = {'contigs': [c_novel], 'mob_results': mob_results, 'contig_results': mob_contigs}

    output_files = wf._generate_outputs(
        assignment_results, novel_results,
        {"c_assigned": c_assigned, "c_novel": c_novel},
        str(tmp_path), "sampleX"
    )

    assert 'plasmids' in output_files
    assert (tmp_path / "sampleX_plasmid_subcomm_1_contigs.fasta").exists()
    assert (tmp_path / "sampleX_plasmid_novel_cl_1_contigs.fasta").exists()


def test_process_plasmids_replicon_filtering(base_config):
    """Test that _process_plasmids filters out plasmids whose rep_type is not found in the sample."""
    from accio.core.plasmid import Plasmid
    wf = PlasmidAnalysisWorkflow(base_config)

    c1 = Contig("c1", SeqRecord(Seq("ATGC" * 100), id="c1"))
    c1.add_rep_type(pd.DataFrame([{'sseqid': 'IncFII_1'}]))
    c2 = Contig("c2", SeqRecord(Seq("GGCC" * 100), id="c2"))  # rep_types is None by default

    p_match = Plasmid("p_match", "c1", "pt1", 1000, "db.msh", rep_type="IncFII")
    p_untyped = Plasmid("p_untyped", "c2", "pt2", 1000, "db.msh", rep_type=None)
    p_mismatch = Plasmid("p_mismatch", "c3", "pt3", 1000, "db.msh", rep_type="IncX3")

    plasmids = {"p_match": p_match, "p_untyped": p_untyped, "p_mismatch": p_mismatch}
    blast_results = {
        'blast': pd.DataFrame([{'qseqid': 'c1', 'sseqid': 'p_match'}]),
        'mash': pd.DataFrame()
    }

    matches = wf._process_plasmids({"c1": c1, "c2": c2}, plasmids, blast_results, coverage_results=None)
    assert "p_match" in plasmids
    assert "p_untyped" in plasmids
    assert "p_mismatch" not in plasmids
    assert matches == ["p_match"]


def test_get_median_coverage_and_load_plasmids(base_config, tmp_path):
    """Test _get_median_coverage length-weighted calculation and _load_plasmids from CSV."""
    wf = PlasmidAnalysisWorkflow(base_config)

    assert wf._get_median_coverage(pd.DataFrame()) == 30.0

    cov_df = pd.DataFrame([
        {'contig': 'c1', 'length': 100000, 'meandepth': 45.0},
        {'contig': 'c2', 'length': 20000, 'meandepth': 120.0},
        {'contig': 'c3', 'length': 10000, 'meandepth': 90.0},
    ])
    # c1 has >50% of total length (100k / 130k), so weighted median is 45.0
    assert wf._get_median_coverage(cov_df) == 45.0

    # Test _load_plasmids
    info_df = pd.DataFrame([{
        'sample_id': 'plasmid_1',
        'pling_type': 'subcomm_A',
        'primary_cluster_id': 'AA100',
        'size': 55000,
        'rep_type': 'IncFIA',
        'amr_genes': 'blaKPC-2',
        'subcommunity_amr_genes': 'blaKPC-2,tet(A)'
    }])
    info_df.to_csv(tmp_path / "plasmidDB_info.csv", index=False)

    plasmids = wf._load_plasmids(str(tmp_path))
    assert "plasmid_1" in plasmids
    assert plasmids["plasmid_1"].rep_type == "IncFIA"
    assert plasmids["plasmid_1"].length == 55000




