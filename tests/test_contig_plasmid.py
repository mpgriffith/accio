"""Tests for Contig and Plasmid classes in accio.core."""

import pytest
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from accio.core.contig import Contig
from accio.core.plasmid import Plasmid
from accio.config import AnalysisConfig


def test_contig_initialization():
    seq = SeqRecord(Seq("ATGC" * 100), id="contig_1", description="circular=true")
    c = Contig("contig_1", seq, circular=True)

    assert c.id == "contig_1"
    assert c.length == 400
    assert c.circular is True
    assert c.type is True  # circular contigs default to True
    assert c.assigned is False


def test_contig_add_amr_genes():
    seq = SeqRecord(Seq("ATGC" * 50), id="contig_2")
    c = Contig("contig_2", seq)
    c.add_amr_genes(["blaCTX-M-15", "aac(6')-Ib"])

    assert c.amr_genes == ["blaCTX-M-15", "aac(6')-Ib"]
    summary = c.get_summary()
    assert "blaCTX-M-15" in summary['amr_genes']


def test_contig_add_rep_type():
    seq = SeqRecord(Seq("ATGC" * 50), id="contig_3")
    c = Contig("contig_3", seq)
    
    rep_df = pd.DataFrame([
        {'sseqid': 'IncFII_1_pKPX-1', 'pident': 99.0, 'length': 200, 'slen': 200}
    ])
    c.add_rep_type(rep_df)
    assert c.rep_types == ['IncFII']


def test_plasmid_initialization_and_summary():
    p = Plasmid(
        id="plasmid_ref_1",
        cluster_id="cluster_10",
        pling_type="type_X",
        length=50000,
        plasmid_mash_db="db.msh",
        rep_type="IncFII,IncFIA",
        amr_genes="blaKPC-2,tet(A)"
    )

    assert p.id == "plasmid_ref_1"
    assert p.cluster_id == "cluster_10"
    assert p.pling_type == "type_X"
    assert p.length == 50000
    assert p.rep_type == "IncFII,IncFIA"

    summary = p.get_summary()
    assert summary['plasmid_id'] == "plasmid_ref_1"
    assert summary['cluster_id'] == "cluster_10"
    assert summary['pling_type'] == "type_X"
    assert summary['length'] == 50000


def test_plasmid_add_contigs_and_score():
    p = Plasmid(
        id="plasmid_ref_1",
        cluster_id="c1",
        pling_type="t1",
        length=5000,
        plasmid_mash_db="db.msh"
    )

    seq = SeqRecord(Seq("ATGC" * 1250), id="c1")
    contig = Contig("c1", seq)
    
    # Add blast hit covering entire plasmid
    contig.blast_data = pd.DataFrame([
        {
            'qseqid': 'c1',
            'sseqid': 'plasmid_ref_1',
            'pident': 100.0,
            'length': 5000,
            'mismatch': 0,
            'gapopen': 0,
            'qstart': 1,
            'qend': 5000,
            'sstart': 1,
            'send': 5000,
            'evalue': 0.0,
            'bitscore': 9500.0,
            'qlen': 5000,
            'slen': 5000,
            'stitle': 'Ref Plasmid',
            'qcovs': 100.0,
            'qframe': 1,
            'sframe': 1,
            'nident': 5000,
            'use_hit': True
        }
    ])

    p.add_contigs(contig)
    assert "c1" in p.contigs
    assert p.scores['blast_coverage'] > 90
    assert p.scores['length_score'] > 0.9
    assert p.score > 0


def test_mash_temp_fasta_uses_tempfile(monkeypatch):
    """Test that Plasmid.calc_mash_score uses NamedTemporaryFile in tempdir and cleans up."""
    import os
    import tempfile
    p = Plasmid(
        id="plasmid_test",
        cluster_id="c1",
        pling_type="pt1",
        length=5000,
        plasmid_mash_db="dummy.msh"
    )
    contig = Contig("c1", SeqRecord(Seq("ATGC" * 100), id="c1"))
    p.contigs["c1"] = contig

    created_files = []
    real_named_temp = tempfile.NamedTemporaryFile

    def tracking_tempfile(*args, **kwargs):
        tf = real_named_temp(*args, **kwargs)
        created_files.append(tf.name)
        return tf

    monkeypatch.setattr(tempfile, "NamedTemporaryFile", tracking_tempfile)
    monkeypatch.setattr(
        "subprocess.check_output",
        lambda *args, **kwargs: b"0.98\t450/1000\t1\t0.0\tplasmid_test\tplasmid_test\n"
    )

    score = p.calc_mash_score()
    assert abs(score - 0.98) < 1e-4
    assert len(created_files) == 1
    temp_path = created_files[0]
    assert "accio_mash_plasmid_test_" in os.path.basename(temp_path)
    assert not os.path.exists(temp_path), "Temporary fasta should be deleted in finally block"


def test_contig_minus_strand_repetitive_data():
    """Test that minus-strand BLAST hits (sstart > send) are properly counted in add_repetitive_data."""
    seq = SeqRecord(Seq("ATGC" * 250), id="c_rep")  # 1000 bp
    c = Contig("c_rep", seq)
    rep_df = pd.DataFrame([
        {
            'qseqid': 'c_rep',
            'sseqid': 'IS26',
            'pident': 99.0,
            'length': 820,
            'qstart': 100,
            'qend': 919,
            'sstart': 820,  # minus strand: sstart > send
            'send': 1,
            'slen': 820
        }
    ])
    c.add_repetitive_data(rep_df)
    assert c.contig_type == 'repetitive'
    assert not c.rep_data.empty
    assert c.rep_data.iloc[0]['rep_elem_cov'] > 95.0


def test_repetitive_contig_assign_and_repr():
    """Test that assigning a repetitive contig does not retain placeholder (None, None) or crash __repr__."""
    seq = SeqRecord(Seq("ATGC" * 200), id="c_rep")
    c = Contig("c_rep", seq)
    c.contig_type = 'repetitive'

    assert "unassigned" in repr(c)
    c.assign("p1", "cl1")
    c.assign("p2", "cl2")
    assert c.plas_data == [("p1", "cl1"), ("p2", "cl2")]
    assert "p1, p2" in repr(c)
    summary = c.get_summary()
    assert summary['assigned plasmid(s)'] == ["p1", "p2"]


def test_contig_add_blast_data_and_copy_num():
    """Test Contig.add_blast_data filtering (linear vs circular) and add_copy_num threshold."""
    seq = SeqRecord(Seq("ATGC" * 250), id="c_lin")  # 1000 bp
    c_lin = Contig("c_lin", seq, circular=False)

    df = pd.DataFrame([
        # High coverage hit -> kept
        {'qseqid': 'c_lin', 'sseqid': 'p1', 'pident': 95.0, 'qcovs': 85.0, 'qstart': 100, 'qend': 950, 'qlen': 1000, 'bitscore': 800, 'sstart': 1, 'send': 850},
        # Low coverage hit near contig start (qstart < 50) -> kept
        {'qseqid': 'c_lin', 'sseqid': 'p2', 'pident': 95.0, 'qcovs': 20.0, 'qstart': 10, 'qend': 210, 'qlen': 1000, 'bitscore': 300, 'sstart': 1, 'send': 200},
        # Low coverage hit in middle -> filtered out
        {'qseqid': 'c_lin', 'sseqid': 'p3', 'pident': 95.0, 'qcovs': 20.0, 'qstart': 300, 'qend': 500, 'qlen': 1000, 'bitscore': 300, 'sstart': 1, 'send': 200},
        # Below MIN_IDENTITY -> filtered out
        {'qseqid': 'c_lin', 'sseqid': 'p4', 'pident': 50.0, 'qcovs': 95.0, 'qstart': 1, 'qend': 950, 'qlen': 1000, 'bitscore': 300, 'sstart': 1, 'send': 950},
    ])
    c_lin.add_blast_data(df)
    assert set(c_lin.blast_data['sseqid']) == {'p1', 'p2'}

    # Copy number above HIGH_COPY_THRESHOLD (1.15) marks contig repetitive
    c_lin.add_copy_num(2.5)
    assert c_lin.contig_type == 'repetitive'


def test_contig_add_plasme_and_genomad_data():
    """Test Contig.add_plasme_data (ambiguous region & Transformer checks) and add_genomad_data."""
    seq = SeqRecord(Seq("ATGC" * 250), id="c1")
    c_amb = Contig("c1", seq)
    # >50% ambiguous region -> rejected
    c_amb.add_plasme_data(pd.DataFrame([{'score': 0.95, 'length': 1000, 'amb_region': '1-600', 'evidence': 'BLAST'}]), threshold=0.75)
    assert c_amb.type is None

    # Transformer with score <= 0.75 -> rejected even if threshold is 0.5
    c_trans = Contig("c1", seq)
    c_trans.add_plasme_data(pd.DataFrame([{'score': 0.70, 'length': 1000, 'amb_region': '', 'evidence': 'Transformer'}]), threshold=0.5)
    assert c_trans.type is None

    # Valid PLASMe score -> accepted
    c_ok = Contig("c1", seq)
    c_ok.add_plasme_data(pd.DataFrame([{'score': 0.90, 'length': 1000, 'amb_region': '1-100', 'evidence': 'Transformer'}]), threshold=0.75)
    assert c_ok.type == 'plasmid'
    assert c_ok.classifier_scores['plasme'] == 0.90

    # Valid geNomad score -> accepted
    c_gen = Contig("c1", seq)
    c_gen.add_genomad_data(pd.DataFrame([{'plasmid_score': 0.85}]), threshold=0.75)
    assert c_gen.type == 'plasmid'
    assert c_gen.classifier_scores['genomad'] == 0.85


def test_plasmid_remove_contigs_preserves_config_and_resets_scores(monkeypatch):
    """Test that remove_contigs keeps self.config settings and resets scores when all contigs are removed."""
    cfg = AnalysisConfig()
    cfg.USE_MINIMAP = True
    cfg.MAPPING_PRESET = 'sr'
    cfg.READS_SUPPORT_ONLY = True

    p = Plasmid("p1", "cl1", "pt1", 1000, "db.msh")
    monkeypatch.setattr(p, "calc_mash_score", lambda: 0.0)

    c1 = Contig("c1", SeqRecord(Seq("ATGC" * 250), id="c1"))
    c1.blast_data = pd.DataFrame([{
        'qseqid': 'c1', 'sseqid': 'p1', 'pident': 100.0, 'length': 1000,
        'mismatch': 0, 'gapopen': 0, 'qstart': 1, 'qend': 1000,
        'sstart': 1, 'send': 1000, 'evalue': 0.0, 'bitscore': 2000.0,
        'qlen': 1000, 'slen': 1000, 'stitle': 'p1', 'qcovs': 100.0,
        'qframe': 1, 'sframe': 1, 'nident': 1000, 'use_hit': True
    }])

    p.add_contigs(c1, config=cfg)
    assert p.use_minimap is True
    assert p.reads_support_only is True
    assert p.scores['length_score'] == 1.0

    p.remove_contigs("c1")
    assert len(p.contigs) == 0
    assert p.use_minimap is True
    assert p.reads_support_only is True
    assert p.scores['blast_coverage'] == 0
    assert p.scores['length_score'] == 0


def test_plasmid_synteny_and_read_support_calculations():
    """Test block/nucmer synteny scoring, find_gaps, and calculate_read_support."""
    p = Plasmid("p1", "cl1", "pt1", 1000, "db.msh")
    c1 = Contig("c1", SeqRecord(Seq("ATGC" * 250), id="c1"), length=1000)

    nucmer_df = pd.DataFrame([
        {
            'sstart': 1, 'send': 500, 'qstart': 1, 'qend': 500,
            'ref_aln_len': 500, 'query_aln_len': 500, 'pident': 99.0,
            'slen': 1000, 'qlen': 1000, 'scov': 50.0, 'qcovs': 50.0,
            'sframe': 1, 'qframe': 1, 'sseqid': 'p1', 'qseqid': 'c1'
        },
        {
            'sstart': 501, 'send': 1000, 'qstart': 501, 'qend': 1000,
            'ref_aln_len': 500, 'query_aln_len': 500, 'pident': 99.0,
            'slen': 1000, 'qlen': 1000, 'scov': 50.0, 'qcovs': 50.0,
            'sframe': 1, 'qframe': 1, 'sseqid': 'p1', 'qseqid': 'c1'
        }
    ])
    c1.add_nucmer_data(nucmer_df, pd.DataFrame())
    p.contigs['c1'] = c1

    p.calc_synteny_score(synteny_type='blocks')
    assert p.scores['synteny_score'] > 0
    assert p.scores['length_score'] == 1.0

    p.calc_synteny_score(synteny_type='nucmer')
    assert p.scores['synteny_score'] > 0

    # Mock read_cov_data for calculate_read_support
    class MockBam:
        def count_coverage(self, ref_id):
            # Depth of 12 across all 1000 positions (3+3+3+3)
            return ([3] * 1000, [3] * 1000, [3] * 1000, [3] * 1000)

    p.read_cov_data = MockBam()
    support_short = p.calculate_read_support(use_minimap=False)
    assert support_short['read_support'] is True
    assert support_short['read_support_pct'] == 100.0
    assert support_short['read_support_mean_depth'] == 12.0

    support_long = p.calculate_read_support(use_minimap=True, mapping_preset='map-ont')
    assert support_long['read_support'] is True



