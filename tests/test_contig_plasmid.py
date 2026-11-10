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

