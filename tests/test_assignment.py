"""Tests for plasmid assignment algorithm in accio.core.assignment."""

import pytest
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from accio.core.assignment import assign_plasmids
from accio.core.contig import Contig
from accio.core.plasmid import Plasmid
from accio.config import AnalysisConfig


def test_assign_plasmids_empty():
    assigned_plasmids, assigned_contigs = assign_plasmids({}, {})
    assert assigned_plasmids == []
    assert assigned_contigs == {}


def test_assign_plasmids_single_match():
    config = AnalysisConfig()
    config.MIN_PLASMID_SCORE = 0.5
    config.MIN_BLAST_SIMILARITY = 50

    plasmid = Plasmid(
        id="pRef1",
        cluster_id="cl1",
        pling_type="pt1",
        length=2000,
        plasmid_mash_db="db.msh"
    )

    contig = Contig("c1", SeqRecord(Seq("ATGC" * 500), id="c1"))
    contig.blast_data = pd.DataFrame([
        {
            'qseqid': 'c1',
            'sseqid': 'pRef1',
            'pident': 100.0,
            'length': 2000,
            'mismatch': 0,
            'gapopen': 0,
            'qstart': 1,
            'qend': 2000,
            'sstart': 1,
            'send': 2000,
            'evalue': 0.0,
            'bitscore': 4000.0,
            'qlen': 2000,
            'slen': 2000,
            'stitle': 'pRef1',
            'qcovs': 100.0,
            'qframe': 1,
            'sframe': 1,
            'nident': 2000,
            'use_hit': True
        }
    ])

    plasmid.add_contigs(contig)
    
    assigned_plasmids, assigned_contigs = assign_plasmids(
        {"pRef1": plasmid}, {"c1": contig}, config=config
    )

    assert len(assigned_plasmids) == 1
    assert assigned_plasmids[0].id == "pRef1"
    assert assigned_contigs["c1"] == "pRef1"
    assert contig.assigned is True

