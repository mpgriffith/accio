"""Pytest configuration and shared fixtures for Accio tests."""

import pytest
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from accio.config import AnalysisConfig


@pytest.fixture
def base_config():
    """Returns a default AnalysisConfig instance."""
    return AnalysisConfig()


@pytest.fixture
def sample_seqrecord():
    """Returns a simple SeqRecord for testing."""
    seq = "ATGC" * 2500  # 10,000 bp
    return SeqRecord(Seq(seq), id="plasmid_test_1", description="test sequence circular=true")


@pytest.fixture
def circular_overlap_seqrecord():
    """Returns a SeqRecord with a 50bp terminal overlap to simulate circular assembly."""
    overlap = "ACGT" * 12 + "AC"  # 50 bp
    middle = "TTAA" * 1000
    seq = overlap + middle + overlap
    return SeqRecord(Seq(seq), id="circ_overlap_1", description="overlap test")


@pytest.fixture
def sample_blast_df():
    """Returns a mock BLAST results DataFrame."""
    return pd.DataFrame([
        {
            'qseqid': 'contig_1',
            'sseqid': 'ref_plasmid_A',
            'pident': 99.5,
            'length': 5000,
            'mismatch': 5,
            'gapopen': 0,
            'qstart': 1,
            'qend': 5000,
            'sstart': 100,
            'send': 5100,
            'evalue': 0.0,
            'bitscore': 9200.0,
            'qlen': 5000,
            'slen': 50000,
            'stitle': 'Reference Plasmid A',
            'qcovs': 100.0,
            'qframe': 1,
            'sframe': 1,
            'nident': 4975
        }
    ])

