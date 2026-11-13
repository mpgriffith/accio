"""Tests for filtering utilities in accio.utils.filters."""

import pytest
import pandas as pd
from accio.utils.filters import (
    filter_overlaps,
    filter_blast_results,
    filter_by_similarity_score,
    filter_plasmid_candidates,
    remove_low_quality_alignments,
)


def test_filter_overlaps_empty():
    empty_df = pd.DataFrame()
    result = filter_overlaps(empty_df)
    assert result.empty


def test_filter_overlaps_non_overlapping():
    df = pd.DataFrame([
        {
            'qseqid': 'contig_1', 'sseqid': 'ref1',
            'qstart': 1, 'qend': 1000,
            'bitscore': 500, 'pident': 99.0
        },
        {
            'qseqid': 'contig_1', 'sseqid': 'ref2',
            'qstart': 1500, 'qend': 2500,
            'bitscore': 450, 'pident': 98.0
        }
    ])
    result = filter_overlaps(df, by='qseqid')
    assert len(result) == 2


def test_filter_overlaps_removes_worse_overlapping_hit():
    df = pd.DataFrame([
        {
            'qseqid': 'contig_1', 'sseqid': 'ref1',
            'qstart': 100, 'qend': 2000,
            'bitscore': 1000, 'pident': 99.5
        },
        {
            'qseqid': 'contig_1', 'sseqid': 'ref2',
            'qstart': 150, 'qend': 1950,
            'bitscore': 800, 'pident': 95.0
        }
    ])
    result = filter_overlaps(df, by='qseqid', overlap_threshold=0.8)
    assert len(result) == 1
    assert result.iloc[0]['sseqid'] == 'ref1'
    assert result.iloc[0]['bitscore'] == 1000


def test_filter_overlaps_grouped_by_target():
    df = pd.DataFrame([
        {
            'qseqid': 'c1', 'sseqid': 'ref_A',
            'sstart': 100, 'send': 1000,
            'bitscore': 600, 'pident': 99.0
        },
        {
            'qseqid': 'c2', 'sseqid': 'ref_A',
            'sstart': 100, 'send': 1000,
            'bitscore': 900, 'pident': 100.0
        }
    ])
    result = filter_overlaps(df, by='sseqid', pos_cols=['sstart', 'send'])
    assert len(result) == 1
    assert result.iloc[0]['qseqid'] == 'c2'


def test_filter_blast_results():
    assert filter_blast_results(pd.DataFrame()).empty

    df = pd.DataFrame([
        {'qseqid': 'c1', 'pident': 95.0, 'qcovs': 90.0, 'length': 1200},
        {'qseqid': 'c2', 'pident': 70.0, 'qcovs': 95.0, 'length': 1500},  # low identity
        {'qseqid': 'c3', 'pident': 90.0, 'qcovs': 50.0, 'length': 1500},  # low coverage
        {'qseqid': 'c4', 'pident': 92.0, 'qcovs': 85.0, 'length': 300},   # short length
    ])
    res = filter_blast_results(df, min_identity=80.0, min_coverage=80.0, min_length=500)
    assert list(res['qseqid']) == ['c1']


def test_filter_by_similarity_score():
    assert filter_by_similarity_score(pd.DataFrame()).empty

    df = pd.DataFrame([
        {'qseqid': 'c1', 'pident': 100.0, 'length': 900, 'slen': 1000},  # sim_score = 90
        {'qseqid': 'c2', 'pident': 80.0, 'length': 500, 'slen': 1000},   # sim_score = 40
    ])
    res = filter_by_similarity_score(df, min_score=80.0)
    assert list(res['qseqid']) == ['c1']
    assert abs(res.iloc[0]['sim_score'] - 90.0) < 1e-5

    # Invalid formula falls back to returning original DataFrame without crashing
    fallback = filter_by_similarity_score(df, score_formula='nonexistent_col * 2')
    assert len(fallback) == 2


def test_filter_plasmid_candidates():
    assert filter_plasmid_candidates(pd.DataFrame()).empty

    df = pd.DataFrame([
        {'contig_id': 'c1', 'length': 50000, 'circular': True},
        {'contig_id': 'c2', 'length': 100000, 'circular': False},
        {'contig_id': 'c3', 'length': 500000, 'circular': True},
    ])
    res_all = filter_plasmid_candidates(df, max_length=300000, require_circular=False)
    assert list(res_all['contig_id']) == ['c1', 'c2']

    res_circ = filter_plasmid_candidates(df, max_length=300000, require_circular=True)
    assert list(res_circ['contig_id']) == ['c1']


def test_remove_low_quality_alignments():
    assert remove_low_quality_alignments(pd.DataFrame()).empty

    df = pd.DataFrame([
        {'qseqid': 'c1', 'bitscore': 500.0, 'evalue': 1e-50},
        {'qseqid': 'c2', 'bitscore': 50.0, 'evalue': 1e-50},
        {'qseqid': 'c3', 'bitscore': 500.0, 'evalue': 0.5},
    ])
    res = remove_low_quality_alignments(df, min_bitscore=100.0, max_evalue=1e-5)
    assert list(res['qseqid']) == ['c1']


