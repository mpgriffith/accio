"""Tests for filtering utilities in accio.utils.filters."""

import pytest
import pandas as pd
from accio.utils.filters import filter_overlaps


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

