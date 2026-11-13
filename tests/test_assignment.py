"""Tests for plasmid assignment algorithm in accio.core.assignment."""

import pytest
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from accio.core.assignment import (
    assign_plasmids,
    _has_cluster_conflict,
    calculate_assignment_confidence,
    resolve_assignment_conflicts,
    validate_assignments,
)
from accio.core.contig import Contig
from accio.core.plasmid import Plasmid
from accio.config import AnalysisConfig


def test_assign_plasmids_empty():
    assigned_plasmids, assigned_contigs = assign_plasmids({}, {})
    assert assigned_plasmids == []
    assert assigned_contigs == {}


def test_assign_plasmids_single_match(monkeypatch):
    config = AnalysisConfig()
    config.MIN_PLASMID_SCORE = 0.5
    config.MIN_BLAST_SIMILARITY = 50
    config.MIN_MASH_SCORE = 0.90

    plasmid = Plasmid(
        id="pRef1",
        cluster_id="cl1",
        pling_type="pt1",
        length=2000,
        plasmid_mash_db="db.msh"
    )
    monkeypatch.setattr(plasmid, "calc_mash_score", lambda: 0.99)

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
    assert contig.assigned == 'plasmid'


def test_select_best_plasmid_respects_min_score():
    from accio.core.assignment import _select_best_plasmid
    p_high = Plasmid(id="p_high", cluster_id="c1", pling_type="typeA", length=1000, plasmid_mash_db="db.msh")
    p_high.score = 0.75
    p_low = Plasmid(id="p_low", cluster_id="c2", pling_type="typeB", length=1000, plasmid_mash_db="db.msh")
    p_low.score = 0.60

    plasmids = {"p_high": p_high, "p_low": p_low}
    # In _select_best_plasmid, pling_ratios = current_series / all_pling_counts
    # So typeB (1/1 = 1.0) has a higher relative ratio than typeA (1/10 = 0.1)
    pling_counts = pd.Series({"typeA": 10, "typeB": 1})

    # With min_score=0.70, only p_high qualifies (score 0.75 >= 0.70)
    best_with_70 = _select_best_plasmid(["p_high", "p_low"], plasmids, pling_counts, min_score=0.70)
    assert best_with_70 == "p_high"

    # With default min_score=0.80, neither qualifies, so fallback considers both and picks typeB (p_low) by pling ratio
    best_with_80 = _select_best_plasmid(["p_high", "p_low"], plasmids, pling_counts, min_score=0.80)
    assert best_with_80 == "p_low"


def test_has_cluster_conflict_and_none_novel_handling():
    """Verify cluster conflict blocks matching clusters but ignores None, 'novel', and 'unknown'."""
    clusters_assigned = {'pling_type': [None, 'novel', 'unknown', 'type_1'], 'mob_cluster': [None, 'novel', 'cluster_A']}

    p_conflict_pling = Plasmid("p1", cluster_id="cluster_B", pling_type="type_1", length=1000, plasmid_mash_db="db.msh")
    assert _has_cluster_conflict(p_conflict_pling, clusters_assigned, []) is True

    p_conflict_mob = Plasmid("p2", cluster_id="cluster_A", pling_type="type_2", length=1000, plasmid_mash_db="db.msh")
    assert _has_cluster_conflict(p_conflict_mob, clusters_assigned, []) is True

    p_novel = Plasmid("p3", cluster_id="novel", pling_type="novel", length=1000, plasmid_mash_db="db.msh")
    assert _has_cluster_conflict(p_novel, clusters_assigned, []) is False

    p_none = Plasmid("p4", cluster_id=None, pling_type=None, length=1000, plasmid_mash_db="db.msh")
    assert _has_cluster_conflict(p_none, clusters_assigned, []) is False


def test_repetitive_contig_shared_across_plasmids(monkeypatch):
    """Verify non-repetitive contigs are exclusive while repetitive contigs can be shared across plasmids."""
    config = AnalysisConfig()
    config.MIN_ASSIGNMENT_SCORE = 0.3
    config.MIN_BLAST_SIMILARITY = 40
    config.MIN_MASH_SCORE = 0.90
    monkeypatch.setattr(Plasmid, "calc_mash_score", lambda self: 0.99)

    c_unique = Contig("c_unique", SeqRecord(Seq("ATGC" * 500), id="c_unique"))
    c_rep = Contig("c_rep", SeqRecord(Seq("GGCC" * 250), id="c_rep"))
    c_rep.contig_type = 'repetitive'

    def make_hit(qseqid, sseqid, length, sstart=1, send=None, slen=2000):
        if send is None:
            send = sstart + length - 1
        return {
            'qseqid': qseqid, 'sseqid': sseqid, 'pident': 100.0, 'length': length,
            'mismatch': 0, 'gapopen': 0, 'qstart': 1, 'qend': length,
            'sstart': sstart, 'send': send, 'evalue': 0.0, 'bitscore': length * 2.0,
            'qlen': length, 'slen': slen, 'stitle': sseqid, 'qcovs': 100.0,
            'qframe': 1, 'sframe': 1, 'nident': length, 'use_hit': True
        }

    # On p1 (3000bp) and p2 (2000bp), c_unique covers 1..2000 and c_rep covers a non-overlapping region 2001..3000 (or 1..1000 on p2)
    c_unique.blast_data = pd.DataFrame([
        make_hit("c_unique", "p1", 2000, sstart=1, send=2000, slen=3000),
        make_hit("c_unique", "p2", 1000, sstart=1, send=1000, slen=2000)
    ])
    c_rep.blast_data = pd.DataFrame([
        make_hit("c_rep", "p1", 1000, sstart=2001, send=3000, slen=3000),
        make_hit("c_rep", "p2", 1000, sstart=1001, send=2000, slen=2000)
    ])

    p1 = Plasmid("p1", cluster_id="cl1", pling_type="pt1", length=3000, plasmid_mash_db="db.msh")
    p2 = Plasmid("p2", cluster_id="cl2", pling_type="pt2", length=2000, plasmid_mash_db="db.msh")

    p1.add_contigs([c_unique, c_rep])
    p2.add_contigs([c_unique, c_rep])

    assigned_plasmids, assigned_contigs = assign_plasmids(
        {"p1": p1, "p2": p2},
        {"c_unique": c_unique, "c_rep": c_rep},
        config=config
    )

    assert len(assigned_plasmids) == 2
    by_id = {p.id: p for p in assigned_plasmids}
    assert "c_unique" in by_id["p1"].contigs
    assert "c_unique" not in by_id["p2"].contigs  # removed from p2 after p1 claimed it
    assert "c_rep" in by_id["p1"].contigs and "c_rep" in by_id["p2"].contigs  # repetitive contig shared by both
    assert c_rep.plas_data == [("p1", "cl1"), ("p2", "cl2")]


def test_calculate_assignment_confidence_and_conflict_resolution():
    p_strong = Plasmid("p_strong", "c1", "pt1", 5000, "db.msh")
    p_strong.scores.update({
        'blast_similarity': 95.0,
        'synteny_score': 0.9,
        'length_score': 0.95,
        'read_gaps_score': 0.9
    })
    p_strong.mash_score = 0.998
    p_strong.read_cov_data = object()

    p_weak = Plasmid("p_weak", "c2", "pt2", 5000, "db.msh")
    p_weak.scores.update({
        'blast_similarity': 70.0,
        'synteny_score': 0.2,
        'length_score': 0.5,
        'read_gaps_score': 0.0
    })
    p_weak.mash_score = 0.90

    conf_strong = calculate_assignment_confidence(p_strong)
    conf_weak = calculate_assignment_confidence(p_weak)
    assert conf_strong == 1.0
    assert conf_weak < 0.4

    # Test resolve_assignment_conflicts
    contig = Contig("c1", SeqRecord(Seq("ATGC" * 100), id="c1"))
    resolved = resolve_assignment_conflicts(
        {"c1": "p_strong"},
        {"p_strong": p_strong, "p_weak": p_weak},
        {"c1": contig}
    )
    assert resolved["c1"] == "p_strong"


def test_validate_assignments_warnings():
    p1 = Plasmid("p1", "c1", "pt1", 10000, "db.msh")
    c_small = Contig("c_small", SeqRecord(Seq("ATGC" * 250), id="c_small"), length=1000)
    c_small.copy_num = 8.0  # > 5x high copy number

    warnings = validate_assignments({"c_small": "p1"}, {"p1": p1}, {"c_small": c_small})
    assert len(warnings) == 2
    assert any("Large size mismatch" in w for w in warnings)
    assert any("High copy number" in w for w in warnings)

    c_low = Contig("c_low", SeqRecord(Seq("ATGC" * 2500), id="c_low"), length=10000)
    c_low.copy_num = 0.2  # < 0.5x low copy number
    warnings_low = validate_assignments({"c_low": "p1"}, {"p1": p1}, {"c_low": c_low})
    assert len(warnings_low) == 1
    assert "Low copy number" in warnings_low[0]


