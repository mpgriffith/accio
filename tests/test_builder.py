"""Tests for DatabaseBuilder in accio.core.builder."""

import pytest
import os
import pandas as pd
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from accio.core.builder import DatabaseBuilder


@pytest.fixture
def builder():
    return DatabaseBuilder(threads=1)


def test_clean_sequence_id(builder):
    assert builder._clean_sequence_id("seq1|plasmid:A;B,C extra") == "seq1_plasmid_A_B_C"
    assert builder._clean_sequence_id("simple_id") == "simple_id"


def test_validate_sequence(builder):
    valid_rec = SeqRecord(Seq("ATGC" * 100), id="valid")
    assert builder._validate_sequence(valid_rec) is True

    # High N content (>5%)
    n_rec = SeqRecord(Seq("ATGC" * 10 + "NNNNNNNNNNNNNNNNNNNN"), id="high_n")
    assert builder._validate_sequence(n_rec) is False

    # Invalid bases
    invalid_rec = SeqRecord(Seq("ATGC123XYZ!"), id="invalid")
    assert builder._validate_sequence(invalid_rec) is False


def test_check_circularity_description(builder):
    rec = SeqRecord(Seq("ATGC" * 100), id="c1", description="assembly circular=true contig")
    assert builder._check_circularity(rec) is True

    rec_linear = SeqRecord(Seq("ATGC" * 100), id="c2", description="linear contig")
    assert builder._check_circularity(rec_linear) is False


def test_check_circularity_terminal_overlap(builder):
    overlap = "ATGCATGCATGCATGCATGC"  # 20 bp
    middle = "AAAA" * 50
    seq = overlap + middle + overlap
    rec = SeqRecord(Seq(seq), id="c_overlap", description="assembly contig")
    assert builder._check_circularity(rec) is True


def test_deduplicate_sequences(builder):
    rec1 = SeqRecord(Seq("ATGC" * 50), id="seq1")
    rec2 = SeqRecord(Seq("ATGC" * 50), id="seq2")  # duplicate sequence
    rec3 = SeqRecord(Seq("GGCC" * 50), id="seq3")

    seqs = {"seq1": rec1, "seq2": rec2, "seq3": rec3}
    deduped = builder._deduplicate_sequences(seqs)
    assert len(deduped) == 2
    assert "seq1" in deduped
    assert "seq3" in deduped
    assert "seq2" not in deduped


def test_create_plasmid_communities_premade_dir_missing(builder, tmp_path):
    output_dir = tmp_path / "db_out"
    non_existent = tmp_path / "does_not_exist_pling"
    meta_df = pd.DataFrame(index=["p1"], columns=["plasmid_id"], data=[["p1"]])

    with pytest.raises(FileNotFoundError, match="Provided pling directory does not exist"):
        builder._create_plasmid_communities(
            sequences={},
            metadata_df=meta_df,
            output_dir=str(output_dir),
            premade_pling_dir=str(non_existent)
        )


def test_create_plasmid_communities_premade_dir_success(builder, tmp_path):
    # Setup pre-made pling directory structure
    premade_dir = tmp_path / "my_premade_pling"
    premade_dir.mkdir(parents=True)
    
    # 1. typing.tsv in pling_dcj_length_subcommunities/objects/
    subcomm_objects = premade_dir / "pling_dcj_length_subcommunities" / "objects"
    subcomm_objects.mkdir(parents=True)
    typing_df = pd.DataFrame({"plasmid": ["p1", "p2"], "typing": ["type_A", "type_B"]})
    typing_df.to_csv(subcomm_objects / "typing.tsv", sep="\t", index=False)

    # 2. communities.tsv in containment/containment_communities/objects/
    comm_objects = premade_dir / "containment" / "containment_communities" / "objects"
    comm_objects.mkdir(parents=True)
    communities_df = pd.DataFrame({"plasmid": ["p1", "p2"], "comm": ["comm_1", "comm_2"]})
    communities_df.to_csv(comm_objects / "communities.tsv", sep="\t", index=False)

    # 3. hub_plasmids.csv in dcj_thresh_4_graph/objects/
    dcj_objects = premade_dir / "dcj_thresh_4_graph" / "objects"
    dcj_objects.mkdir(parents=True)
    hub_df = pd.DataFrame(index=["p1"], data={"hub": [True]})
    hub_df.to_csv(dcj_objects / "hub_plasmids.csv")

    output_dir = tmp_path / "db_out"
    meta_df = pd.DataFrame(
        index=["p1", "p2"],
        data={
            "plasmid_id": ["p1", "p2"],
            "amr_genes": ["blaTEM-1", "tet(A),blaTEM-1"]
        }
    )

    result_df = builder._create_plasmid_communities(
        sequences={},
        metadata_df=meta_df,
        output_dir=str(output_dir),
        premade_pling_dir=str(premade_dir)
    )

    # Check that directory was copied to output_dir / pling
    target_pling = output_dir / "pling"
    assert target_pling.exists()
    assert (target_pling / "pling_dcj_length_subcommunities" / "objects" / "typing.tsv").exists()

    # Check metadata was annotated with communities and pling_types
    assert "pling_type" in result_df.columns
    assert "community" in result_df.columns
    assert "subcommunity_amr_genes" in result_df.columns
    # p1 was a hub plasmid, so it gets hub suffix
    assert result_df.loc["p1", "pling_type"] == "comm_1_hub0"


def test_create_plasmid_communities_flat_csv_data(builder, tmp_path):
    # Minimal flat directory layout with typing.csv, all_plasmid_distances.csv, hub_plasmids.csv
    premade_dir = tmp_path / "flat_pling"
    premade_dir.mkdir(parents=True)

    # 1. typing.csv
    typing_df = pd.DataFrame({"plasmid": ["p1", "p2"], "typing": ["type_1", "type_2"]})
    typing_df.to_csv(premade_dir / "typing.csv", index=False)

    # 2. hub_plasmids.csv
    hub_df = pd.DataFrame(index=["p1"], data={"hub": [True]})
    hub_df.to_csv(premade_dir / "hub_plasmids.csv")

    # 3. all_plasmid_distances.csv
    dist_df = pd.DataFrame({
        "plasmid_1": ["p1"],
        "plasmid_2": ["p2"],
        "distance": [0.02]
    })
    dist_df.to_csv(premade_dir / "all_plasmid_distances.csv", index=False)

    output_dir = tmp_path / "flat_db_out"
    meta_df = pd.DataFrame(
        index=["p1", "p2"],
        data={
            "plasmid_id": ["p1", "p2"],
            "amr_genes": ["blaTEM-1", "tet(A)"]
        }
    )

    result_df = builder._create_plasmid_communities(
        sequences={},
        metadata_df=meta_df,
        output_dir=str(output_dir),
        premade_pling_dir=str(premade_dir)
    )

    assert "pling_type" in result_df.columns
    # p1 is in hub_plasmids, so it receives hub suffix
    assert result_df.loc["p1", "pling_type"] == "type_1_hub0"
    assert result_df.loc["p2", "pling_type"] == "type_2"

    # Verify all_plasmids_distances.tsv was normalized
    assert (output_dir / "pling" / "all_plasmids_distances.tsv").exists()


def test_load_previous_pling_data_and_annotate(builder, tmp_path):
    pling_dir = tmp_path / "test_load_pling"
    pling_dir.mkdir()

    typing_df = pd.DataFrame({"plasmid_id": ["p10", "p20"], "typing": ["T1", "T2"]})
    typing_df.to_csv(pling_dir / "typing.tsv", sep="\t", index=False)

    hub_df = pd.DataFrame(index=["p10"], data={"hub": [True]})
    hub_df.to_csv(pling_dir / "hub_plasmids.tsv", sep="\t")

    pling_types, hub_plasmids, communities = builder._load_previous_pling_data(pling_dir)

    assert not pling_types.empty
    assert "pling_type" in pling_types.columns
    assert len(hub_plasmids) == 1
    assert communities.empty

    meta_df = pd.DataFrame(
        index=["p10", "p20"],
        data={"plasmid_id": ["p10", "p20"], "amr_genes": ["blaKPC-2", "blaOXA-48"]}
    )

    annotated = builder._annotate_metadata_with_pling(
        meta_df, pling_types, hub_plasmids, communities
    )

    assert annotated.loc["p10", "pling_type"] == "T1_hub0"
    assert annotated.loc["p20", "pling_type"] == "T2"
    assert annotated.loc["p10", "subcommunity_amr_genes"] == "blaKPC-2"


def test_create_plasmid_communities_missing_comm_pkl_handled_gracefully(builder, tmp_path):
    # Directory has distances but NO typing and NO communities.pkl
    premade_dir = tmp_path / "no_pkl_pling"
    premade_dir.mkdir()
    (premade_dir / "all_plasmids_distances.tsv").write_text("plasmid_1\tplasmid_2\tdistance\np1\tp2\t0.1\n")

    output_dir = tmp_path / "no_pkl_db_out"
    meta_df = pd.DataFrame(index=["p1", "p2"], data={"plasmid_id": ["p1", "p2"]})

    # Should not crash with ExternalToolError or FileNotFoundError trying to call plasnet/edit_dcj
    result_df = builder._create_plasmid_communities(
        sequences={},
        metadata_df=meta_df,
        output_dir=str(output_dir),
        premade_pling_dir=str(premade_dir)
    )
    assert isinstance(result_df, pd.DataFrame)


def test_create_plasmid_communities_run_edit_dcj_explicit_flags(builder, tmp_path):
    premade_dir = tmp_path / "explicit_flags_pling"
    premade_dir.mkdir()
    (premade_dir / "typing.csv").write_text("plasmid,typing\np1,type_existing\n")
    (premade_dir / "all_plasmids_distances.tsv").write_text("plasmid_1\tplasmid_2\tdistance\np1\tp2\t0.1\n")
    meta_df = pd.DataFrame(index=["p1"], data={"plasmid_id": ["p1"]})
    output_dir = tmp_path / "explicit_out"

    # 1. run_edit_dcj=True should raise FileNotFoundError because communities.pkl is missing
    with pytest.raises(FileNotFoundError, match="Cannot run edit_dcj"):
        builder._create_plasmid_communities(
            sequences={},
            metadata_df=meta_df,
            output_dir=str(output_dir),
            premade_pling_dir=str(premade_dir),
            run_edit_dcj=True
        )

    # 2. run_edit_dcj=False should succeed and use pre-existing typing without error
    res_skip = builder._create_plasmid_communities(
        sequences={},
        metadata_df=meta_df,
        output_dir=str(output_dir),
        premade_pling_dir=str(premade_dir),
        run_edit_dcj=False
    )
    assert res_skip.loc["p1", "pling_type"] == "type_existing"


def test_create_plasmid_communities_new_pling_skips_edit_dcj(builder, tmp_path, monkeypatch):
    """Verify that when premade_pling_dir is None and run_edit_dcj=False, edit_dcj is skipped."""
    output_dir = tmp_path / "new_pling_out"
    pling_dir = output_dir / "pling"
    dcj_graph_objects = pling_dir / "dcj_thresh_4_graph" / "objects"
    dcj_graph_objects.mkdir(parents=True)

    # Fake typing from pling cluster
    (dcj_graph_objects / "typing.tsv").write_text("plasmid\ttyping\np1\tnew_type_1\n")
    (dcj_graph_objects / "hub_plasmids.csv").write_text("plasmid\np1\n")

    seq1 = SeqRecord(Seq("ATGC" * 100), id="p1")
    meta_df = pd.DataFrame(index=["p1"], data={"plasmid_id": ["p1"]})

    # Mock run_pling to do nothing (files simulated above)
    edit_dcj_called = []
    from accio.wrappers import PlingRunner
    monkeypatch.setattr(PlingRunner, "run_pling", lambda self, *args, **kwargs: None)
    monkeypatch.setattr(PlingRunner, "edit_dcj", lambda self, *args, **kwargs: edit_dcj_called.append(True))

    res = builder._create_plasmid_communities(
        sequences=[seq1],
        metadata_df=meta_df,
        output_dir=str(output_dir),
        premade_pling_dir=None,
        run_edit_dcj=False
    )

    # edit_dcj must NOT have been called
    assert len(edit_dcj_called) == 0
    # Pling type loaded from dcj_thresh_4_graph typing.tsv
    assert "new_type_1" in res.loc["p1", "pling_type"]


def test_pling_runner_edit_dcj_handles_length_and_size_columns(tmp_path, monkeypatch):
    """Verify PlingRunner.edit_dcj safely reads length/size without KeyError."""
    from accio.wrappers import PlingRunner
    runner = PlingRunner()
    runner.output_dir = str(tmp_path)

    # Write dummy all_plasmids_distances.tsv
    (tmp_path / "all_plasmids_distances.tsv").write_text(
        "plasmid_1\tplasmid_2\tdistance\np1\tp2\t1.0\n"
    )

    # Mock run_plasnet to return a path without running subprocess
    monkeypatch.setattr(runner, "run_plasnet", lambda *args, **kwargs: str(tmp_path / "out"))

    # metadata dataframe only has 'length', NOT 'size'
    meta_df = pd.DataFrame(
        index=["p1", "p2"],
        data={"plasmid_id": ["p1", "p2"], "length": [50000, 60000]}
    )

    res = runner.edit_dcj(meta_df)
    assert res == str(tmp_path / "out")
    # Distance by length was created
    assert (tmp_path / "pling_dcj_dists_by_length.tsv").exists()




def test_write_database_files(builder, tmp_path):
    output_dir = tmp_path / "write_db_test"
    seq1 = SeqRecord(Seq("ATGC" * 100), id="seq1")
    sequences = {"seq1": seq1}
    metadata_df = pd.DataFrame([{"plasmid_id": "seq1", "cluster_id": "c1", "length": 400}])

    builder._write_database_files(sequences, metadata_df, str(output_dir))

    assert (output_dir / "plasmidDB.fasta").exists()
    assert (output_dir / "plasmidDB_info.csv").exists()
    assert (output_dir / "database_stats.txt").exists()


def test_repetitive_db_nhr_detection(tmp_path, monkeypatch):
    """Test that download_repetitive_db detects repetitive.dna.fas.nhr and does not rebuild."""
    builder = DatabaseBuilder(threads=4)
    rep_fasta = tmp_path / 'repetitive.dna.fas'
    rep_nhr = tmp_path / 'repetitive.dna.fas.nhr'
    rep_fasta.write_text(">rep\nATCG\n")
    rep_nhr.write_text("dummy nhr")

    called = []
    monkeypatch.setattr(builder.blast_runner, 'run_makeblastdb', lambda *args, **kwargs: called.append(True))
    builder.download_repetitive_db(tmp_path)
    assert len(called) == 0, "run_makeblastdb should not have been called when repetitive.dna.fas.nhr exists"


def test_builder_threads_propagation():
    """Test that DatabaseBuilder propagates threads to all internal runners."""
    builder = DatabaseBuilder(threads=8)
    assert builder.threads == 8
    assert builder.config.THREADS == 8
    assert builder.blast_runner.config.THREADS == 8
    assert builder.mash_runner.config.THREADS == 8
    assert builder.mob_runner.config.THREADS == 8
    assert builder.skani_runner.config.THREADS == 8
    assert builder.pling_runner.config.THREADS == 8
    assert builder.amr_runner.config.THREADS == 8

