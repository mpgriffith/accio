"""Tests for accio.cli argument parsing and commands."""

import argparse
import pytest
from accio.cli import create_main_parser, handle_find_command, handle_create_command, handle_multi_command


def test_main_parser_create_subcommand():
    parser = create_main_parser()
    
    # Test create with --pling_dir
    args = parser.parse_args(['create', 'plasmids.fasta', '-o', 'out_db', '--pling_dir', '/path/to/pling'])
    assert args.command == 'create'
    assert args.input == ['plasmids.fasta']
    assert args.output == 'out_db'
    assert args.pling_dir == '/path/to/pling'
    assert args.metadata is None
    assert args.force is False

    # Test create with kebab-case alias --pling-dir and --metadata
    args2 = parser.parse_args(['create', 'p1.fasta', 'p2.fasta', '-o', 'out_db', '--pling-dir', '/my/pling', '--metadata', 'meta.csv'])
    assert args2.command == 'create'
    assert args2.input == ['p1.fasta', 'p2.fasta']
    assert args2.pling_dir == '/my/pling'
    assert args2.metadata == 'meta.csv'
    assert args2.skip_edit_dcj is False
    assert args2.run_edit_dcj is False

    # Test create with --skip_edit_dcj and aliases
    for flag in ['--skip_edit_dcj', '--skip-edit-dcj', '--no_edit_dcj', '--no-edit-dcj']:
        args_skip = parser.parse_args(['create', 'plasmids.fasta', '-o', 'out_db', flag])
        assert args_skip.skip_edit_dcj is True

    # Test create with --run_edit_dcj and alias
    for flag in ['--run_edit_dcj', '--run-edit-dcj']:
        args_run = parser.parse_args(['create', 'plasmids.fasta', '-o', 'out_db', flag])
        assert args_run.run_edit_dcj is True


def test_main_parser_find_subcommand_skip_novel_plasmids():
    parser = create_main_parser()

    # Default find has skip_novel_plasmids = False, default classifier is genomad, min_blast_similarity is 85, genomad_splits is None
    default_args = parser.parse_args(['find', 'assembly.fasta', '-d', 'my_db', '-o', 'res'])
    assert default_args.command == 'find'
    assert default_args.skip_novel_plasmids is False
    assert default_args.classifier == 'genomad'
    assert default_args.min_blast_similarity == 85.0
    assert default_args.genomad_splits is None

    # Test --min_blast_similarity and --splits
    custom_args = parser.parse_args(['find', 'assembly.fasta', '-d', 'my_db', '-o', 'res', '--min_blast_similarity', '90', '--splits', '8'])
    assert custom_args.min_blast_similarity == 90.0
    assert custom_args.genomad_splits == 8

    # Test --classifier none
    none_args = parser.parse_args(['find', 'assembly.fasta', '-d', 'my_db', '-o', 'res', '--classifier', 'none'])
    assert none_args.classifier == 'none'

    # Test --skip_novel_plasmids flag
    args1 = parser.parse_args(['find', 'assembly.fasta', '-d', 'my_db', '-o', 'res', '--skip_novel_plasmids'])
    assert args1.skip_novel_plasmids is True

    # Test aliases
    args2 = parser.parse_args(['find', 'assembly.fasta', '-d', 'my_db', '-o', 'res', '--skip-novel-plasmids'])
    assert args2.skip_novel_plasmids is True

    args3 = parser.parse_args(['find', 'assembly.fasta', '-d', 'my_db', '-o', 'res', '--skip_novel'])
    assert args3.skip_novel_plasmids is True

    args4 = parser.parse_args(['find', 'assembly.fasta', '-d', 'my_db', '-o', 'res', '--skip-novel'])
    assert args4.skip_novel_plasmids is True


def test_main_parser_multi_subcommand():
    parser = create_main_parser()
    args = parser.parse_args(['multi', '--assemblies', 'asm_dir', '-d', 'my_db', '-o', 'out_dir', '--skip_novel_plasmids'])
    assert args.command == 'multi'
    assert args.assemblies == ['asm_dir']
    assert args.skip_novel_plasmids is True

    # Test -a and --output aliases
    args2 = parser.parse_args(['multi', '-a', 'a1.fa', 'a2.fa', '-d', 'my_db', '--output', 'out_dir'])
    assert args2.assemblies == ['a1.fa', 'a2.fa']
    assert args2.output_dir == 'out_dir'


def test_main_parser_check_subcommand():
    parser = create_main_parser()
    args = parser.parse_args(['check', '--tools'])
    assert args.command == 'check'
    assert args.tools is True

    # Test -d alias
    args2 = parser.parse_args(['check', '-d', 'my_db'])
    assert args2.databases == 'my_db'


def test_main_parser_create_options_batch1():
    parser = create_main_parser()

    # Default min_length is 10000, check_circular is False, containment_distance is 0.3, cluster_threshold is 99.95
    args_default = parser.parse_args(['create', 'plasmids.fasta', '-o', 'out_db'])
    assert args_default.min_length == 10000
    assert args_default.check_circular is False
    assert args_default.containment_distance == 0.3
    assert args_default.cluster_threshold == 99.95
    assert args_default.recursive is False
    assert args_default.validate_sequences is False
    assert args_default.no_plasmid_check is False

    # Explicit flags and aliases
    args = parser.parse_args([
        'create', 'input_dir', '-o', 'out_db',
        '--recursive', '--validate_sequences', '--min_length', '5000',
        '--check_circular', '--containment_distance', '0.25', '--no_plasmid_check'
    ])
    assert args.recursive is True
    assert args.validate_sequences is True
    assert args.min_length == 5000
    assert args.check_circular is True
    assert args.containment_distance == 0.25
    assert args.no_plasmid_check is True


def test_main_parser_query_subcommand(tmp_path):
    parser = create_main_parser()
    args = parser.parse_args(['query', 'my_db', '--plasmid_id', 'p123', '--pling_type', 'type_1', '-o', 'out.csv'])
    assert args.command == 'query'
    assert args.database == 'my_db'
    assert args.plasmid_id == 'p123'
    assert args.pling_type == 'type_1'
    assert args.output == 'out.csv'

    # Test query with -d flag and replicon alias
    args2 = parser.parse_args(['query', '-d', 'my_db', '--replicon', 'IncFIB'])
    assert args2.database_opt == 'my_db'
    assert args2.replicon_type == 'IncFIB'


def test_handle_query_command(tmp_path):
    import pandas as pd
    from accio.cli import handle_query_command
    import argparse

    db_dir = tmp_path / "mock_db"
    db_dir.mkdir()
    df = pd.DataFrame({
        "plasmid_id": ["p1", "p2", "p3"],
        "length": [5000, 20000, 45000],
        "cluster_id": ["c1", "c1", "c2"],
        "rep_types": ["IncFIB", "IncI1", "ColRNAI"],
        "pling_type": ["type_A", "type_A", "type_B"],
        "community": ["comm_1", "comm_1", "comm_2"]
    })
    df.to_csv(db_dir / "plasmidDB_info.csv", index=False)

    # 1. Query by plasmid_id
    args1 = argparse.Namespace(database=str(db_dir), database_opt=None, plasmid_id="p1", pling_type=None, community=None, replicon_type=None, output=None)
    assert handle_query_command(args1) == 0

    # 2. Query by pling_type with output CSV and -d flag
    out_csv = tmp_path / "query_out.csv"
    args2 = argparse.Namespace(database=None, database_opt=str(db_dir), plasmid_id=None, pling_type="type_A", community=None, replicon_type=None, output=str(out_csv))
    assert handle_query_command(args2) == 0
    res_df = pd.read_csv(out_csv)
    assert len(res_df) == 2
    assert set(res_df["plasmid_id"]) == {"p1", "p2"}

    # 3. Query with non-existent database returns 1
    args_bad = argparse.Namespace(database=str(tmp_path / "non_existent"), database_opt=None, plasmid_id="p1", pling_type=None, community=None, replicon_type=None, output=None)
    assert handle_query_command(args_bad) == 1


def test_main_module_import():
    import accio.__main__
    assert hasattr(accio.__main__, 'main')


def test_essential_tool_check_in_cli(monkeypatch, tmp_path):
    """Test that missing tools with available=False correctly fail find and create commands."""
    mock_availability = {
        'blastn': {'available': False, 'install_hint': 'conda install blast'},
        'mash': {'available': True, 'install_hint': 'conda install mash'},
        'nucmer': {'available': True, 'install_hint': 'conda install mummer'},
        'show-coords': {'available': True, 'install_hint': 'conda install mummer'},
        'samtools': {'available': True, 'install_hint': 'conda install samtools'},
        'makeblastdb': {'available': False, 'install_hint': 'conda install blast'},
    }
    monkeypatch.setattr("accio.cli.check_tool_availability", lambda: mock_availability)

    dummy_fa = tmp_path / "dummy.fa"
    dummy_fa.write_text(">c1\nATGC\n")

    # find command
    find_args = argparse.Namespace(
        assembly=str(dummy_fa),
        database=str(tmp_path),
        output=str(tmp_path),
        reads=None,
        minimap2=False,
        log_file=None
    )
    assert handle_find_command(find_args) == 1

    # create command
    create_args = argparse.Namespace(
        input=[str(dummy_fa)],
        output=str(tmp_path),
        force=True,
        threads=4,
        min_length=10000,
        max_length=1000000
    )
    assert handle_create_command(create_args) == 1


def test_multi_thread_division(monkeypatch, tmp_path):
    """Test that handle_multi_command divides total threads among parallel workers."""
    # Create 4 dummy assembly files
    sample_dir = tmp_path / "assemblies"
    sample_dir.mkdir()
    for i in range(4):
        (sample_dir / f"sample_{i}.fasta").write_text(f">s{i}\nATGC\n")

    captured_jobs = []

    def mock_map(func, jobs):
        captured_jobs.extend(jobs)
        return [0] * len(jobs)

    class MockPool:
        def __init__(self, processes):
            self.processes = processes
        def __enter__(self):
            return self
        def __exit__(self, *args):
            pass
        def map(self, func, jobs):
            return mock_map(func, jobs)

    monkeypatch.setattr("multiprocessing.Pool", MockPool)
    monkeypatch.setattr("accio.cli.handle_combine_command", lambda args: 0)

    multi_args = argparse.Namespace(
        assemblies=str(sample_dir),
        assembly=None,
        database=str(tmp_path),
        output_dir=str(tmp_path / "multi_out"),
        reads_dir=None,
        threads=16
    )

    handle_multi_command(multi_args)
    assert len(captured_jobs) == 4
    # With 16 threads and 4 jobs, num_parallel = 16 // 8 = 2, so threads_per_sample = 16 // 2 = 8
    for job in captured_jobs:
        assert job.threads == 8


def test_cli_setup_mob_parser():
    """Test setup-mob parser arguments."""
    parser = create_main_parser()
    args = parser.parse_args(['setup-mob', '--env-name', 'custom_mob'])
    assert args.command == 'setup-mob'
    assert args.env_name == 'custom_mob'


def test_cli_find_and_create_mob_env_args():
    """Test --mob_env and --mob-bin-dir in find and create parsers."""
    parser = create_main_parser()
    args_find = parser.parse_args(['find', 'asm.fasta', '-d', 'db', '-o', 'out', '--mob_env', 'test_env', '--mob_bin_dir', '/path/to/bins'])
    assert args_find.mob_env == 'test_env'
    assert args_find.mob_bin_dir == '/path/to/bins'

    args_create = parser.parse_args(['create', 'p1.fasta', '-o', 'out_db', '--mob-env', 'create_env', '--mob-bin-dir', '/other/bins'])
    assert args_create.mob_env == 'create_env'
    assert args_create.mob_bin_dir == '/other/bins'


def test_handle_setup_mob_command_runs_without_name_error(monkeypatch):
    """Test that handle_setup_mob_command has shutil and subprocess imported and executes cleanly."""
    from accio.cli import handle_setup_mob_command
    called_cmd = []
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mamba" if name == "mamba" else None)
    monkeypatch.setattr("subprocess.check_call", lambda cmd: called_cmd.append(cmd))

    rc = handle_setup_mob_command(argparse.Namespace(env_name="mob_suite_test"))
    assert rc == 0
    assert len(called_cmd) == 1
    assert "mob_suite_test" in called_cmd[0]


def test_handle_combine_command_basic_and_db_remapping(tmp_path):
    """Test accio collate (handle_combine_command) both with and without database pling_type remapping."""
    from accio.cli import handle_combine_command
    import pandas as pd

    # Set up two sample output directories
    db_old = tmp_path / "db_v1"
    db_new = tmp_path / "db_v2"
    db_old.mkdir()
    db_new.mkdir()

    # Collation database (db_v2) maps old pling_type -> new pling_type
    pd.DataFrame([
        {'plasmid_id': 'p1', 'pling_type': 'new_subcomm_1', 'db_v1_pling_type': 'old_subcomm_1'},
        {'plasmid_id': 'p2', 'pling_type': 'new_subcomm_2', 'db_v1_pling_type': 'old_subcomm_2'},
    ]).to_csv(db_new / "plasmidDB_info.csv", index=False)

    s1_dir = tmp_path / "results" / "sampleA"
    s2_dir = tmp_path / "results" / "sampleB"
    s1_dir.mkdir(parents=True)
    s2_dir.mkdir(parents=True)

    pd.DataFrame([
        {'ACCIO_Plasmids': 'chosen_plasmid', 'plasmid_id': 'p1', 'cluster_id': 'c1', 'pling_type': 'old_subcomm_1'}
    ]).to_csv(s1_dir / "sampleA_plasmids_chosen.csv", index=False)
    (s1_dir / "sampleA_accio.log").write_text(f"INFO - Using database: {db_old}\n")

    pd.DataFrame([
        {'ACCIO_Plasmids': 'chosen_plasmid', 'plasmid_id': 'p2', 'cluster_id': 'c2', 'pling_type': 'old_subcomm_2'}
    ]).to_csv(s2_dir / "sampleB_plasmids_chosen.csv", index=False)
    (s2_dir / "sampleB_accio.log").write_text(f"INFO - Using database: {db_old}\n")

    out_csv = tmp_path / "collated_summary.csv"
    collate_args = argparse.Namespace(
        input_pattern=str(tmp_path / "results" / "*" / "*_plasmids_chosen.csv"),
        output=str(out_csv),
        database=str(db_new)
    )
    rc = handle_combine_command(collate_args)
    assert rc == 0
    assert out_csv.exists()

    summary_df = pd.read_csv(out_csv)
    assert len(summary_df) == 2
    assert set(summary_df['Sample']) == {'sampleA', 'sampleB'}
    assert set(summary_df['pling_type']) == {'new_subcomm_1', 'new_subcomm_2'}


def test_cli_version_flag():
    """Test that -v / --version prints version and exits with code 0."""
    parser = create_main_parser()
    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args(['--version'])
    assert exc_info.value.code == 0



