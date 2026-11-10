"""Tests for accio.cli argument parsing and commands."""

import pytest
from accio.cli import create_main_parser


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

    # Default find has skip_novel_plasmids = False and default classifier is genomad
    default_args = parser.parse_args(['find', 'assembly.fasta', '-d', 'my_db', '-o', 'res'])
    assert default_args.command == 'find'
    assert default_args.skip_novel_plasmids is False
    assert default_args.classifier == 'genomad'

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

    # Default min_length is 10000, check_circular is False, containment_distance is 0.3
    args_default = parser.parse_args(['create', 'plasmids.fasta', '-o', 'out_db'])
    assert args_default.min_length == 10000
    assert args_default.check_circular is False
    assert args_default.containment_distance == 0.3
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

