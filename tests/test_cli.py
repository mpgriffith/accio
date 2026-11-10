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

    # Default find has skip_novel_plasmids = False
    default_args = parser.parse_args(['find', 'assembly.fasta', '-d', 'my_db', '-o', 'res'])
    assert default_args.command == 'find'
    assert default_args.skip_novel_plasmids is False

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


def test_main_parser_check_subcommand():
    parser = create_main_parser()
    args = parser.parse_args(['check', '--tools'])
    assert args.command == 'check'
    assert args.tools is True

