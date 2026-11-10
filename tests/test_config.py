"""Tests for accio.config."""

import pytest
from accio.config import AnalysisConfig, get_config, validate_config, FilePatterns


def test_get_config_defaults():
    config = get_config()
    assert isinstance(config, AnalysisConfig)
    assert config.MIN_COVERAGE == 75
    assert config.MIN_IDENTITY == 75
    assert config.FIND_NOVEL_PLASMIDS is True
    assert config.THREADS == 4
    assert config.CLASSIFIER == 'plasme'


def test_validate_config_valid():
    config = get_config()
    assert validate_config(config) is True


def test_validate_config_invalid_coverage():
    config = get_config()
    config.MIN_COVERAGE = 150
    with pytest.raises(ValueError, match="MIN_COVERAGE must be between 0 and 100"):
        validate_config(config)

    config.MIN_COVERAGE = -5
    with pytest.raises(ValueError, match="MIN_COVERAGE must be between 0 and 100"):
        validate_config(config)


def test_validate_config_invalid_identity():
    config = get_config()
    config.MIN_IDENTITY = 105
    with pytest.raises(ValueError, match="MIN_IDENTITY must be between 0 and 100"):
        validate_config(config)

    config.MIN_IDENTITY = -1
    with pytest.raises(ValueError, match="MIN_IDENTITY must be between 0 and 100"):
        validate_config(config)


def test_validate_config_invalid_scoring_weights():
    config = get_config()
    config.SCORING_WEIGHTS_WITH_READS = {'mash_score': 1.5}
    with pytest.raises(ValueError, match="Scoring weights must be between 0 and 1"):
        validate_config(config)


def test_file_patterns_formatting():
    assert FilePatterns.BLAST_OUTPUT.format(sample="sample1") == "sample1_plasmid_blast.tsv"
    assert FilePatterns.PLASMID_DATA_OUTPUT.format(sample="sample1") == "sample1_plasmids_chosen.csv"
    assert FilePatterns.NOVEL_PLASMIDS_OUTPUT.format(sample="sample1") == "sample1_novel_plasmid_contigs.fasta"

