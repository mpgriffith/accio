"""Configuration constants and default parameters for plasmid analysis."""

from typing import List, Dict, Any

# BLAST output format and field definitions
DEFAULT_OUTFMT = '"6 std qlen slen stitle qcovs qframe sframe nident"'
DEFAULT_FIELDS = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch',
    'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
]
EXTRA_FIELDS = ['qlen', 'slen', 'stitle', 'qcovs', 'qframe', 'sframe', 'nident']

# Nucmer output column definitions
NUCMER_COORD_COLS = [
    'sstart', 'send', 'qstart', 'qend', 
    'ref_aln_len', 'query_aln_len', 'pident', 'slen', 'qlen', 
    'scov', 'qcovs', 'sframe', 'qframe',
    'sseqid', 'qseqid'
]

NUCMER_SNP_COLS = [
    'ref_pos', 'ref_nucl', 'query_nucl', 'query_pos', 'buff', 'dist_from_end',
    'ref_repeat', 'query_repeat', 'ref_strand', 'query_strand',
    'ref', 'query'
]

MOB_TYPER_COLS = ['sample_id', 'num_contigs', 'size', 'gc', 'md5',
 'rep_type(s)', 'rep_type_accession(s)', 'relaxase_type(s)', 'relaxase_type_accession(s)',
 'mpf_type', 'mpf_type_accession(s)','orit_type(s)', 'orit_accession(s)', 'predicted_mobility',
 'mash_nearest_neighbor', 'mash_neighbor_distance', 'mash_neighbor_identification',
 'primary_cluster_id', 'secondary_cluster_id',
 'predicted_host_range_overall_rank', 'predicted_host_range_overall_name',
 'observed_host_range_ncbi_rank', 'observed_host_range_ncbi_name',
 'reported_host_range_lit_rank', 'reported_host_range_lit_name',
 'associated_pmid(s)']

# Analysis parameters
MAX_PLASMID_LENGTH = 300000

# Default thresholds and parameters
class AnalysisConfig:
    """Configuration class containing all analysis parameters."""
    
    # BLAST filtering parameters
    MIN_COVERAGE = 75
    MIN_IDENTITY = 75
    
    # Repetitive element detection
    REP_MIN_COV = 80
    REP_MIN_IDENT = 80
    
    # Overlap filtering
    OVERLAP_THRESHOLD = 0.95
    
    # Synteny analysis
    GAP_THRESHOLD = 200
    BIN_SIZE = 500
    
    # Copy number thresholds
    HIGH_COPY_THRESHOLD = 1.15
    
    # Scoring weights (with read data)
    SCORING_WEIGHTS_WITH_READS = {
        'mash_score': 0.24,
        'read_gaps_score': 0.22,
        'length_score': 0.17,
        'blast_coverage': 0.16,
        'synteny_score': 0.12,
        'blast_identity': 0.09
    }
    
    # Scoring weights (without read data)
    SCORING_WEIGHTS_NO_READS = {
        'mash_score': 0.26,
        'read_gaps_score': 0,
        'length_score': 0.24,
        'blast_coverage': 0.20,
        'synteny_score': 0.17,
        'blast_identity': 0.12
    }
    
    # Assignment thresholds
    MIN_BLAST_SIMILARITY = 85
    MIN_PLASMID_SCORE = 0.80
    MIN_MASH_SCORE = 0.990
    
    # Use Pling community counts for choosing best plasmid
    USE_PLING_COMMUNITY_COUNTS = False

    # External tool parameters
    THREADS = 4
    FASTANI_FRAG_LEN = 500
    NUCMER_DIAGDIFF = 20
    NUCMER_BREAKLEN = 500
    MIN_READ_DEPTH_FACTOR = 0.25
    DCJ_THRESHOLD=0.05
    DCJ_THRESHOLD_LEN = 120000
    DCJ_THRESHOLD_NUM = 4
    
    OVERWRITE_FILES = False

# File naming patterns
class FilePatterns:
    """Standard file naming patterns for outputs."""
    
    BLAST_OUTPUT = "{sample}_plasmid_blast.tsv"
    PLASMIDFINDER_OUTPUT = "{sample}_plasmidfinder.tsv"
    REP_BLAST_OUTPUT = "{sample}_rep_blast.tsv"
    MASH_OUTPUT = "{sample}_plasmid_mash.csv"
    NUCMER_PREFIX = "{sample}_vs_plasmidDB"
    FASTANI_OUTPUT = "{sample}_fastani.csv"
    CONTIG_COV_OUTPUT = "{sample}_contig_cov.tsv"
    PLASMID_DATA_OUTPUT = "{sample}_plasmids_chosen.csv"
    CONTIG_DATA_OUTPUT = "{sample}_contig_data.csv"
    NOVEL_PLASMIDS_OUTPUT = "{sample}_novel_plasmid_contigs.fasta"

# Database file requirements
class DatabaseFiles:
    """Expected database files and their purposes."""
    
    PLASMID_DB = "plasmidDB.fasta"
    PLASMIDFINDER_DB = "PlasmidFinder"
    RESFINDER_DB = "ResFinder"
    MOB_DB = "mob_cluster/references_updated.fasta"
    REP_DB = "repetitive.dna.fas"
    PLASMID_INFO = "plasmidDB_info.csv"
    MASH_DB = "plasmidDB.fasta.msh"

def get_config() -> AnalysisConfig:
    """Get the default analysis configuration."""
    return AnalysisConfig()

def validate_config(config: AnalysisConfig) -> bool:
    """Validate configuration parameters."""
    if config.MIN_COVERAGE < 0 or config.MIN_COVERAGE > 100:
        raise ValueError("MIN_COVERAGE must be between 0 and 100")
    
    if config.MIN_IDENTITY < 0 or config.MIN_IDENTITY > 100:
        raise ValueError("MIN_IDENTITY must be between 0 and 100")
    
    if not all(0 <= w <= 1 for w in config.SCORING_WEIGHTS_WITH_READS.values()):
        raise ValueError("Scoring weights must be between 0 and 1")
    
    return True