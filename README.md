# 🪄 Accio - Plasmid Analysis Tool
## (Assembly-based Circular Contig Identification for Outbreaks) 
An Assembly-Based Tool Enabling Plasmid Detection for Hospital Epidemiology 

> *"Accio plasmids!"* - Summon plasmids from bacterial genome assemblies

A bioinformatics tool for assessing the plasmid content of clinical bacterial isolates using input genome assemblies (short-, long-, or hybrid assemblies) using a custom plasmid database.

##  Features

- Plasmid Discovery : Combined scoring combining BLAST, Mash, synteny, and read coverage analysis
- - assign plasmids to contigs in a genome assembly

- Database: Create and manage custom plasmid reference databases
- - Check sequences for replicons and circularity
- - Cluster plasmid sequences using SKANI
- - Group plasmids using Pling
- - Create BLAST, MASH databases with representative plasmid sequences

##  Quick Start

### Installation

```bash

# From source
git clone https://github.com/mpgriffith/accio.git
cd accio
pip install -e .
```

### Find Plasmids in Assembly

```bash
# Basic plasmid discovery
accio find assembly.fasta -d /path/to/database -o results

# With read data for enhanced accuracy
accio find assembly.fasta -r reads_R1.fastq reads_R2.fastq -d database -o results

# Long read analysis with minimap2
accio find assembly.fasta -r nanopore.fastq --minimap2 -d database -o results
```

### Create Custom Database

```bash
# Create database from plasmid sequences
accio create plasmids.fasta -o my_database

# With metadata for enhanced analysis
accio create plasmids.fasta -o my_database --metadata plasmid_info.csv

# From directory of FASTA files
accio create /path/to/plasmid_fastas/ -o my_database --recursive
```

##  Commands

### `accio find` - Plasmid Discovery and Assignment

Identify and assign plasmids in bacterial genome assemblies.

**Basic Usage:**
```bash
accio find <assembly.fasta> -d <database_dir> -o <output_dir> [options]
```

**Key Options:**
- `-r, --reads`: Read files for coverage analysis
- `--minimap2`: Use minimap2 for long read mapping
- `--min_identity`: BLAST identity threshold (default: 75)
- `--min_coverage`: BLAST coverage threshold (default: 75) 
- `--min_score`: Minimum assignment score (default: 0.80)
- `--skip_novel_plasmids`: Skip finding and analyzing novel plasmid contigs (skips PLASMe/geNomad and MOB-recon)
- `-t, --threads`: Number of threads (default: 4)

**Examples:**
```bash
# High sensitivity analysis
accio find assembly.fasta -d databases -o results --min_identity 70 --min_score 0.70

# Skip novel plasmid detection
accio find assembly.fasta -d databases -o results --skip_novel_plasmids

# Paired-end reads with custom thresholds
accio find assembly.fasta -r R1.fastq R2.fastq -d databases -o results -t 8

# Oxford Nanopore long reads
accio find assembly.fasta -r nanopore.fastq --minimap2 --mapping_preset map-ont -d databases -o results
```

### `accio create` - Database Creation

Build custom plasmid reference databases from sequences.

**Basic Usage:**
```bash
accio create <input> -o <database_dir> [options]
```

**Key Options:**
- `--metadata`: CSV file with plasmid metadata
- `--pling_dir`: Path to pre-made pling directory to reuse instead of running pling
- `--containment_distance`: Pling containment distance threshold (default: 0.3)
- `--skip_edit_dcj`: Skip pling `edit_dcj` step and use base typing directly (works on both new pling runs and with `--pling_dir`)
- `--run_edit_dcj`: Force running pling `edit_dcj` (requires `communities.pkl` and distance files)
- `--recursive`: Search directories recursively for FASTA files
- `--min_length`: Minimum sequence length (default: 10000)
- `--check_circular`: Detect and mark circular sequences (default: off)
- `--no_plasmid_check`: Do not filter sequences by circularity / plasmid replicon checks
- `--cluster_threshold`: Skani similarity percentage threshold for clustering (default: 99.95)
- `--no_clustering`: Skip clustering step (skani) (false by default)

**Examples:**
```bash
# Simple database creation
accio create my_plasmids.fasta -o database

# With quality control and metadata
accio create plasmids/ -o database --recursive --deduplicate --validate_sequences --metadata info.csv

# With a pre-made pling directory (runs edit_dcj by default if communities.pkl is present, or falls back to pre-existing typing)
accio create plasmids.fasta -o database --pling_dir /path/to/premade_pling

# With a pre-made pling directory and explicitly skip running edit_dcj
accio create plasmids.fasta -o database --pling_dir /path/to/premade_pling --skip_edit_dcj

# New pling run from scratch, but skipping edit_dcj
accio create plasmids.fasta -o database --skip_edit_dcj

# Update existing database
accio create new_plasmids.fasta -o updated_database --update previous_database
```

### `accio check` - System Verification

Verify tool availability and database integrity.

```bash
# Check external tools
accio check --tools

# Check database files
accio check --databases /path/to/database

# Check everything
accio check --all --databases /path/to/database
```

##  Output Files

### Plasmid Discovery (`accio find`)

- **`{sample}_plasmids_assigned.csv`**: Assigned plasmids with scores and metadata
- **`{sample}_contig_analysis.csv`**: Per-contig analysis results  
- **`{sample}_novel_plasmids.fasta`**: Novel plasmid sequences
- **`{sample}_summary.json`**: Analysis summary and statistics

### Database Creation (`accio create`)

- **`plasmidDB.fasta`**: Combined plasmid sequences
- **`plasmidDB_info.csv`**: Plasmid metadata and clustering information
- **`plasmidDB.*`**: BLAST database files (`.nhr`, `.nin`, `.nsq`)
- **`plasmidDB.fasta.msh`**: Mash sketch database
- **`database_stats.txt`**: Database statistics and summary

##  Prerequisites

### Required External Tools
- **BLAST+** (blastn, makeblastdb)
- **Mash** (mash)
- **MUMmer** (nucmer, show-coords)
- **SAMtools** (samtools)
- **BWA** or **minimap2** (for read mapping)
- **MOB-suite** (plasmid analysis)
- **PLASme** (ML-based plasmid prediction) or **geNomad**
- **skani** (ANI calculation for clustering plasmids)
- **AMRFinder** 


Install tools via conda:
```bash
conda install -c bioconda blast mash mummer samtools bwa minimap2 pysam mob_suite skani amrfinder
```

## 🧬 Algorithm Overview

### 1. Sequence Analysis
- **BLAST**: Identity and coverage against reference databases
- **Mash**: Rapid k-mer based similarity screening  
- **Nucmer**: Synteny and structural analysis
- **PLASMe** or **geNomad**: Machine learning-based classification

### 2. Read Analysis (Optional)
- **Read Mapping**: BWA or minimap2 alignment
- **Coverage**: Depth and copy number estimation
- **Gap Analysis**: Coverage of alignment gaps

### 3. Scoring Scoring
Combines multiple evidence types:
- BLAST similarity (identity × coverage)
- Mash k-mer similarity
- Synteny conservation score
- Length compatibility
- Read coverage support

### 4. Assignment of Plasmids
-  Assigns highest-scoring plasmids first
- allows for repetitive contigs to match to multiple plasmids/multiple regions on a plasmid


##  Python API

```python
from accio import AccioWorkflow, get_config

# Initialize with custom configuration
config = get_config()
config.MIN_IDENTITY = 80
config.MIN_COVERAGE = 80

# Create workflow
workflow = AccioWorkflow(config)

# Find plasmids
results = workflow.find_plasmids(
    assembly_file="assembly.fasta",
    reads=["R1.fastq", "R2.fastq"],
    database_dir="databases",
    output_dir="results"
)

# Access results
for plasmid in results.assigned_plasmids:
    print(f"Assigned: {plasmid.id} (score: {plasmid.score:.3f})")
```

## 📊 Configuration

### Scoring Weights

Customize scoring weights for different data types:

```python
# With read data
config.SCORING_WEIGHTS_WITH_READS = {
    'mash_score': 0.3,
    'read_gaps_score': 0.2,
    'blast_similarity': 0.2,
    'synteny_score': 0.15,
    'length_score': 0.15
}

# Assembly only
config.SCORING_WEIGHTS_NO_READS = {
    'blast_similarity': 0.4,
    'mash_score': 0.3,
    'synteny_score': 0.2,
    'length_score': 0.1
}
```

### Analysis Thresholds

```python
config.MIN_IDENTITY = 75        # BLAST identity threshold
config.MIN_COVERAGE = 75        # BLAST coverage threshold  
config.MIN_PLASMID_SCORE = 0.80 # Assignment score threshold
config.MIN_MASH_SCORE = 0.90    # Mash similarity threshold
```

## Running Tests

Run the test suite using `pytest`:

```bash
# Install test dependencies
pip install -e .[test]

# Run all tests
pytest
```

## Citation
ACCIO is currently available as a preprint, please cite: [Raabe NJ, Griffith MP, Rangachar Srinivasa V, Waggle KD, Sundermann AJ, Pless LL, Snyder GM, Brooks MM, Van Tyne D, Harrison LH. ACCIO: an assembly-based tool enabling plasmid detection. Microb Genom. 2026 Oct;12(10). doi: 10.1099/mgen.0.001843. PMID: 42932221: https://doi.org/10.1101/2025.10.30.25338662](https://www.medrxiv.org/content/10.1101/2025.10.30.25338662).

##  License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

_“The truth." Dumbledore sighed. "It is a beautiful and terrible thing—and should therefore be treated with great caution.”—Harry Potter and the Sorcerer's Stone._