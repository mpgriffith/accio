# 🪄 Accio - Plasmid Analysis Tool

> *"Accio plasmids!"* - Summon plasmids from bacterial genome assemblies with magical precision.

A comprehensive Python package for identifying, analyzing, and assigning plasmids in bacterial genome assemblies using multi-modal similarity scoring and intelligent assignment algorithms.

## ✨ Features

- **🔍 Smart Plasmid Discovery**: Multi-modal scoring combining BLAST, Mash, synteny, and read coverage analysis
- **🧬 Intelligent Assignment**: Greedy algorithm prevents conflicts and ensures biological consistency  
- **🗄️ Database Management**: Create and manage custom plasmid reference databases
- **🔄 Circular Detection**: Enhanced analysis for circular contigs and plasmids
- **📊 Comprehensive Reporting**: Detailed results with confidence scores and validation
- **🛠️ Tool Integration**: Seamless integration with BLAST, Mash, Nucmer, and read mapping tools

## 🚀 Quick Start

### Installation

```bash
# From PyPI (recommended)
pip install accio-plasmids

# From source
git clone https://github.com/yourusername/accio.git
cd accio
pip install -e .
```

### Find Plasmids in Assembly

```bash
# Basic plasmid discovery
accio find assembly.fasta -d /path/to/databases -o results

# With read data for enhanced accuracy
accio find assembly.fasta -r reads_R1.fastq reads_R2.fastq -d databases -o results

# Long read analysis with minimap2
accio find assembly.fasta -r nanopore.fastq --minimap2 -d databases -o results
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

## 📋 Commands

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
- `-t, --threads`: Number of threads (default: 4)

**Examples:**
```bash
# High sensitivity analysis
accio find assembly.fasta -d databases -o results --min_identity 70 --min_score 0.70

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
- `--recursive`: Search directories recursively for FASTA files
- `--deduplicate`: Remove duplicate sequences
- `--cluster_threshold`: Mash similarity for clustering (default: 0.95)
- `--min_length`: Minimum sequence length (default: 1000)
- `--check_circular`: Detect circular sequences

**Examples:**
```bash
# Simple database creation
accio create my_plasmids.fasta -o database

# With quality control and metadata
accio create plasmids/ -o database --recursive --deduplicate --validate_sequences --metadata info.csv

# Update existing database
accio create new_plasmids.fasta -o existing_database --update
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

## 📁 Output Files

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

## 🔧 Prerequisites

### Required External Tools
- **BLAST+** (blastn, makeblastdb)
- **Mash** (mash)
- **MUMmer** (nucmer, show-coords)
- **SAMtools** (samtools)
- **BWA** or **minimap2** (for read mapping)
- **MOB-suite** (plasmid analysis)


### Optional Tools
- **PLASMe** (ML-based plasmid prediction)

Install tools via conda:
```bash
conda install -c bioconda blast mash mummer samtools bwa minimap2
```

### Database Requirements

Accio requires several reference databases:

```
database_dir/
├── plasmidDB.fasta          # Reference plasmid sequences
├── plasmidDB_info.csv       # Plasmid metadata
├── plasmidDB.*              # BLAST database files  
├── plasmidDB.fasta.msh      # Mash sketch
├── PlasmidFinder/           # PlasmidFinder database
└── repetitive_elements.fasta # Repetitive elements
```

## 🧬 Algorithm Overview

### 1. Sequence Analysis
- **BLAST**: Identity and coverage against reference databases
- **Mash**: Rapid k-mer based similarity screening  
- **Nucmer**: Synteny and structural analysis
- **PLASMe**: Machine learning-based classification

### 2. Read Analysis (Optional)
- **Read Mapping**: BWA or minimap2 alignment
- **Coverage**: Depth and copy number estimation
- **Gap Analysis**: Coverage of alignment gaps

### 3. Multi-Modal Scoring
Combines multiple evidence types:
- BLAST similarity (identity × coverage)
- Mash k-mer similarity
- Synteny conservation score
- Length compatibility
- Read coverage support

### 4. Intelligent Assignment
- **Greedy Algorithm**: Assigns highest-scoring plasmids first
- **Conflict Resolution**: Prevents duplicate assignments
- **Cluster Awareness**: Avoids assigning multiple plasmids from same cluster
- **Biological Validation**: Ensures assignments make biological sense

## 🐍 Python API

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

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md).

### Development Setup
```bash
git clone https://github.com/yourusername/accio.git
cd accio
pip install -e ".[dev]"
pytest tests/
```

## 📚 Citation

If you use Accio in your research, please cite:

```bibtex
@software{accio_plasmid_tool,
  title = {Accio: Multi-modal plasmid identification and assignment},
  author = {Your Name},
  url = {https://github.com/yourusername/accio},
  year = {2025}
}
```

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🆘 Support & Help

- **Issues**: [GitHub Issues](https://github.com/yourusername/accio/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/accio/discussions)
- **Documentation**: [Wiki](https://github.com/yourusername/accio/wiki)

## 🏆 Acknowledgments

- BLAST+ team for sequence alignment tools
- Mash developers for rapid similarity screening
- MUMmer team for synteny analysis tools
- BioPython community for sequence handling utilities

---

*🪄 "The magic is not in the wand, but in the algorithm that wields it."*