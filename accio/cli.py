"""Accio - Command line interface for plasmid analysis and database creation."""

import argparse
import sys
import os
import logging
import glob
import multiprocessing
from pathlib import Path
from typing import Optional, List
import pandas as pd
from .config import AnalysisConfig, get_config, validate_config
from .wrappers.external_tools import check_tool_availability, validate_databases
from .core.workflow import PlasmidAnalysisWorkflow
from .core.builder import DatabaseBuilder


def setup_logging(log_level: str = 'INFO', log_file: Optional[str] = None) -> None:
    """Set up logging configuration."""
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {log_level}')
    
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    if log_file:
        logging.basicConfig(
            level=numeric_level,
            format=log_format,
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
    else:
        logging.basicConfig(level=numeric_level, format=log_format)


def create_main_parser() -> argparse.ArgumentParser:
    """Create the main argument parser."""
    parser = argparse.ArgumentParser(
        prog='accio',
        description='🪄 Accio - Summon plasmids from bacterial genome assemblies',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Available commands:
    find        Find and assign plasmids in genome assemblies
    create      Create plasmid databases from reference sequences
    check       Check tool and database availability
    collate     Combine results from multiple genomes
    query       Get info about a plasmid in a database
    multi       Find and assign plasmids in many genomes 

Examples:
    # Find plasmids in an assembly
    accio find assembly.fasta -d /path/to/databases -o output_dir
    
    # Create a new database from plasmid sequences
    accio create plasmids.fasta -o database_dir --metadata plasmid_info.csv
    
    # Check system requirements
    accio check --tools --databases /path/to/databases
        """
    )
    
    # Global options
    parser.add_argument('--log_level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='Logging level (default: INFO)'
    )
    
    parser.add_argument('--log_file', help='Log file path (default: stdout only)')
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Add subcommands
    add_find_parser(subparsers)
    add_create_parser(subparsers)
    add_check_parser(subparsers)
    add_collate_parser(subparsers)
    add_multi_parser(subparsers)
    add_query_parser(subparsers)

    return parser


def add_find_parser(subparsers) -> None:
    """Add the 'find' subcommand parser."""
    find_parser = subparsers.add_parser(
        'find',
        help='Find and assign plasmids in genome assemblies',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='🔍 Find plasmids in bacterial genome assemblies using multi-modal scoring',
        epilog="""
Examples:
    # Basic plasmid assignment
    accio find assembly.fasta -d /path/to/databases -o results
    
    # With short read data
    accio find assembly.fasta -r reads_R1.fastq reads_R2.fastq -d databases -o results
    
    # With long read data using minimap2
    accio find assembly.fasta -r nanopore_reads.fastq --minimap2 -d databases -o results
    
    # High sensitivity analysis
    accio find assembly.fasta -d databases -o results --min_identity 70 --min_coverage 70
        """
    )
    
    # Required arguments
    find_parser.add_argument(
        'assembly',
        help='Input genome assembly FASTA file'
    )
    
    find_parser.add_argument(
        '-d', '--database',
        required=True,
        help='Directory containing Accio plasmid databases'
    )
    
    find_parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory for results'
    )
    
    # Optional inputs
    find_parser.add_argument(
        '-r', '--reads',
        nargs='+',
        help='Read files for coverage analysis (FASTQ format)'
    )
    add_find_analyses_parsers(find_parser)

def add_find_analyses_parsers(parser):   
    # Analysis parameters
    analysis_group = parser.add_argument_group('Analysis Parameters')
    
    analysis_group.add_argument(
        '-t', '--threads',
        type=int,
        default=4,
        help='Number of threads to use (default: 4)'
    )
    
    analysis_group.add_argument(
        '--min_identity',
        type=float,
        default=75,
        help='Minimum BLAST identity threshold (default: 75)'
    )
    
    analysis_group.add_argument(
        '--min_coverage',
        type=float,
        default=75,
        help='Minimum BLAST coverage threshold (default: 75)'
    )
    
    analysis_group.add_argument(
        '--min_score',
        type=float,
        default=0.80,
        help='Minimum plasmid assignment score (default: 0.80)'
    )
    
    analysis_group.add_argument(
        '--circular_only',
        action='store_true',
        help='Only consider circular contigs as plasmid candidates'
    )
    analysis_group.add_argument('--use_pling_type_counts', action='store_true',
                                help='Use counts of matching pling types to assign plasmids')
    # Read mapping options
    mapping_group = parser.add_argument_group('Read Mapping Options')
    
    mapping_group.add_argument(
        '--minimap2',
        action='store_true',
        help='Use minimap2 instead of BWA for read mapping (recommended for long reads)'
    )
    
    mapping_group.add_argument(
        '--mapping_preset',
        choices=['sr', 'map-ont', 'map-pb', 'asm5', 'asm10', 'asm20'],
        help='Minimap2 preset (sr=short reads, map-ont=Oxford Nanopore, map-pb=PacBio)'
    )
    
    # Output options
    output_group = parser.add_argument_group('Output Options')
    
    output_group.add_argument(
        '-k', '--keep_intermediate',
        action='store_true',
        help='Keep intermediate analysis files'
    )
    
    output_group.add_argument(
        '--output_format',
        choices=['csv', 'tsv', 'json'],
        default='csv',
        help='Output format for results (default: csv)'
    )
    
    output_group.add_argument(
        '--sample_name',
        help='Sample name for output files (default: derived from assembly filename)'
    )


def add_create_parser(subparsers) -> None:
    """Add the 'create' subcommand parser."""
    create_parser = subparsers.add_parser(
        'create',
        help='Create plasmid databases from reference sequences',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='🏗️  Create Accio plasmid databases from reference plasmid sequences',
        epilog="""
Examples:
    # Create database from FASTA file
    accio create plasmids.fasta -o database_dir
    
    # Create database with metadata
    accio create plasmids.fasta -o database_dir --metadata plasmid_info.csv
    
    # Create database from directory of FASTA files
    accio create /path/to/plasmid_fastas/ -o database_dir --recursive
    
    # Update existing database
    accio create new_plasmids.fasta -o existing_db_dir --update
        """
    )
    
    # Input arguments
    create_parser.add_argument(
        'input',
        nargs='+',
        help='Input plasmid sequences (one or more FASTA files, a directory, or a file listing paths).'
    )
    
    create_parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory for database files'
    )
    
    # Database options
    db_group = create_parser.add_argument_group('Database Options')
    
    
    db_group.add_argument(
        '--update',
        action='store_true',
        help='Update existing database instead of creating new one'
    )
    
    db_group.add_argument(
        '--force',
        action='store_true',
        help='Overwrite existing database files'
    )
    
    # Processing options
    process_group = create_parser.add_argument_group('Processing Options')
    
    process_group.add_argument(
        '-t', '--threads',
        type=int,
        default=4,
        help='Number of threads for database building (default: 4)'
    )
    
    process_group.add_argument(
        '--min_length',
        type=int,
        default=10000,
        help='Minimum plasmid length to include (default: 1000 bp)'
    )
    
    process_group.add_argument(
        '--max_length',
        type=int,
        default=1000000,
        help='Maximum plasmid length to include (default: 1,000,000 bp)'
    )
    
    process_group.add_argument(
        '--cluster_threshold',
        type=float,
        default=0.95,
        help='Similarity threshold for clustering (default: 0.95)'
    )
    process_group.add_argument('--no_plasmid_check',
                                action='store_true', default=False,
                                help='Do not check for circularity/plasmid replicons')
    process_group.add_argument('--no_clustering',
                               action='store_true', default=False,
                               help='Do not cluster plasmid sequences')

    # Quality control
    qc_group = create_parser.add_argument_group('Quality Control')
    
    qc_group.add_argument(
        '--check_circular',
        action='store_true',
        help='Check for circular sequences and mark accordingly'
    )
    
    qc_group.add_argument(
        '--deduplicate',
        action='store_true',
        help='Remove duplicate sequences based on 100%% identity'
    )
    
    qc_group.add_argument(
        '--validate_sequences',
        action='store_true',
        help='Validate sequence quality and remove problematic sequences'
    )


def add_check_parser(subparsers) -> None:
    """Add the 'check' subcommand parser."""
    check_parser = subparsers.add_parser(
        'check',
        help='Check tool and database availability',
        description='🔧 Check system requirements and database integrity'
    )
    
    check_parser.add_argument(
        '--tools',
        action='store_true',
        help='Check availability of external tools'
    )
    
    check_parser.add_argument(
        '--databases',
        metavar='DIR',
        help='Check database files in specified directory'
    )
    
    check_parser.add_argument(
        '--all',
        action='store_true',
        help='Check both tools and databases (requires --databases)'
    )

def add_collate_parser(subparsers) -> None:
    """Add the 'combine' subcommand parser"""
    collate_parser = subparsers.add_parser('collate', help='Combine results from multiple genomes')
    collate_parser.add_argument('input_pattern', help='Glob pattern for input files (e.g., "results/*/*_plasmids_chosen.csv")')
    collate_parser.add_argument('-o', '--output', required=True, help='Output file for combined results (e.g., summary.csv)')

def add_query_parser(subparsers) -> None:
    """Add the 'query' subcommand parser"""
    query_parser = subparsers.add_parser('query', help='Get info about a database plasmid')
    query_parser.add_argument('database', help='Path to the database directory')
    query_parser.add_argument('--plasmid_id', help='ID of the plasmid to query')
    query_parser.add_argument('--pling_type', help='Pling type to query')
    query_parser.add_argument('--community', help='Plasmid community to query')
    query_parser.add_argument('--replicon_type', help='Replicon type to query')

def add_multi_parser(subparsers) -> None:
    """Add the 'multi' subcommand parser"""
    multi_parser = subparsers.add_parser('multi', help='Find and assign plasmids in many genomes')
    multi_parser.add_argument('--assemblies', nargs='+', required=True, help="Directory containing FASTA files.")
    multi_parser.add_argument('--reads_dir', help="Directory containing read files.")
    multi_parser.add_argument('-d', '--database', required=True, help="Path to the database directory.")
    multi_parser.add_argument('-o', '--output_dir', required=True, help="Parent directory for all output.")
    multi_parser.add_argument('--collate', action='store_true', help="Automatically collate results into a summary file.")

    add_find_analyses_parsers(multi_parser) 


def handle_find_command(args: argparse.Namespace) -> int:
    """Handle the 'find' subcommand."""
    logger = logging.getLogger(__name__)
    
    try:
        # Validate inputs
        if not os.path.exists(args.assembly):
            logger.error(f"Assembly file does not exist: {args.assembly}")
            return 1
            
        if args.reads:
            for read_file in args.reads:
                if not os.path.exists(read_file):
                    logger.error(f"Read file does not exist: {read_file}")
                    return 1
                    
        if not os.path.exists(args.database):
            logger.error(f"Database directory does not exist: {args.database}")
            return 1
        
        # Check prerequisites
        logger.info("Checking system requirements...")
        tool_availability = check_tool_availability()
        essential_tools = ['blastn', 'mash', 'nucmer', 'show-coords', 'samtools']
        
        if args.minimap2:
            essential_tools.append('minimap2')
        # else:
        #     essential_tools.extend(['bwa', 'bwa-mem2'])
            
        missing_tools = [tool for tool in essential_tools 
                        if not tool_availability.get(tool, False)]
        
        if missing_tools:
            logger.error(f"Missing essential tools: {', '.join(missing_tools)}")
            return 1
            
        # Validate databases
        # db_status = validate_databases(args.database)
        # missing_dbs = [db for db, exists in db_status.items() if not exists]
        
        # if missing_dbs:
        #     logger.error(f"Missing database files: {', '.join(missing_dbs)}")
        #     return 1
        
        # Create output directory
        os.makedirs(args.output, exist_ok=True)
        
        # Create configuration
        config = get_config()
        config.MIN_IDENTITY = args.min_identity
        config.MIN_COVERAGE = args.min_coverage
        config.MIN_PLASMID_SCORE = args.min_score
        config.THREADS = args.threads
        config.USE_PLING_COMMUNITY_COUNTS = args.use_pling_type_counts
        
        if args.circular_only:
            config.CIRCULAR_ONLY = True
            
        validate_config(config)
        
        # Determine sample name
        sample_name = args.sample_name
        if not sample_name:
            sample_name = Path(args.assembly).stem
        
        # Initialize workflow
        logger.info(f"🔍 Starting plasmid analysis for sample: {sample_name}")
        workflow = PlasmidAnalysisWorkflow(config)
        
        # Run analysis
        results = workflow.run(
            fasta_file=args.assembly,
            reads=args.reads,
            db_dir=args.database,
            output_dir=args.output,
            use_minimap=args.minimap2,
            keep_intermediate=args.keep_intermediate,
            sample_name=sample_name,
            output_format=args.output_format
        )
        
        # Report results
        logger.info("✅ Analysis completed successfully!")
        
        if results.get('plasmids_assigned'):
            num_plasmids = len(results['plasmids_assigned'])
            logger.info(f"🧬 Assigned {num_plasmids} plasmid(s)")
            
            # Show assigned plasmids
            for plasmid in results['plasmids_assigned'][:5]:  # Show first 5
                logger.info(f"   └─ {plasmid.id} (score: {plasmid.score:.3f})")
            
            if num_plasmids > 5:
                logger.info(f"   └─ ... and {num_plasmids - 5} more")
                
        if results.get('novel_plasmid_contigs'):
            num_novel = len(results['novel_plasmid_contigs'])
            logger.info(f"🆕 Found {num_novel} novel plasmid contig(s)")
            
        logger.info(f"📁 Results written to: {args.output}")
        
        return 0
        
    except Exception as e:
        logger.error(f"❌ Analysis failed: {str(e)}")
        logger.debug("Exception details:", exc_info=True)
        return 1


def handle_create_command(args: argparse.Namespace) -> int:
    """Handle the 'create' subcommand."""
    logger = logging.getLogger(__name__)
    
    try:
        # Validate inputs
        for path in args.input:
            if not os.path.exists(path):
                # An exception for a file of file paths, which is handled later
                if len(args.input) == 1 and not any(path.endswith(ext) for ext in ['.fasta', '.fa', '.fna']):
                    continue
                logger.error(f"Input path does not exist: {path}")
                return 1
            
        # Check if output directory exists
        if os.path.exists(args.output) and not args.update and not args.force:
            logger.error(f"Output directory exists: {args.output}")
            logger.error("Use --update to add to existing database or --force to overwrite")
            return 1
            
        # Create output directory
        os.makedirs(args.output, exist_ok=True)
        
        # Check for required tools
        tool_availability = check_tool_availability()
        required_tools = ['makeblastdb', 'mash']
        missing_tools = [tool for tool in required_tools 
                        if not tool_availability.get(tool, False)]
        
        if missing_tools:
            logger.error(f"Missing required tools: {', '.join(missing_tools)}")
            return 1
        
        # Initialize database builder
        logger.info("🏗️  Starting database creation...")
        builder = DatabaseBuilder(
            threads=args.threads,
            min_length=args.min_length,
            max_length=args.max_length
        )
        
        # Collect input files
        input_files = set()
        fasta_extensions = ['.fasta', '.fa', '.fna']

        for item in args.input:
            path = Path(item)
            if path.is_dir():
                logger.info(f"Searching for FASTA files in directory: {path}")
                glob_method = path.rglob if args.recursive else path.glob
                for ext in fasta_extensions:
                    input_files.update(glob_method(f"*{ext}"))
            elif path.is_file():
                # If it's a FASTA file, add it directly
                if path.suffix.lower() in fasta_extensions:
                    input_files.add(path)
                # Otherwise, assume it's a file containing a list of paths
                else:
                    logger.info(f"Reading file paths from: {path}")
                    try:
                        with open(path, 'r') as f:
                            for line in f:
                                file_path = Path(line.strip())
                                if file_path.exists() and file_path.is_file():
                                    input_files.add(file_path)
                                else:
                                    logger.warning(f"Path from list file does not exist or is not a file: {file_path}")
                    except Exception as e:
                        logger.error(f"Could not read file list from {path}: {e}")
                        return 1

        if not input_files:
            logger.error("No FASTA files found in input")
            return 1
            
        logger.info(f"📄 Found {len(input_files)} unique FASTA file(s) to process.")
        
        # Build database
        results = builder.build_database(
            input_files=[str(p) for p in input_files],
            output_dir=args.output,
            update_existing=args.update,
            check_circular=args.check_circular,
            no_cluster=args.no_clustering,
            deduplicate=args.deduplicate,
            use_all=args.no_plasmid_check,
            validate_sequences=args.validate_sequences,
            cluster_threshold=args.cluster_threshold
        )
        
        # Report results
        logger.info("✅ Database creation completed!")
        logger.info(f"📊 Processed {results['total_sequences']} sequences")
        logger.info(f"🧬 Created database with {results['final_sequences']} plasmids")
        
        if results.get('clusters_created'):
            logger.info(f"🔗 Organized into {results['clusters_created']} clusters")
            
        if results.get('duplicates_removed'):
            logger.info(f"🗑️  Removed {results['duplicates_removed']} duplicate sequences")
            
        logger.info(f"📁 Database files written to: {args.output}")
        
        return 0
        
    except Exception as e:
        logger.error(f"❌ Database creation failed: {str(e)}")
        logger.debug("Exception details:", exc_info=True)
        return 1


def handle_check_command(args: argparse.Namespace) -> int:
    """Handle the 'check' subcommand."""
    exit_code = 0
    
    if args.tools or args.all:
        print("🔧 Checking external tool availability...")
        tool_availability = check_tool_availability()
        
        # Group tools by category
        essential_tools = ['blastn', 'makeblastdb', 'mash', 'nucmer', 'show-coords', 'samtools']
        mapping_tools = ['bwa', 'minimap2'] # bwa-mem2 is often separate
        optional_tools = ['PLASMe.py', 'mob_recon', 'bwa-mem2']
        
        def check_tool_group(tools: List[str], group_name: str) -> bool:
            print(f"\n  {group_name}:")
            all_available = True
            for tool in tools:
                info = tool_availability.get(tool)
                if info:
                    status = "✅" if info['available'] else "❌"
                    print(f"    {status} {tool}")
                    if not info['available']:
                        print(f"       └─ Install with: {info['install_hint']}")
                        if tool in essential_tools:
                            all_available = False
            return all_available
        
        essential_ok = check_tool_group(essential_tools, "Essential tools")
        check_tool_group(mapping_tools, "Read mapping tools")
        check_tool_group(optional_tools, "Optional tools")
        
        if not essential_ok:
            print(f"\n❌ Some essential tools are missing!")
            exit_code = 1
        else:
            print(f"\n✅ All essential tools are available!")
    
    if args.databases or args.all:
        if not args.databases:
            print("\n❌ Database directory not specified (use --databases DIR)")
            return 1
            
        print(f"\n🗄️  Checking database files in {args.databases}...")
        
        if not os.path.exists(args.databases):
            print(f"❌ Database directory does not exist: {args.databases}")
            return 1
            
        db_status = validate_databases(args.databases)
        
        all_available = True
        for db, exists in db_status.items():
            status = "✅" if exists else "❌"
            print(f"  {status} {db}")
            if not exists:
                all_available = False
                
        if not all_available:
            print(f"\n❌ Some database files are missing!")
            exit_code = 1
        else:
            print(f"\n✅ All database files are available!")
    
    if not args.tools and not args.databases and not args.all:
        print("❓ No check specified. Use --tools, --databases, or --all")
        return 1
        
    return exit_code

def handle_query_command(args: argparse.Namespace) -> int:
    pass

def handle_multi_command(args: argparse.Namespace) -> int:
    """Handle the 'multi' subcommand to run analysis on multiple samples."""
    logger = logging.getLogger(__name__)
    
    # Find all assembly files
    assembly_files = args.assemblies
    if not assembly_files:
        logger.error(f"No assembly files found in {args.assemblies}")
        return 1

    logger.info(f"Found {len(assembly_files)} assemblies to process.")
    logger.info(f"Running up to {args.threads} samples in parallel.")
    
    # Create the main output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    jobs_to_run = []
    for assembly_file in assembly_files:
        sample_name = Path(assembly_file).stem
        logger.info(f"--- Processing sample: {sample_name} ---")
        
        sample_output_dir = Path(args.output_dir) / sample_name
        sample_output_dir.mkdir(exist_ok=True)

        # Prepare arguments for the 'find' command for this specific sample
        job_args = argparse.Namespace(**vars(args))
        job_args.assembly = assembly_file
        job_args.output = str(sample_output_dir)
        job_args.sample_name = sample_name

        # Find corresponding read files if a reads directory is provided
        if args.reads_dir:
            # Assumes reads are named like: {sample_name}_R1.fastq.gz, {sample_name}_R2.fastq.gz, etc.
            read_files = sorted(glob.glob(f"{args.reads_dir}/{sample_name}*"))
            if read_files:
                logger.info(f"Found {len(read_files)} read file(s) for {sample_name}")
                job_args.reads = read_files
            else:
                logger.warning(f"No read files found for {sample_name} in {args.reads_dir}")
                job_args.reads = None
        else:
            job_args.reads = None

        jobs_to_run.append(job_args)
    num_parallel = int(args.threads / 8) if args.threads >=8 else 1
    # Run jobs in parallel
    with multiprocessing.Pool(processes=num_parallel) as pool:
        pool.map(handle_find_command, jobs_to_run)
    
    logger.info("--- Multi-analysis complete! ---")    
    # Automatically collate results if requested
    if args.collate:
        logger.info("Collating results...")
        collate_pattern = f"{args.output_dir}/*/{Path(FilePatterns.PLASMID_DATA_OUTPUT).name.format(sample='*')}"
        collate_output = Path(args.output_dir) / "summary_report.csv"
        
        collate_args = argparse.Namespace(
            input_pattern=collate_pattern,
            output=str(collate_output)
        )
        handle_combine_command(collate_args)
    else:
        logger.info(f"To combine all results into a single summary, run: accio collate '{args.output_dir}/*/{Path(FilePatterns.PLASMID_DATA_OUTPUT).name.format(sample='*')}' -o {args.output_dir}/summary_report.csv")
    return 0

def handle_combine_command(args: argparse.Namespace) -> int:
    """Collates results from an 'accio multi' run."""
    logger = logging.getLogger(__name__)
    results_files = glob.glob(args.input_pattern, recursive=True)
    
    all_dfs = []
    for f in results_files:
        try:
            df = pd.read_csv(f)
            sample_name = f.split(os.sep)[-1].split('_')[0]
            df['sample'] = sample_name

            all_dfs.append(df)
        except pd.errors.EmptyDataError:
            logger.warning(f"Skipping empty file {f}")
            continue
            
    if not all_dfs:
        logger.warning("No result files found to collate.")
        return 1

    summary_df = pd.concat(all_dfs, ignore_index=True)
    summary_df.to_csv(args.output, index=False)
    logger.info(f"Summary report written to {args.output}")
    return 0


def main() -> int:
    """Main entry point for Accio CLI."""
    parser = create_main_parser()
    args = parser.parse_args()
    
    # Handle case where no command is provided
    if not args.command:
        parser.print_help()
        return 1
    
    # Set up logging
    setup_logging(args.log_level, args.log_file)
    logger = logging.getLogger(__name__)
    
    # Print banner
    print("🪄 Accio - Plasmid Analysis Tool")
    print("   Summoning plasmids from bacterial genomes...")
    print()
    
    try:
        # Route to appropriate command handler
        if args.command == 'find':
            return handle_find_command(args)
        elif args.command == 'create':
            return handle_create_command(args)
        elif args.command == 'check':
            return handle_check_command(args)
        elif args.command == 'collate':
            return handle_combine_command(args)
        elif args.command == 'query':
            return handle_query_command(args)
        elif args.command == 'multi':
            return handle_multi_command(args)
        else:
            logger.error(f"Unknown command: {args.command}")
            return 1
            
    except KeyboardInterrupt:
        logger.info("\n⏸️  Analysis interrupted by user")
        return 130
    except Exception as e:
        logger.error(f"❌ Unexpected error: {str(e)}")
        logger.debug("Exception details:", exc_info=True)
        return 1


if __name__ == '__main__':
    sys.exit(main())