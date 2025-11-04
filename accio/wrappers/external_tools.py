"""External tool wrappers for plasmid analysis."""

import subprocess
import shlex
import os
import io
import pandas as pd
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Any
from Bio import SeqIO

from ..config import (
    DEFAULT_OUTFMT, DEFAULT_FIELDS, EXTRA_FIELDS, 
    NUCMER_COORD_COLS, NUCMER_SNP_COLS, AnalysisConfig
)


class ExternalToolError(Exception):
    """Exception raised when external tools fail."""
    pass


class BLASTRunner:
    """Wrapper for BLAST operations."""
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        self.config = config or AnalysisConfig()
        
    def run_blastn(self, query: str, database: str, 
                   output_file: Optional[str] = None,
                   additional_args: str = '',
                   sample_name: Optional[str] = None) -> pd.DataFrame:
        """
        Run blastn against a database.
        
        Args:
            query: Path to query FASTA file
            database: Path to BLAST database
            output_file: Optional output file path
            additional_args: Additional BLAST arguments
            sample_name: Sample name for tracking
            
        Returns:
            DataFrame containing BLAST results
            
        Raises:
            ExternalToolError: If BLAST execution fails
        """
        # Prepare output format
        if additional_args and '-outfmt' in additional_args:
            output_fmt = ''
            fields = []
        else:
            output_fmt = f"-outfmt {DEFAULT_OUTFMT}"
            fields = DEFAULT_FIELDS + EXTRA_FIELDS
            
        # Prepare output redirection
        output_redirect = f'-out {output_file}' if output_file else ''
        
        # Build command
        cmd = f'blastn -query {query} -db {database} {output_redirect} {output_fmt} {additional_args}'
        
        try:
            if not output_file:
                # Capture output directly
                result = subprocess.check_output(shlex.split(cmd), stderr=subprocess.PIPE)
                result_io = io.StringIO(result.decode())
                df = pd.read_csv(result_io, sep='\t', names=fields, dtype={'qseqid': str}, header=None)
            else:
                # Run with file output
                subprocess.check_call(shlex.split(cmd))
                df = pd.read_csv(output_file, sep='\t', names=fields, dtype={'qseqid': str}, header=None)
                
            if sample_name:
                df['Sample'] = sample_name
                
            return df
            
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.decode() if e.stderr else str(e)
            raise ExternalToolError(f"BLAST failed: {error_msg}")
        except Exception as e:
            raise ExternalToolError(f"BLAST execution error: {str(e)}")
    
    def run_makeblastdb(self, input_fasta: str, output_file: str,
                        dbtype: Optional[str] = 'nucl',
                        out: Optional[str] = None,
                        parse_seqids: Optional[bool] = True,
                        title: Optional[str] = None):
        """
        Run makeblastdb to create a BLAST database.
        """
        cmd = ['makeblastdb', '-in', input_fasta, 
                '-out', output_file, '-dbtype', dbtype]
        if out:
            cmd += ['-out', out]
        if parse_seqids:
            cmd += ['-parse_seqids']
        if title:
            cmd += ['-title', title]
        try:
            subprocess.run(cmd)
        except:
            raise ExternalToolError(f"makeblastdb failed: {str(e)}")


class MashRunner:
    """Wrapper for Mash operations."""
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        self.config = config or AnalysisConfig()

    def run_sketch(self, input_fasta: str, output_file: str,
                   k: Optional[int] = 21, s: Optional[int] = 1000,
                   reads: Optional[bool] = False,
                   multi: Optional[bool] = False):
        """
        Run mash sketch.
        
        Args:
            input_fasta: Path to input FASTA file
            output_file: Output file path
            
        Raises:
            ExternalToolError: If Mash execution fails
        """
        cmd = ['mash', 'sketch',
            '-o', output_file,
            '-k', str(k),
            '-s', str(s),
            '-p', str(self.config.THREADS),
            input_fasta
        ]
        if reads:
            cmd.append('-r')
        if multi:
            cmd.append('-i')
        try:
            subprocess.run(cmd, check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(f"Mash sketch failed: {str(e)}")
        
    def run_screen(self, database: str, query: str, 
                   output_file: Optional[str] = None) -> pd.DataFrame:
        """
        Run mash screen against a database.
        
        Args:
            database: Path to Mash database
            query: Path to query sequences
            output_file: Optional output file path
            
        Returns:
            DataFrame containing Mash results
            
        Raises:
            ExternalToolError: If Mash execution fails
        """
        cmd = ['mash', 'screen', database, query]
        
        try:
            result = subprocess.run(cmd, capture_output=True, check=True)
            result_io = io.StringIO(result.stdout.decode())
            
            df = pd.read_csv(result_io, sep='\t', header=None, names=[
                'pident', 'shared-hashes', 'median-mult', 'pvalue', 'sseqid', 'stitle'
            ])
            
            if output_file:
                df.to_csv(output_file, index=False)
                
            return df
            
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.decode() if e.stderr else str(e)
            raise ExternalToolError(f"Mash screen failed: {error_msg}")
        
class NucmerRunner:
    """Wrapper for Nucmer operations."""
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        self.config = config or AnalysisConfig()
        
    def run_nucmer(self, reference: str, query: str, 
                   prefix: str) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, str]]:
        """
        Run nucmer alignment and process results.
        
        Args:
            reference: Path to reference sequences
            query: Path to query sequences  
            prefix: Output file prefix
            
        Returns:
            Tuple of (coordinates_df, snps_df, output_files)
            
        Raises:
            ExternalToolError: If Nucmer execution fails
        """
        try:
            # Run nucmer
            nucmer_cmd = (
                f'nucmer --diagdiff {self.config.NUCMER_DIAGDIFF} '
                f'--breaklen {self.config.NUCMER_BREAKLEN} '
                f'--maxmatch -p {prefix} {reference} {query}'
            )
            subprocess.check_call(shlex.split(nucmer_cmd))
            delta_file = f'{prefix}.delta'
            
            # Filter delta file
            filtered_delta_file = f'{prefix}.1delta'
            delta_filter_cmd = f'delta-filter -qr {delta_file}'
            with open(filtered_delta_file, 'w') as outfile:
                subprocess.check_call(shlex.split(delta_filter_cmd), stdout=outfile)
                
            # Extract coordinates
            coords_file = f'{prefix}.coords'
            coords_cmd = f'show-coords -TrcldH -I 80 {delta_file}'
            with open(coords_file, 'w') as outfile:
                subprocess.check_call(shlex.split(coords_cmd), stdout=outfile)
                
            # Extract SNPs
            snps_file = f'{prefix}.snps'
            snps_cmd = f'show-snps -TrH {delta_file}'
            with open(snps_file, 'w') as outfile:
                subprocess.check_call(shlex.split(snps_cmd), stdout=outfile)
                
            # Read results
            coords_df = pd.read_csv(
                coords_file, sep='\t', header=None, names=NUCMER_COORD_COLS
            )
            snps_df = pd.read_csv(
                snps_file, sep='\t', header=None, names=NUCMER_SNP_COLS
            )
            
            # Ensure proper data types
            coords_df['qseqid'] = coords_df['qseqid'].astype(str)
            snps_df['query'] = snps_df['query'].astype(str)

            output_files = {
                "delta": delta_file,
                "filtered_delta": filtered_delta_file,
                "coords": coords_file,
                "snps": snps_file,
            }
            return coords_df, snps_df, output_files
        
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(f"Nucmer execution failed: {str(e)}")
        except Exception as e:
            raise ExternalToolError(f"Nucmer processing error: {str(e)}")
    


class PLASMeRunner:
    """Wrapper for PLASMe operations."""
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        self.config = config or AnalysisConfig()
        
    def run_plasme(self, input_fasta: str, output_dir: str) -> Tuple[pd.DataFrame, Dict[str, str]]:
        """
        Run PLASMe plasmid prediction.
        
        Args:
            input_fasta: Path to input FASTA file
            output_dir: Output directory
            
        Returns:
            Tuple of (DataFrame with PLASMe results, dictionary of output files)
            
        Raises:
            ExternalToolError: If PLASMe execution fails
        """
        output_fasta = os.path.join(output_dir, 'plasme.fasta')
        cmd = ['PLASMe.py', '--temp', output_dir, input_fasta, output_fasta]
        
        try:
            subprocess.check_call(cmd)
            
            report_file = f'{output_fasta}_report.csv'
            if os.path.exists(report_file):
                df = pd.read_csv(report_file, sep='\t')
                output_files = {"report": report_file, "fasta": output_fasta}
                return df, output_files
            else:
                return pd.DataFrame(), {}
                
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(f"PLASMe execution failed: {str(e)}")
        except Exception as e:
            raise ExternalToolError(f"PLASMe processing error: {str(e)}")




class ReadMappingRunner:
    """Wrapper for read mapping operations."""
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        self.config = config or AnalysisConfig()
        
    def map_reads(self, reference: str, reads: List[str], 
                  output_prefix: str,  use_minimap: bool = False) -> str:
        """
        Map reads to reference and create sorted BAM file.
        
        Args:
            reference: Path to reference sequences
            reads: List of read file paths
            output_prefix: Output file prefix
            use_minimap: Whether to use minimap2 instead of BWA
            
        Returns:
            Path to output BAM file
            
        Raises:
            ExternalToolError: If mapping fails
        """
        bam_file = f'{output_prefix}.bam'
        if os.path.isfile(bam_file) and not self.config.OVERWRITE_FILES:
            return bam_file
        try:
            if not use_minimap:
                # Index reference with BWA
                if not os.path.isfile(f'{reference}.bwt'):
                    index_cmd = ['bwa', 'index', reference]
                    subprocess.check_call(index_cmd)
                
                # # Map with BWA MEM
                map_cmd = ['bwa', 'mem', '-t', str(self.config.THREADS), reference] + reads
            else:
                # Map with minimap2
                map_cmd = [
                    'minimap2', '-a', '-P', '-x', 'map-ont', 
                    '-t', str(self.config.THREADS), reference
                ] + reads
                
            # Sort with samtools
            sort_cmd = ['samtools', 'sort', '-o', bam_file, '-']
            
            # Create pipeline
            map_proc = subprocess.Popen(map_cmd, stdout=subprocess.PIPE)
            sort_proc = subprocess.Popen(sort_cmd, stdin=map_proc.stdout)
            
            map_proc.stdout.close()
            sort_proc.wait()
            
            if sort_proc.returncode != 0:
                raise subprocess.CalledProcessError(sort_proc.returncode, 'samtools sort')
                
            # Index BAM file
            index_cmd = ['samtools', 'index', bam_file]
            subprocess.check_call(index_cmd)
            
            return bam_file
            
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(f"Read mapping failed: {str(e)}")
    
    def calculate_coverage(self, bam_file: str, output_file: str) -> Tuple[pd.DataFrame, str]:
        """
        Calculate coverage statistics from BAM file.
        
        Args:
            bam_file: Path to BAM file
            output_file: Output file path
            
        Returns:
            Tuple of (DataFrame with coverage statistics, path to output file)
            
        Raises:
            ExternalToolError: If coverage calculation fails
        """
        try:
            cmd = ['samtools', 'coverage', bam_file]
            with open(output_file, 'w') as outfile:
                subprocess.check_call(cmd, stdout=outfile)
                
            cov_df = pd.read_csv(output_file, sep='\t')
            cov_df['contig'] = cov_df['#rname'].astype(str)
            cov_df['length'] = cov_df['endpos'] - cov_df['startpos']
            
            return cov_df, output_file
            
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(f"Coverage calculation failed: {str(e)}")
        

class MOBRunner:
    """Wrapper for MOB-suite operations."""
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        self.config = config or AnalysisConfig()
        
    def run_mob_recon(self, input_fasta: str, output_dir: str, plasmid_db: str = None) -> Dict[str, str]:
        """
        Run MOB-recon for plasmid reconstruction.
        
        Args:
            input_fasta: Path to input FASTA file
            output_dir: Output directory
            plasmid_db: Path to plasmid database (optional)
            
        Raises:
            ExternalToolError: If MOB-recon execution fails. Returns dict of output files.
        """
        cmd = [
            'mob_recon', '-i', input_fasta, 
            '-o', output_dir, '--force', '-u'
        ]
     
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(f"MOB-recon failed: {str(e)}")
        
        output_files = {
            "mobtyper_results": os.path.join(output_dir, 'mobtyper_results.txt'),
            "contig_report": os.path.join(output_dir, 'contig_report.txt'),
        }
        return output_files
    
    def run_mob_typer(self, input_fasta: str, output_file: str, 
                      multi: Optional[bool] = False,
                      sample_id: Optional[str] = None,
                      primary_cluster_dist: Optional[float] = 0.06,
                      secondary_cluster_dist: Optional[float] = 0.025,
                      plasmid_mash_db: Optional[str] = None,
                      plasmid_meta: Optional[str] = None) -> str:
        """
        Run MOB-typer for plasmid type identification.
        """
        cmd = ['mob_typer', 
               '-i', input_fasta, 
               '-o', output_file,
               '-n', str(self.config.THREADS)]
        if multi:
            cmd += ['-x']
        if sample_id:
            cmd += ['-s', sample_id]
        if plasmid_meta:
            cmd +=['-m', plasmid_meta,
                   '--plasmid_mash_db', plasmid_mash_db]
        
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(f"MOB-typer failed: {str(e)}")
        return output_file

    def run_mob_cluster(self, input_fasta: str, output_dir: str,
                         mob_typer_results: str,
                         taxonomy_file: str,
                         previous_mob_db: Optional[str] = None, previous_mob_fastas: Optional[str] = None,
                        primary_cluster_dist: Optional[float] = 0.06, secondary_cluster_dist: Optional[float] = 0.025):
        """
        Run MOB-cluster for plasmid clustering.
        """
        cmd = ['mob_cluster', 
                   '-m' 'build' if  not previous_mob_db else 'update',
                   '-t', taxonomy_file, 
                   '-p', mob_typer_results, 
                   '-o', output_dir, 
                   '-f', input_fasta,
                   '--primary_cluster_dist', str(primary_cluster_dist),
                   '--secondary_cluster_dist', str(secondary_cluster_dist),
                   ]
        if previous_mob_db:
            cmd += ['-c', previous_mob_db,
                    '-r', previous_mob_fastas]
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(f"MOB-cluster failed: {str(e)}")



class SkaniRunner:
    """Wrapper for skani operations."""
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        self.config = config or AnalysisConfig()
    
    def run_ska_dist(self, plasmids, output_file):
        cmd = [
            'skani', 'dist', '--ri', plasmids, '--qi', plasmids,
            '--medium', '-m', '200', '-o', str(output_file)
        ]
        
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(f"MOB-typer failed: {str(e)}")
                
class AMRRunner:
    """Wrapper for AMRFinder operations."""
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        self.config = config or AnalysisConfig()
        
    def run_amrfinder(self, input_fasta: str, output_file: str) -> str:
        """
        Run AMRFinder for antimicrobial resistance gene identification.
        
        Args:
            input_fasta: Path to input FASTA file
            output_file: Output file path
            
        Returns:
            Path to output file
            
        Raises:
            ExternalToolError: If AMRFinder execution fails
        """
        cmd = [
            'amrfinder', '-n', input_fasta,
            '-o', output_file,
            '--threads', str(self.config.THREADS)
        ]
        
        try:
            subprocess.check_call(cmd)
            return output_file
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(f"AMRFinder failed: {str(e)}")

def check_tool_availability() -> Dict[str, bool]:
    """
    Check availability of external tools.
    
    Returns:
        Dictionary mapping tool names to availability status
    """
    tools = {
        'blastn': 'blastn -version',
        'mash': 'mash --version',
        'nucmer': 'nucmer --version',
        'show-coords': 'show-coords -h',
        'show-snps': 'show-snps -h',
        'delta-filter': 'delta-filter -h',
        'bwa': 'bwa',
        'bwa-mem2': 'bwa',
        'minimap2': 'minimap2 --version',
        'samtools': 'samtools --version',
        'PLASMe.py': 'PLASMe.py --help',
        'mob_recon': 'mob_recon --version',
        'makeblastdb': 'makeblastdb -h',
        
    }
    
    availability = {}
    
    for tool, cmd in tools.items():
        try:
            subprocess.check_output(
                shlex.split(cmd), 
                stderr=subprocess.DEVNULL,
                timeout=10
            )
            availability[tool] = True
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
            availability[tool] = False
            
    return availability


def validate_databases(db_dir: str) -> Dict[str, bool]:
    """
    Validate that required database files exist.
    
    Args:
        db_dir: Database directory path
        
    Returns:
        Dictionary mapping database names to existence status
    """
    from ..config import DatabaseFiles
    
    required_files = {
        'plasmid_db': os.path.join(db_dir, DatabaseFiles.PLASMID_DB),
        'plasmidfinder_db': os.path.join(db_dir, DatabaseFiles.PLASMIDFINDER_DB),
        'resfinder_db': os.path.join(db_dir, DatabaseFiles.RESFINDER_DB),
        'mob_db': os.path.join(db_dir, DatabaseFiles.MOB_DB),
        'rep_db': os.path.join(db_dir, DatabaseFiles.REP_DB),
        'plasmid_info': os.path.join(db_dir, DatabaseFiles.PLASMID_INFO),
        'mash_db': os.path.join(db_dir, DatabaseFiles.MASH_DB)
    }
    
    existence = {}
    for name, path in required_files.items():
        existence[name] = os.path.exists(path)
        
    return existence