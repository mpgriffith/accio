"""External tool wrappers for plasmid analysis."""

import subprocess
import shlex
import shutil
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
        # Build command list safely for paths with spaces
        cmd = ['blastn', '-query', str(query), '-db', str(database), '-num_threads', str(self.config.THREADS)]
        if output_file:
            cmd.extend(['-out', str(output_file)])
            
        if additional_args and '-outfmt' in additional_args:
            fields = []
        else:
            cmd.extend(['-outfmt', DEFAULT_OUTFMT.strip('"').strip("'")])
            fields = DEFAULT_FIELDS + EXTRA_FIELDS
            
        if additional_args:
            cmd.extend(shlex.split(additional_args))
        
        try:
            if not output_file:
                # Capture output directly
                result = subprocess.check_output(cmd, stderr=subprocess.PIPE)
                result_io = io.StringIO(result.decode())
                df = pd.read_csv(result_io, sep='\t', names=fields, dtype={'qseqid': str}, header=None)
            else:
                # Run with file output
                subprocess.check_call(cmd)
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
            subprocess.run(cmd, check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.decode() if e.stderr else str(e)
            raise ExternalToolError(f"makeblastdb failed: {error_msg}")
        except Exception as e:
            raise ExternalToolError(f"makeblastdb execution error: {str(e)}")


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
        cmd = ['mash', 'screen', '-p', str(self.config.THREADS), str(database), str(query)]
        
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
            nucmer_cmd = [
                'nucmer',
                '--diagdiff', str(self.config.NUCMER_DIAGDIFF),
                '--breaklen', str(self.config.NUCMER_BREAKLEN),
                '--maxmatch',
                '-p', str(prefix),
                str(reference),
                str(query)
            ]
            subprocess.check_call(nucmer_cmd)
            delta_file = f'{prefix}.delta'
            
            # Filter delta file
            filtered_delta_file = f'{prefix}.1delta'
            delta_filter_cmd = ['delta-filter', '-qr', delta_file]
            with open(filtered_delta_file, 'w') as outfile:
                subprocess.check_call(delta_filter_cmd, stdout=outfile)
                
            # Extract coordinates
            coords_file = f'{prefix}.coords'
            coords_cmd = ['show-coords', '-TrcldH', '-I', '80', delta_file]
            with open(coords_file, 'w') as outfile:
                subprocess.check_call(coords_cmd, stdout=outfile)
                
            # Extract SNPs
            snps_file = f'{prefix}.snps'
            snps_cmd = ['show-snps', '-TrH', delta_file]
            with open(snps_file, 'w') as outfile:
                subprocess.check_call(snps_cmd, stdout=outfile)
                
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
                preset = getattr(self.config, 'MAPPING_PRESET', None) or 'map-ont'
                map_cmd = [
                    'minimap2', '-a', '-P', '-x', preset, 
                    '-t', str(self.config.THREADS), reference
                ] + reads
                
            # Sort with samtools
            sort_cmd = ['samtools', 'sort', '-o', bam_file, '-']
            
            # Create pipeline
            map_proc = subprocess.Popen(map_cmd, stdout=subprocess.PIPE)
            sort_proc = subprocess.Popen(sort_cmd, stdin=map_proc.stdout)
            
            map_proc.stdout.close()
            sort_proc.wait()
            map_proc.wait()
            
            if map_proc.returncode != 0:
                raise subprocess.CalledProcessError(map_proc.returncode, map_cmd[0])
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
    """Wrapper for MOB-suite operations with multi-environment support."""
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        self.config = config or AnalysisConfig()
        
    def _resolve_tool(self, tool_name: str) -> List[str]:
        """
        Resolve the command prefix for a MOB-suite tool.
        Checks in order:
        1. Custom directory configured in config.MOB_BIN_DIR or ACCIO_MOB_BIN_DIR
        2. Direct executable on PATH (shutil.which)
        3. Configured Conda environment (config.MOB_CONDA_ENV or ACCIO_MOB_ENV, default 'mob_suite')
           via conda/mamba/micromamba run
        4. Fallback to direct tool name
        """
        bin_dir = getattr(self.config, 'MOB_BIN_DIR', None) or os.environ.get('ACCIO_MOB_BIN_DIR')
        if bin_dir:
            cand = os.path.join(bin_dir, tool_name)
            if shutil.which(cand):
                return [cand]
            if os.name == 'nt' and shutil.which(cand + '.exe'):
                return [cand + '.exe']

        if shutil.which(tool_name):
            return [tool_name]

        env_name = getattr(self.config, 'MOB_CONDA_ENV', None) or os.environ.get('ACCIO_MOB_ENV', 'mob_suite')
        if env_name:
            conda_exe = shutil.which('mamba') or shutil.which('conda') or shutil.which('micromamba')
            if conda_exe:
                return [conda_exe, 'run', '--no-capture-output', '-n', env_name, tool_name]

        return [tool_name]

    def _format_mob_error(self, tool_name: str, err: Exception) -> str:
        """Provide detailed diagnostic error message for MOB-suite failures."""
        env_name = getattr(self.config, 'MOB_CONDA_ENV', None) or os.environ.get('ACCIO_MOB_ENV', 'mob_suite')
        return (
            f"{tool_name} execution failed ({str(err)}).\n"
            f"Note: MOB-suite requires pandas <= 1.5.3 while Pling/Plasnet requires pandas >= 2.0.\n"
            f"If MOB-suite is not installed in your PATH, you can set up a dedicated Conda environment:\n"
            f"    conda create -n {env_name} -c bioconda mob_suite\n"
            f"Accio will automatically detect and run MOB-suite inside the '{env_name}' environment.\n"
            f"(Or pass --mob_env <env_name> if using a different environment name, or --skip_novel_plasmids)."
        )

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
        cmd = self._resolve_tool('mob_recon') + [
            '-i', input_fasta, 
            '-o', output_dir, '--force', '-u'
        ]
     
        try:
            subprocess.check_call(cmd)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            raise ExternalToolError(self._format_mob_error('MOB-recon', e))
        
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
        cmd = self._resolve_tool('mob_typer') + [
            '-i', input_fasta, 
            '-o', output_file,
            '-n', str(self.config.THREADS)
        ]
        if multi:
            cmd += ['-x']
        if sample_id:
            cmd += ['-s', sample_id]
        if plasmid_meta:
            cmd += ['-m', plasmid_meta,
                    '--plasmid_mash_db', plasmid_mash_db]
        
        try:
            subprocess.check_call(cmd)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            raise ExternalToolError(self._format_mob_error('MOB-typer', e))
        return output_file

    def run_mob_cluster(self, input_fasta: str, output_dir: str,
                         mob_typer_results: str,
                         taxonomy_file: str,
                         previous_mob_db: Optional[str] = None, previous_mob_fastas: Optional[str] = None,
                        primary_cluster_dist: Optional[float] = 0.06, secondary_cluster_dist: Optional[float] = 0.025):
        """
        Run MOB-cluster for plasmid clustering.
        """
        cmd = self._resolve_tool('mob_cluster') + [
            '-m', 'build' if not previous_mob_db else 'update',
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
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            raise ExternalToolError(self._format_mob_error('MOB-cluster', e))



class SkaniRunner:
    """Wrapper for skani operations."""
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        self.config = config or AnalysisConfig()
    
    def run_ska_dist(self, plasmids, output_file):
        cmd = [
            'skani', 'dist', '--ri', str(plasmids), '--qi', str(plasmids),
            '--medium', '-m', '200', '-t', str(self.config.THREADS), '-o', str(output_file)
        ]
        
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(f"skani dist failed: {str(e)}")
                
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

class GenomadRunner:
    """Wrapper for Genomad operations."""
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        self.config = config or AnalysisConfig()
        
    def run_genomad(self, input_fasta: str, output_file: str) -> str:
        """
        Run Genomad for plasmid detection.
        
        Args:
            input_fasta: Path to input FASTA file
            output_file: Output file path
            
        Returns:
            Path to output file
            
        Raises:
            ExternalToolError: If Genomad execution fails
        """
        cmd = [
            'genomad', 'end-to-end', 
            '--threads', str(self.config.THREADS),
        ]

        if self.config.GENOMAD_DISABLE_NN_PRED:
            cmd.append('--disable-nn-classification')

        splits = getattr(self.config, 'GENOMAD_SPLITS', None)
        if splits:
            cmd.extend(['--splits', str(splits)])

        cmd.extend([input_fasta, output_file, self.config.GENOMAD_DB_PATH])
        
        try:
            subprocess.check_call(cmd)
            return output_file
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(f"Genomad failed: {str(e)}")
    
    def read_genomad_results(self, output_file: str) -> pd.DataFrame:
        """
        Read Genomad results from output file.
        
        Args:
            output_file: Path to Genomad output file
            
        Returns:
            DataFrame with Genomad results
        """
        try:
            df = pd.read_csv(output_file, sep='\t')
            return df
        except Exception as e:
            raise ExternalToolError(f"Failed to read Genomad results: {str(e)}")

def check_tool_availability() -> Dict[str, Dict[str, Any]]:
    """
    Check availability of external tools.
    
    Returns:
        Dictionary mapping tool names to their status and installation hints.
    """
    tools = {
        'blastn': {'cmd': 'blastn -version', 'install': 'conda install -c bioconda blast'},
        'makeblastdb': {'cmd': 'makeblastdb -version', 'install': 'conda install -c bioconda blast'},
        'mash': {'cmd': 'mash --version', 'install': 'conda install -c bioconda mash'},
        'nucmer': {'cmd': 'nucmer --version', 'install': 'conda install -c bioconda mummer4'},
        'show-coords': {'cmd': 'show-coords -h', 'install': 'conda install -c bioconda mummer4'},
        'show-snps': {'cmd': 'show-snps -h', 'install': 'conda install -c bioconda mummer4'},
        'delta-filter': {'cmd': 'delta-filter -h', 'install': 'conda install -c bioconda mummer4'},
        'bwa': {'cmd': 'bwa', 'install': 'conda install -c bioconda bwa'},
        'bwa-mem2': {'cmd': 'bwa-mem2', 'install': 'conda install -c bioconda bwa-mem2'},
        'minimap2': {'cmd': 'minimap2 --version', 'install': 'conda install -c bioconda minimap2'},
        'samtools': {'cmd': 'samtools --version', 'install': 'conda install -c bioconda samtools'},
        'skani': {'cmd': 'skani --version', 'install': 'conda install -c bioconda skani'},
        'amrfinder': {'cmd': 'amrfinder --version', 'install': 'conda install -c bioconda ncbi-amrfinderplus'},
        'mob_recon': {'cmd': 'mob_recon --version', 'install': 'conda create -n mob_suite -c bioconda mob_suite (or accio setup-mob)'},
        'mob_typer': {'cmd': 'mob_typer --version', 'install': 'conda create -n mob_suite -c bioconda mob_suite (or accio setup-mob)'},
        'mob_cluster': {'cmd': 'mob_cluster --help', 'install': 'conda create -n mob_suite -c bioconda mob_suite (or accio setup-mob)'},
        'pling': {'cmd': 'pling --help', 'install': 'conda install -c bioconda pling'},
        'plasnet': {'cmd': 'plasnet --help', 'install': 'pip install plasnet'},
        'genomad': {'cmd': 'genomad --version', 'install': 'pip install genomad'},
        'PLASMe.py': {'cmd': 'PLASMe.py --help', 'install': 'See PLASMe GitHub for installation instructions.'},
    }
    
    availability = {}
    
    for tool, info in tools.items():
        availability[tool] = {
            'available': False,
            'install_hint': info['install']
        }
        cmd_parts = shlex.split(info['cmd'])
        binary = cmd_parts[0]
        tool_note = None

        if binary in ['mob_recon', 'mob_typer', 'mob_cluster']:
            mob_runner = MOBRunner()
            resolved = mob_runner._resolve_tool(binary)
            if resolved != [binary]:
                cmd_parts = resolved + cmd_parts[1:]
                binary = cmd_parts[0]
                env_used = getattr(mob_runner.config, 'MOB_CONDA_ENV', 'mob_suite')
                tool_note = f"via conda env {env_used}"

        if not shutil.which(binary):
            continue

        try:
            res = subprocess.run(
                cmd_parts,
                capture_output=True,
                timeout=10
            )
            if res.returncode == 0:
                availability[tool]['available'] = True
                if tool_note:
                    availability[tool]['note'] = tool_note
            elif binary in ['bwa', 'bwa-mem2']:
                combined = res.stdout + res.stderr
                if b'Program: bwa' in combined or b'Usage: bwa' in combined:
                    availability[tool]['available'] = True
            elif binary in ['show-coords', 'show-snps', 'delta-filter']:
                combined = res.stdout + res.stderr
                if binary.encode() in combined or b'USAGE' in combined.upper():
                    availability[tool]['available'] = True
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
            pass
            
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
    db_path = Path(db_dir)

    required_files = {
        'plasmid_info': db_path / DatabaseFiles.PLASMID_INFO,
        'plasmid_fasta': db_path / DatabaseFiles.PLASMID_DB,
        'plasmid_blast': db_path / DatabaseFiles.PLASMID_BLAST,
        'mash_db': db_path / DatabaseFiles.MASH_DB,
        'plasmidfinder_blast': db_path / DatabaseFiles.PLASMIDFINDER_DB,
        'rep_blast': db_path / DatabaseFiles.REP_DB,
    }

    existence = {}
    for name, path in required_files.items():
        existence[name] = path.exists()

    return existence