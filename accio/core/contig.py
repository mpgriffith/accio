"""Core Contig class for plasmid analysis."""

import pandas as pd
from typing import Optional, List, Dict, Any, Union
from Bio.SeqRecord import SeqRecord

from ..config import AnalysisConfig, MAX_PLASMID_LENGTH
from ..utils.filters import filter_overlaps


class Contig:
    """Represents a genomic contig with associated analysis data."""
    
    def __init__(self, id: str, seq: SeqRecord, length: Optional[int] = None, circular: bool = False):
        """
        Initialize a Contig object.
        
        Args:
            id: Contig identifier
            seq: SeqRecord object containing sequence data
            length: Length of the contig (computed from seq if not provided)
            circular: Whether the contig is circular
        """
        self.id = str(id)
        self.seq = seq
        self.length = length if length is not None else len(seq)
        self.circular = circular

        # Classification attributes
        self.contig_type: Optional[str] = None
        self.type: Optional[str] = True if self.circular else None
        self.assigned: Union[bool, str] = False
        
        # Analysis data storage
        self.matches: Dict[str, Any] = {}
        self.gene_data: Dict[str, Any] = {}
        self.plasmid_matches: List[str] = []
        self.rep_types: Optional[List[str]] = None
        self.copy_num: Optional[float] = None
        self.plas_data: List[Optional[str]] = [(None, None)]  # [plasmid_id, cluster_id]
        
        # Analysis results DataFrames
        self.blast_data = pd.DataFrame()
        self.rep_data = pd.DataFrame()
        self.fastani_data = pd.DataFrame(columns=[
            'plasmid', 'orthologous_matches', 'num_matches', 'num_breaks', 'match_len'
        ])
        self.nucmer_data = pd.DataFrame()
        self.synteny_data = pd.DataFrame(columns=[
            'plasmid', 'orthologous_matches', 'num_matches', 'num_breaks', 'match_len'
        ])
        self.mash_data = pd.DataFrame()

    def add_blast_data(self, df: pd.DataFrame, config: Optional[AnalysisConfig] = None) -> None:
        """
        Add and filter BLAST results for this contig.
        
        Args:
            df: DataFrame containing BLAST results
            config: Analysis configuration (uses default if None)
        """
        if config is None:
            config = AnalysisConfig()
            
        if df.empty:
            self.blast_data = pd.DataFrame()
            return
            
        keep_idx = []
        contig_cover = {}
        
        for idx, row in df.iterrows():
            if row['pident'] < config.MIN_IDENTITY:
                continue
                
            if self.circular:
                plas = row['sseqid']
                if plas not in contig_cover:
                    contig_cover[plas] = [0] * self.length
                    
                cstart, cend = int(row['qstart']), int(row['qend'])
                pct_cover_contig = sum(1 for i in contig_cover[plas][cstart:cend] if i >= 1) / abs(cend - cstart)
                
                if pct_cover_contig < 0.9:
                    keep_idx.append(idx)
                    for i in range(cstart, cend):
                        if i < len(contig_cover[plas]):
                            contig_cover[plas][i] += 1
                            
            elif row['qcovs'] >= config.MIN_COVERAGE:
                keep_idx.append(idx)
            elif row['qstart'] < 50 or row['qlen'] - row['qend'] <= 50:
                keep_idx.append(idx)
        
        if keep_idx:
            filtered_df = filter_overlaps(df.loc[keep_idx], by='sseqid')
            self.blast_data = filtered_df
        else:
            self.blast_data = pd.DataFrame()

    def add_rep_type(self, df: pd.DataFrame) -> None:
        """
        Add replicon type information from PlasmidFinder results.
        
        Args:
            df: DataFrame containing PlasmidFinder results
        """
        if not df.empty:
            self.type = 'plasmid'
            self.plas_data = [(None, None)]
            self.rep_types = list(df['sseqid'].unique())

    def add_plasme_data(self, plasme_df: pd.DataFrame) -> None:
        """
        Add PLASMe classification results.
        
        Args:
            plasme_df: DataFrame containing PLASMe results
        """
        if plasme_df.empty:
            return
            
        row = plasme_df.iloc[0]
        
        # Check for ambiguous regions
        amb_regions = row.get('amb_region', '')
        if isinstance(amb_regions, str) and amb_regions:
            amb_region_len = 0
            for r in amb_regions.split(','):
                if '-' in r:
                    s, e = r.split('-')
                    try:
                        amb_region_len += int(e) - int(s)
                    except ValueError:
                        continue
                        
            if amb_region_len / row['length'] > 0.5:
                return
                
        # Check transformer evidence
        if (row.get('evidence') == 'Transformer' and 
            row.get('score', 0) <= 0.75):
            return
            
        self.type = 'plasmid'
        self.plas_data = [(None, None)]

    def add_repetitive_data(self, rep_df: pd.DataFrame, config: Optional[AnalysisConfig] = None) -> None:
        """
        Add repetitive element analysis results.
        
        Args:
            rep_df: DataFrame containing repetitive element BLAST results
            config: Analysis configuration
        """
        if config is None:
            config = AnalysisConfig()
            
        if rep_df.empty:
            self.rep_data = pd.DataFrame()
            return
            
        # Filter by identity threshold
        rep_df = rep_df[rep_df['pident'] > config.REP_MIN_IDENT].copy()
        
        if rep_df.empty:
            self.rep_data = pd.DataFrame()
            return
            
        rep_df['rep_elem_cov'] = 0.0
        rep_df['contig_cov'] = 0.0
        contig_cov = [0] * (int(self.length) + 1)
        
        for rep_elem, df_group in rep_df.groupby('sseqid'):
            rep_elem_cov = [0] * (int(df_group.iloc[0]['slen']) + 1)
            
            for idx, r in df_group.iterrows():
                # Mark covered positions in repeat element
                for i in range(int(r['sstart']), int(r['send'])):
                    if i - 1 < len(rep_elem_cov):
                        rep_elem_cov[i - 1] = 1
                        
                # Mark covered positions in contig
                for i in range(int(r['qstart']), int(r['qend'])):
                    if i - 1 < len(contig_cov):
                        contig_cov[i - 1] = 1
                        
            coverage_pct = sum(rep_elem_cov) / len(rep_elem_cov) * 100
            rep_df.loc[df_group.index, 'rep_elem_cov'] = coverage_pct
            
        contig_coverage_pct = sum(contig_cov) / len(contig_cov) * 100
        rep_df['contig_cov'] = contig_coverage_pct
        
        # Filter by coverage thresholds
        rep_df = rep_df[
            (rep_df['rep_elem_cov'] > config.REP_MIN_COV) | 
            (rep_df['contig_cov'] > config.REP_MIN_COV)
        ]
        
        if not rep_df.empty:
            self.contig_type = 'repetitive'
            
        self.rep_data = rep_df

    def add_copy_num(self, copy_num: float, config: Optional[AnalysisConfig] = None) -> None:
        """
        Add copy number information.
        
        Args:
            copy_num: Estimated copy number from coverage analysis
            config: Analysis configuration
        """
        if config is None:
            config = AnalysisConfig()
            
        self.copy_num = copy_num
        if copy_num > config.HIGH_COPY_THRESHOLD:
            self.contig_type = 'repetitive'

    def add_mash_data(self, mash_df: pd.DataFrame) -> None:
        """Add Mash screening results."""
        self.mash_data = mash_df

    def add_nucmer_data(self, nucmer_df: pd.DataFrame, snps_df: pd.DataFrame) -> pd.DataFrame:
        """
        Add Nucmer alignment results.
        
        Args:
            nucmer_df: DataFrame containing Nucmer coordinate results
            snps_df: DataFrame containing SNP data
            
        Returns:
            The processed nucmer data
        """
        if not nucmer_df.empty:
            # Filter overlapping alignments
            filtered_nucmer = filter_overlaps(
                nucmer_df, 
                by='sseqid', 
                sort_cols=['pident', 'ref_aln_len'], 
                sort_ascending=[False, False]
            )
            filtered_nucmer['plasmid'] = filtered_nucmer['sseqid']
            self.nucmer_data = filtered_nucmer
            self.synteny_data = self.nucmer_data
        else:
            self.nucmer_data = pd.DataFrame()
            self.synteny_data = pd.DataFrame()
            
        return self.nucmer_data

    def assign(self, plasmid_id: str, cluster_id: str) -> None:
        """
        Assign this contig to a plasmid.
        
        Args:
            plasmid_id: ID of the assigned plasmid
            cluster_id: Cluster ID of the assigned plasmid
        """
        self.assigned = 'plasmid'
        if self.plas_data[0] is not None and self.contig_type == 'repetitive':
            self.plas_data = self.plas_data + [(plasmid_id, cluster_id)]
        else:
            self.plas_data = [(plasmid_id, cluster_id)]

    def is_plasmid_candidate(self) -> bool:
        """
        Check if this contig is a candidate for plasmid assignment.
        
        Returns:
            True if contig shows evidence of being a plasmid
        """
        
        return (
            self.type == 'plasmid' or
            self.length <= MAX_PLASMID_LENGTH or  # Size threshold for plasmids
            self.circular or
            bool(self.rep_types) or
            True # true for now while i think about integrated plasmids
        )

    def get_summary(self) -> Dict[str, Any]:
        """
        Get a summary of this contig's properties and analysis results.
        
        Returns:
            Dictionary containing contig summary information
        """
        return {
            'contig_id': self.id,
            'length': self.length,
            'circular': self.circular,
            'contig_type': self.contig_type,
            'type': self.type,
            'assigned': self.assigned,
            'assigned plasmid(s)': [p[0] for p in self.plas_data] if len(self.plas_data) > 0 else None,
            'assigned cluster(s)': [p[1] for p in self.plas_data] if len(self.plas_data) > 0 else None,
            'copy_num': self.copy_num,
            'rep_types': self.rep_types,
            'num_blast_hits': len(self.blast_data) if not self.blast_data.empty else 0,
            'num_plasmid_matches': len(self.plasmid_matches),
            'is_plasmid_candidate': self.is_plasmid_candidate()
        }

    def __repr__(self) -> str:
        """String representation of the contig."""
        status = f"assigned to {', '.join(self.plas_data[i][0] for i in range(len(self.plas_data)))}" if self.assigned else "unassigned"
        return f"Contig({self.id}, {self.length}bp, {status})"