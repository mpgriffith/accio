"""Core Plasmid class for plasmid analysis."""

import pandas as pd
import numpy as np
import math
import subprocess
import io
import os
import time
from typing import Optional, Dict, List, Any, Union, Tuple
from intervaltree import IntervalTree
from Bio import SeqIO
import pysam

from ..config import AnalysisConfig, DEFAULT_FIELDS, EXTRA_FIELDS
from .contig import Contig


class Plasmid:
    """Represents a reference plasmid with scoring and assignment capabilities."""
    
    def __init__(self, id: str, cluster_id: str, pling_type: str, length: int, plasmid_mash_db: str,
                 other: Optional[Dict[str, Any]] = None, rep_type: Optional[str] = None,
                 amr_genes: Optional[str] = None):
        """
        Initialize a Plasmid object.
        
        Args:
            id: Plasmid identifier
            cluster_id: Cluster identifier for this plasmid
            length: Length of the plasmid
            plasmid_mash_db: Path to Mash database for plasmids
            other: Additional metadata
            rep_type: Replicon type(s) associated with this plasmid
        """
        self.id = str(id)
        self.cluster_id = cluster_id
        self.pling_type = pling_type
        self.length = int(length)
        self.plasmid_db_mash = plasmid_mash_db
        self.rep_type = rep_type
        self.amr_genes = amr_genes
        self.other = other or {}
        
        # Scoring attributes
        self.score = -1.0
        self.mash_score = -1.0
        self.scores = {
            'blast_similarity': -1, 'blast_identity': -1, 'blast_coverage': -1,
            'synteny_score': 0, 'mash_score': -1, 'length_score': 0,
            'read_gaps_score': 0
        }
        
        # Associated contigs and analysis data
        self.contigs: Dict[str, Contig] = {}
        self.blast_data = pd.DataFrame(columns=DEFAULT_FIELDS + EXTRA_FIELDS + ['use_hit'])
        self.synteny_data = pd.DataFrame(columns=[
            'contig', 'plasmid', 'num_orthologous_matches', 'num_matches', 
            'num_breaks', 'match_len'
        ])
        self.mash_data = pd.DataFrame()
        self.nucmer_data = pd.DataFrame()
        
        # Analysis structures
        self.tree = IntervalTree()
        self.synteny_score_type = 'blocks'
        
        # Read coverage data
        self.read_cov_data: Optional[pysam.AlignmentFile] = None
        self.median_read_depth = 30

    def add_rep_type(self, rep_type: str) -> None:
        """Add replicon type information."""
        self.rep_type = rep_type

    def add_contigs(self, contigs: Union[Contig, List[Contig]]) -> None:
        """
        Add contigs to this plasmid and update analysis data.
        
        Args:
            contigs: Single contig or list of contigs to add
        """
        if not isinstance(contigs, list):
            contigs = [contigs]
            
        # Collect data from all contigs
        plasmid_blast_dfs = []
        synteny_scores = []
        
        for contig in contigs:
            self.contigs[contig.id] = contig
            
            # Collect BLAST data specific to this plasmid
            plasmid_blast_data = contig.blast_data[contig.blast_data['sseqid'] == self.id] if hasattr(contig, 'blast_data') and not contig.blast_data.empty else pd.DataFrame()
            if not plasmid_blast_data.empty:
                plasmid_blast_dfs.append(plasmid_blast_data)
                
            # Collect synteny data
            synteny_score_c = contig.synteny_data[contig.synteny_data['plasmid'] == self.id] if hasattr(contig, 'synteny_data') and not contig.synteny_data.empty else pd.DataFrame()
            if not synteny_score_c.empty:
                synteny_scores.append(synteny_score_c)
        
        # Combine BLAST data
        all_blast_dfs = [df for df in [self.blast_data] + plasmid_blast_dfs 
                        if not df.dropna(how='all').empty]
        
        if all_blast_dfs:
            self.blast_data = pd.concat(all_blast_dfs, ignore_index=True)
        
        # Combine synteny data
        if synteny_scores:
            self.synteny_data = pd.concat([self.synteny_data] + synteny_scores, ignore_index=True)
            
        # Recalculate scores
        self.calc_blast_score()
        self.calc_scores()

    def remove_contigs(self, contig_ids: Union[str, List[str]]) -> None:
        """
        Remove contigs from this plasmid and update scores.
        
        Args:
            contig_ids: Single contig ID or list of contig IDs to remove
        """
        if isinstance(contig_ids, str):
            contig_ids = [contig_ids]
            
        if not contig_ids:
            return
            
        # Remove from contigs dict
        for cid in contig_ids:
            if cid in self.contigs:
                del self.contigs[cid]
                
        # Filter data to exclude removed contigs
        self.blast_data = self.blast_data[~self.blast_data['qseqid'].isin(contig_ids)]
        self.synteny_data = self.synteny_data[~self.synteny_data['contig'].isin(contig_ids)]
        
        # Recalculate scores
        self.calc_blast_score()
        self.calc_scores()

    def add_mash_data(self, mash_data: pd.DataFrame) -> None:
        """
        Add Mash screening results.
        
        Args:
            mash_data: DataFrame containing Mash results
        """
        if mash_data.empty:
            self.mash_data = pd.DataFrame()
            self.mash_score = -1
            self.scores['mash_score'] = -1
            return
            
        self.mash_data = mash_data.iloc[0]
        self.mash_score = mash_data['pident'].mean()
        self.scores['mash_score'] = self.mash_score
        
        # Extract additional Mash metrics
        if not mash_data.empty:
            row = mash_data.iloc[0]
            self.scores['mash_score_shared_hashes'] = int(row['shared-hashes'].split('/')[0])
            self.scores['mash_score_median_mult'] = int(row['median-mult'])
            self.scores['mash_score_pvalue'] = float(row['pvalue'])

    def add_read_data(self, mash_data: pd.DataFrame, cov_data: Optional[pysam.AlignmentFile], 
                     median_depth: Optional[float] = None) -> None:
        """
        Add read-based analysis data.
        
        Args:
            mash_data: Mash screening results
            cov_data: Coverage data from read alignment
            median_depth: Median read depth for normalization
        """
        if not mash_data.empty:
            self.add_mash_data(mash_data)
            
        if median_depth:
            self.median_read_depth = median_depth
            
        self.read_cov_data = cov_data

    def calc_blast_score(self) -> None:
        """Calculate BLAST-based similarity scores."""
        if self.blast_data.empty:
            self.scores.update({
                'blast_identity': 0,
                'blast_coverage': 0,
                'blast_similarity': 0
            })
            return
            
        plasmid_cover = [0] * (int(self.length) + 1)
        self.blast_data['use_hit'] = False
        keep_rows = []
        
        # Sort by quality metrics
        blast_data = self.blast_data.sort_values(
            ['bitscore', 'pident', 'length'], 
            ascending=[False, False, False]
        )
        
        for idx, row in blast_data.iterrows():
            plasmid_range = sorted([int(row['sstart']), int(row['send'])])
            
            # Check overlap with already covered regions
            if plasmid_range[1] > len(plasmid_cover):
                plasmid_cover.extend([0] * (plasmid_range[1] - len(plasmid_cover)))
                
            pct_covered = sum(plasmid_cover[i-1] for i in range(
                plasmid_range[0], min(plasmid_range[1], len(plasmid_cover))
            )) / (plasmid_range[1] - plasmid_range[0])
            
            if pct_covered > 0.95:
                continue
                
            # Mark region as covered
            for i in range(plasmid_range[0], plasmid_range[1]):
                if i-1 < len(plasmid_cover):
                    plasmid_cover[i-1] = 1
                    
            keep_rows.append(idx)
            
        self.blast_data.loc[keep_rows, 'use_hit'] = True
        self.blast_data = self.blast_data.loc[keep_rows]
        
        # Calculate interval tree for gap finding
        self.calc_intervals()
        
        # Calculate scores
        used_hits = self.blast_data[self.blast_data['use_hit'] == True]
        if not used_hits.empty:
            self.scores['blast_identity'] = (
                sum(used_hits['pident'] * used_hits['length']) / sum(used_hits['length'])
            )
        else:
            self.scores['blast_identity'] = 0
            
        self.scores['blast_coverage'] = (sum(plasmid_cover) / len(plasmid_cover)) * 100
        self.scores['blast_similarity'] = (
            (sum(plasmid_cover) / len(plasmid_cover)) * self.scores['blast_identity']
        )

    def calc_intervals(self, data: Optional[pd.DataFrame] = None) -> IntervalTree:
        """
        Create an interval tree from alignment data.
        
        Args:
            data: DataFrame to create intervals from (uses self.blast_data if None)
            
        Returns:
            IntervalTree object
        """
        if data is None:
            df = self.blast_data[self.blast_data['use_hit'] == True]
            tuples = [sorted((int(r['sstart']), int(r['send']))) for idx, r in df.iterrows()]
            self.tree = IntervalTree.from_tuples(tuples)
            return self.tree
        else:
            tuples = []
            for idx, row in data.iterrows():
                pos = sorted((int(row['sstart']), int(row['send'])))
                tuples.append((pos[0], pos[1], row))
            return IntervalTree.from_tuples(tuples)

    def find_gaps(self) -> List[Tuple[int, int]]:
        """
        Find gaps in plasmid coverage.
        
        Returns:
            List of (start, end) tuples representing uncovered regions
        """
        gaps = []
        last_end = 0
        
        for interval in sorted(self.tree):
            if last_end < interval.begin:
                gaps.append((last_end, interval.begin))
            last_end = max(last_end, interval.end)
            
        if last_end < self.length:
            gaps.append((last_end, self.length))
            
        return gaps

    def calc_synteny_score(self, synteny_type: Optional[str] = None, 
                          gap_threshold: int = 200) -> None:
        """
        Calculate synteny-based similarity scores.
        
        Args:
            synteny_type: Type of synteny analysis ('nucmer' or 'blocks')
            gap_threshold: Gap size threshold for synteny breaks
        """
        if synteny_type is None:
            synteny_type = self.synteny_score_type
            
        # Initialize data structure
        data = {
            'orthologous_matches': 0, 'num_contigs': 0, 'num_breaks': 0,
            'match_len': 0, 'num_matches': 0, 'score': 0
        }
        
        try:
            # Collect nucmer data from all contigs
            nucmer_dfs = [
                c.nucmer_data[c.nucmer_data['sseqid'] == self.id] 
                for c in self.contigs.values() 
                if hasattr(c, 'nucmer_data') and not c.nucmer_data.empty
            ]
            
            if not nucmer_dfs:
                self.scores.update({'synteny_' + k: v for k, v in data.items()})
                self.scores['synteny_score'] = 0
                return
                
            df = pd.concat(nucmer_dfs)
            
        except Exception:
            self.scores.update({'synteny_' + k: v for k, v in data.items()})
            self.scores['synteny_score'] = 0
            return
            
        if synteny_type == 'nucmer':
            self._calc_nucmer_synteny(df, data)
        elif synteny_type == 'blocks':
            self._calc_block_synteny(df, data)
            
        self.scores.update({'synteny_' + k: v for k, v in data.items()})

    def _calc_nucmer_synteny(self, df: pd.DataFrame, data: Dict[str, Any]) -> None:
        """Calculate synteny score using nucmer-style binning approach."""
        from ..utils.filters import filter_overlaps
        
        bin_size = 500
        df = filter_overlaps(
            df, by='sseqid', 
            sort_cols=['pident', 'ref_aln_len'], 
            sort_ascending=[False, False]
        )
        df = df[df['pident'] > 80].reset_index()
        df = df.sort_values(by='sstart')
        
        # Create genomic bins
        num_bins = int(self.length / bin_size)
        bins_covered = np.zeros(num_bins, dtype=bool)
        
        # Mark covered bins
        for _, row in df.iterrows():
            start_bin = int(min(row['sstart'], row['send']) / bin_size)
            end_bin = int(max(row['sstart'], row['send']) / bin_size)
            bins_covered[start_bin:end_bin + 1] = True
            
        block_coverage = bins_covered.sum() / num_bins
        
        # Count contigs and calculate continuity penalty
        num_contigs = df['qseqid'].nunique()
        continuity_penalty = math.log2(num_contigs + 1)
        
        # Calculate synteny conservation
        df_sorted = df.sort_values(by=['qseqid', 'qstart'])
        conserved_adjacent = 0
        total_adjacent = 0
        prev_row = None
        
        for _, row in df_sorted.iterrows():
            if (prev_row is not None and 
                row['qseqid'] == prev_row['qseqid']):
                
                if ((row['sstart'] > prev_row['sstart'] and row['qstart'] > prev_row['qstart']) or
                    (row['sstart'] < prev_row['sstart'] and row['qstart'] < prev_row['qstart'])):
                    conserved_adjacent += 1
                total_adjacent += 1
                
            prev_row = row
            
        adjacency_score = conserved_adjacent / total_adjacent if total_adjacent > 0 else 0.5
        
        # Check for circular features
        wrap_hits = df[
            (df['qstart'] > df['qend']) &
            (df['qstart'] > 0.9 * df['qlen']) &
            (df['qend'] < 0.1 * df['qlen'])
        ]
        circular_bonus = 1.2 if len(wrap_hits) > 0 else 1.0
        
        # Final score calculation
        score = (block_coverage * adjacency_score / continuity_penalty) * circular_bonus
        
        data.update({
            'block_length': block_coverage,
            'num_contigs': num_contigs,
            'gap_penalty': continuity_penalty,
            'orthologous_matches': adjacency_score,
            'score': score / 1.2
        })

    def _calc_block_synteny(self, df: pd.DataFrame, data: Dict[str, Any]) -> None:
        """Calculate synteny score using block-based approach."""
        from ..utils.filters import filter_overlaps
        
        df = filter_overlaps(
            df, by='sseqid',
            sort_cols=['ref_aln_len', 'pident'],
            sort_ascending=[False, False]
        )
        
        tree = self.calc_intervals(df)
        blocks = [(i, i + 500) for i in range(0, self.length, 500)]
        
        prev_block = None
        prev_contigs = set()
        
        for block_start, block_end in blocks:
            intervals = tree[block_start:block_end]
            
            if not intervals:
                continue
                
            data['num_matches'] += 1
            
            for interval in sorted(intervals):
                row_data = interval.data
                
                if isinstance(prev_block, pd.Series) and prev_block.name == row_data.name:
                    data['orthologous_matches'] += 1
                elif isinstance(prev_block, pd.Series):
                    if row_data['qseqid'] in prev_contigs:
                        contig = self.contigs.get(row_data['qseqid'])
                        
                        if (contig and contig.contig_type == 'repetitive' and 
                            row_data.get('qcovs', 0) > 80):
                            data['orthologous_matches'] += 1
                        elif self._check_circular_continuity(prev_block, row_data):
                            data['orthologous_matches'] += 1
                        else:
                            data['num_breaks'] += 1
                            
                prev_block = row_data
                prev_contigs.add(row_data['qseqid'])
                
        # Calculate final scores
        num_blocks = len(blocks)
        if num_blocks > 0:
            coverage_score = data['num_matches'] / num_blocks
            ortho_score = data['orthologous_matches'] / num_blocks
            gap_penalty = data['num_breaks'] / num_blocks
            data['score'] = (coverage_score - gap_penalty + ortho_score) / 2
        
        # Calculate length score
        data['num_contigs'] = df['qseqid'].nunique() if not df.empty else 0
        self._calc_length_score(df, data)

    def _check_circular_continuity(self, prev_block: pd.Series, curr_block: pd.Series) -> bool:
        """Check if two blocks represent circular continuity."""
        if prev_block['qseqid'] != curr_block['qseqid']:
            return False
            
        if prev_block['qframe'] != curr_block['qframe']:
            return False
            
        # Check for wraparound at contig ends
        return (
            (prev_block['qend'] < 10 and curr_block['qstart'] > curr_block['qlen'] - 10) or
            (prev_block['qstart'] > curr_block['qlen'] - 10 and curr_block['qend'] < 10)
        )

    def _calc_length_score(self, df: pd.DataFrame, data: Dict[str, Any]) -> None:
        """Calculate length-based similarity score."""
        if df.empty:
            data['length_score'] = 0
            return
            
        contig_sum = 0
        rel_contig_depths = []
        
        # Calculate relative contig depths
        for contig_id in df['qseqid'].unique():
            contig = self.contigs.get(contig_id)
            if contig and contig.contig_type != 'repetitive' and contig.copy_num:
                rel_contig_depths.append(contig.copy_num)
                
        rel_contig_depth = (
            sum(rel_contig_depths) / len(rel_contig_depths) 
            if rel_contig_depths else 1
        )
        
        # Calculate expected contig contribution
        for contig_id, contig_df in df.groupby('qseqid'):
            contig = self.contigs.get(contig_id)
            if not contig:
                continue
                
            # Calculate coverage
            contig_cov = [0] * contig.length
            num_hits_80 = 0
            
            for _, row in contig_df.iterrows():
                cstart, cend = sorted([int(row['qstart']), int(row['qend'])])
                for i in range(cstart, min(cend, len(contig_cov))):
                    contig_cov[i] = 1
                    
                if row.get('qcovs', 0) > 80:
                    num_hits_80 += 1
                    
            contig_cov_pct = sum(contig_cov) / len(contig_cov) * 100 if contig_cov else 0
            
            if contig_cov_pct > 80:
                if contig.contig_type == 'repetitive' and num_hits_80 > 0:
                    contig_depth = contig.copy_num / rel_contig_depth if contig.copy_num else num_hits_80
                    contig_sum += len(contig_cov) * min(num_hits_80, round(contig_depth))
                else:
                    contig_sum += len(contig_cov)
                    
        # Handle single contig case
        if len(self.contigs) == 1:
            contig_sum = sum(c.length for c in self.contigs.values())
            
        # Calculate length similarity score
        if max(contig_sum, self.length) > 0:
            data['length_score'] = 1 - abs(contig_sum - self.length) / max(contig_sum, self.length)
        else:
            data['length_score'] = 0
            
        self.scores['length_score'] = data['length_score']
        self.other['length_contigs'] = contig_sum


    def calc_mash_score(self) -> float:
        """
        Calculate Mash similarity score by screening contigs against database.
        
        Returns:
            Mash similarity score
        """
        contig_seqs = [c.seq for c in self.contigs.values()]
        if not contig_seqs:
            return 0.0
            
        # Create temporary file for contigs
        temp_filename = f"{self.id}_{int(time.time())}_remaining_contigs.fasta"
        
        try:
            SeqIO.write(contig_seqs, temp_filename, 'fasta')
            
            # Run mash screen
            cmd = ['mash', 'screen', self.plasmid_db_mash, temp_filename]
            result = subprocess.check_output(cmd, stderr=subprocess.DEVNULL)
            
            # Parse results
            result_io = io.StringIO(result.decode())
            df = pd.read_csv(result_io, sep='\t', header=None, names=[
                'pident', 'shared-hashes', 'median-mult', 'pvalue', 'sseqid', 'stitle'
            ])
            
            # Find score for this plasmid
            score_row = df[df['sseqid'] == self.id]
            if score_row.empty:
                self.mash_score = 0
                self.scores['mash_score'] = 0
                self.scores['mash_sim_score'] = 0
                self.scores['mash_score_shared_hashes'] = 0
                self.scores['mash_score_median_mult'] = 0
                self.scores['mash_score_pvalue'] = 0
 
                return 0.0
            score_row = score_row.iloc[0]
            score_row['hash_pct'] = eval(score_row['shared-hashes'])
            self.scores['mash_score'] = score_row.pident 
            self.scores['mash_sim_score'] = score_row.pident * score_row.hash_pct
            self.scores['mash_score_shared_hashes'] = int(score_row['shared-hashes'].split('/')[0])
            self.scores['mash_score_median_mult'] = int(score_row['median-mult'])
            self.scores['mash_score_pvalue'] = float(score_row['pvalue'])
            return float(score_row.pident)
            
        except subprocess.CalledProcessError:
            return 0.0
        finally:
            if os.path.exists(temp_filename):
                os.remove(temp_filename)

    def calc_read_score(self) -> Tuple[float, float, float, float]:
        """
        Calculate read-based gap coverage score.
        
        Returns:
            Tuple of (score, pct_coverage, gap_reads, gap_length)
        """
        if not self.read_cov_data:
            return 0.0, 0.0, 0.0, 0.0
            
        gaps = self.find_gaps()
        gap_data = []
        
        for start, end in gaps:
            gap_len = end - start
            if gap_len == 0:
                continue
                
            try:
                gap_pileup = self.read_cov_data.pileup(self.id, start, end)
                gap_depths = [len(list(pileup.pileups)) for pileup in gap_pileup]
                
                if gap_depths:
                    min_depth = self.median_read_depth * 0.25
                    gap_cov_pct = sum(1 for depth in gap_depths if depth > min_depth) / gap_len
                    avg_depth = sum(gap_depths) / len(gap_depths)
                else:
                    gap_cov_pct = 0.0
                    avg_depth = 0.0
                    
                gap_data.append([start, end, gap_len, gap_cov_pct, avg_depth])
                
            except Exception:
                # Handle cases where pileup fails
                gap_data.append([start, end, gap_len, 0.0, 0.0])
                
        if not gap_data:
            return 1.0, 1.0, 0.0, 0.0
            
        gap_df = pd.DataFrame(gap_data, columns=['start', 'end', 'gap_len', 'cov_pct', 'gap_reads'])
        
        # Calculate weighted averages
        total_gap_len = gap_df['gap_len'].sum()
        if total_gap_len == 0:
            return 1.0, 1.0, 0.0, 0.0
            
        avg_depth = (gap_df['gap_reads'] * gap_df['gap_len']).sum() / total_gap_len
        pct_cov = (gap_df['cov_pct'] * gap_df['gap_len']).sum() / total_gap_len
        
        # Calculate final score
        covered_len = self.length - total_gap_len
        score = (covered_len * 1.0 + pct_cov * total_gap_len) / self.length
        
        return score, pct_cov, avg_depth, total_gap_len

    def calc_scores(self, score_method: Optional[str] = None, 
                   weights: Optional[Dict[str, float]] = None,
                   config: Optional[AnalysisConfig] = None) -> None:
        """
        Calculate overall similarity score for this plasmid.
        
        Args:
            score_method: Scoring method to use
            weights: Custom scoring weights
            config: Analysis configuration
        """
        if config is None:
            config = AnalysisConfig()
            
        # Determine weights based on data availability
        if weights is None:
            if self.read_cov_data is not None:
                weights = config.SCORING_WEIGHTS_WITH_READS
            else:
                weights = config.SCORING_WEIGHTS_NO_READS
                
        # Calculate component scores
        self.calc_synteny_score()
        mash_score = self.calc_mash_score()
        self.scores['mash_score'] = mash_score
        self.mash_score = mash_score
        
        # Calculate read-based scores if available
        if self.read_cov_data is not None:
            read_score, pct_cov, gap_reads, gap_len = self.calc_read_score()
            self.scores.update({
                'read_gaps_score': read_score,
                'num_gap_reads': gap_reads,
                'gap_len_reads': gap_len,
                'gap_pct_cov': pct_cov
            })
        else:
            self.scores.update({
                'read_gaps_score': 0,
                'num_gap_reads': 0,
                'gap_len_reads': 0
            })
            
        # Calculate weighted final score
        self.score = (
            weights['synteny_score'] * self.scores['synteny_score'] +
            weights['mash_score'] * mash_score +
            weights['blast_identity'] * self.scores['blast_identity'] / 100 +
            weights['blast_coverage'] * self.scores['blast_coverage'] / 100 +
            weights['length_score'] * self.scores['length_score'] +
            weights['read_gaps_score'] * self.scores['read_gaps_score']
        )

    def get_contigs_from_blast(self) -> List[str]:
        """
        Get list of contigs with significant BLAST hits.
        
        Returns:
            List of contig IDs with usable BLAST hits
        """
        if self.blast_data.empty:
            return []
        return list(self.blast_data[self.blast_data['use_hit'] == True]['qseqid'].unique())

    def get_summary(self) -> Dict[str, Any]:
        """
        Get comprehensive summary of this plasmid's analysis.
        
        Returns:
            Dictionary containing plasmid summary information
        """
        return {
            'plasmid_id': self.id,
            'pling_type': self.pling_type,
            'cluster_id': self.cluster_id,
            'length': self.length,
            'rep_type': self.rep_type,
            'amr_genes': self.amr_genes,
            'num_contigs': len(self.contigs),
            'contig_ids': list(self.contigs.keys()),
            'final_score': self.score,
            **self.scores,
            **self.other
        }

    def __repr__(self) -> str:
        """String representation of the plasmid."""
        contigs_str = ', '.join(self.contigs.keys()) if self.contigs else 'none'
        return f"Plasmid({self.id}, score={self.score:.3f}, contigs=[{contigs_str}])"