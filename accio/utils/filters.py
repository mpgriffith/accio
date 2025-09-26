"""Filtering utilities for plasmid analysis."""

import pandas as pd
from typing import List, Optional


def filter_overlaps(df: pd.DataFrame, 
                   by: str = 'qseqid', 
                   pos_cols: Optional[List[str]] = None,
                   sort_cols: Optional[List[str]] = None, 
                   sort_ascending: Optional[List[bool]] = None,
                   overlap_threshold: float = 0.95) -> pd.DataFrame:
    """
    Filter overlapping alignments to keep only the best non-overlapping hits.
    
    Args:
        df: DataFrame containing alignment results
        by: Column to group by when filtering overlaps
        pos_cols: Columns containing start/end positions [start, end]
        sort_cols: Columns to sort by for prioritization
        sort_ascending: Sort order for each sort column
        overlap_threshold: Maximum allowed overlap fraction
        
    Returns:
        Filtered DataFrame with overlaps removed
    """
    if df.empty:
        return df
        
    # Set default position columns based on grouping column
    if pos_cols is None:
        if by == 'sseqid':
            pos_cols = ['sstart', 'send']
        else:
            pos_cols = ['qstart', 'qend']
            
    # Set default sort columns and order
    if sort_cols is None:
        sort_cols = ['bitscore', 'pident']
    if sort_ascending is None:
        sort_ascending = [False] * len(sort_cols)
    elif len(sort_cols) != len(sort_ascending):
        sort_ascending = [False] * len(sort_cols)
        
    keep_idx = []
    
    for group_val, group_df in df.groupby(by):
        # Determine coverage length for this group
        if 'qlen' in df.columns and by == 'qseqid':
            cov_len = group_df.iloc[0]['qlen']
        else:
            cov_len = group_df[pos_cols[-1]].max()
            
        if pd.isna(cov_len) or cov_len <= 0:
            cov_len = 1000  # Default fallback
            
        # Initialize coverage tracking
        covered = [0] * (int(cov_len) + 1)
        
        # Sort hits by priority
        group_df = group_df.sort_values(by=sort_cols, ascending=sort_ascending)
        
        for idx, row in group_df.iterrows():
            # Get alignment positions
            try:
                start_pos = int(row[pos_cols[0]])
                end_pos = int(row[pos_cols[1]])
                range_pos = sorted([start_pos, end_pos])
            except (ValueError, KeyError):
                # Skip rows with invalid position data
                continue
                
            hit_len = range_pos[1] - range_pos[0]
            if hit_len <= 0:
                continue
                
            # Extend coverage array if necessary
            if len(covered) < range_pos[1]:
                covered.extend([0] * (range_pos[1] - len(covered)))
                
            # Check overlap with already covered regions
            try:
                overlap_bases = sum(
                    covered[i-1] for i in range(range_pos[0], min(range_pos[1] + 1, len(covered)))
                    if i > 0
                )
                overlap_fraction = overlap_bases / hit_len if hit_len > 0 else 1.0
                
                # Skip if too much overlap
                if overlap_fraction > overlap_threshold:
                    continue
                    
                # Mark region as covered
                for i in range(range_pos[0], min(range_pos[1] + 1, len(covered))):
                    if i > 0:
                        covered[i-1] = 1
                        
                keep_idx.append(idx)
                
            except IndexError:
                # Handle edge cases with position indexing
                continue
                
    return df.loc[keep_idx] if keep_idx else pd.DataFrame(columns=df.columns)


def filter_blast_results(df: pd.DataFrame, 
                        min_identity: float = 75, 
                        min_coverage: float = 75,
                        min_length: Optional[int] = None) -> pd.DataFrame:
    """
    Filter BLAST results by quality metrics.
    
    Args:
        df: DataFrame containing BLAST results
        min_identity: Minimum percent identity threshold
        min_coverage: Minimum query coverage threshold
        min_length: Minimum alignment length threshold
        
    Returns:
        Filtered DataFrame
    """
    if df.empty:
        return df
        
    filtered = df.copy()
    
    # Apply identity filter
    if 'pident' in filtered.columns:
        filtered = filtered[filtered['pident'] >= min_identity]
        
    # Apply coverage filter
    if 'qcovs' in filtered.columns:
        filtered = filtered[filtered['qcovs'] >= min_coverage]
        
    # Apply length filter
    if min_length is not None and 'length' in filtered.columns:
        filtered = filtered[filtered['length'] >= min_length]
        
    return filtered


def filter_by_similarity_score(df: pd.DataFrame, 
                              min_score: float = 80,
                              score_formula: str = 'pident * length / slen') -> pd.DataFrame:
    """
    Filter alignments by a computed similarity score.
    
    Args:
        df: DataFrame containing alignment data
        min_score: Minimum similarity score threshold
        score_formula: Formula for computing similarity score
        
    Returns:
        Filtered DataFrame with similarity scores
    """
    if df.empty:
        return df
        
    try:
        # Calculate similarity score
        df = df.copy()
        df['sim_score'] = df.eval(score_formula)
        
        # Filter by threshold
        filtered = df[df['sim_score'] >= min_score]
        
        return filtered
        
    except Exception:
        # Return original DataFrame if scoring fails
        return df


def filter_plasmid_candidates(df: pd.DataFrame, 
                             max_length: int = 300000,
                             require_circular: bool = False) -> pd.DataFrame:
    """
    Filter contigs that are candidates for plasmid assignment.
    
    Args:
        df: DataFrame containing contig information
        max_length: Maximum length for plasmid candidates
        require_circular: Whether to require circular contigs only
        
    Returns:
        Filtered DataFrame containing only plasmid candidates
    """
    if df.empty:
        return df
        
    filtered = df.copy()
    
    # Apply length filter
    if 'length' in filtered.columns:
        filtered = filtered[filtered['length'] <= max_length]
        
    # Apply circularity filter if required
    if require_circular and 'circular' in filtered.columns:
        filtered = filtered[filtered['circular'] == True]
        
    return filtered


def remove_low_quality_alignments(df: pd.DataFrame,
                                 min_bitscore: Optional[float] = None,
                                 max_evalue: Optional[float] = None) -> pd.DataFrame:
    """
    Remove alignments with poor quality metrics.
    
    Args:
        df: DataFrame containing alignment results
        min_bitscore: Minimum bit score threshold
        max_evalue: Maximum e-value threshold
        
    Returns:
        Filtered DataFrame
    """
    if df.empty:
        return df
        
    filtered = df.copy()
    
    # Filter by bit score
    if min_bitscore is not None and 'bitscore' in filtered.columns:
        filtered = filtered[filtered['bitscore'] >= min_bitscore]
        
    # Filter by e-value
    if max_evalue is not None and 'evalue' in filtered.columns:
        filtered = filtered[filtered['evalue'] <= max_evalue]
        
    return filtered