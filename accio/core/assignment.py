"""Plasmid assignment algorithm."""

import pandas as pd
from typing import Dict, List, Tuple, Any, Optional
import logging

from .contig import Contig
from .plasmid import Plasmid
from ..config import AnalysisConfig


def assign_plasmids(plasmids: Dict[str, Plasmid], 
                   contigs: Dict[str, Contig],
                   config: Optional[AnalysisConfig] = None) -> Tuple[List[Plasmid], Dict[str, str]]:
    """
    Assign contigs to plasmids using a greedy algorithm.
    
    This algorithm prioritizes plasmids by their similarity scores and assigns
    contigs while avoiding conflicts and ensuring biological consistency.
    
    Args:
        plasmids: Dictionary of plasmid ID to Plasmid objects
        contigs: Dictionary of contig ID to Contig objects  
        config: Analysis configuration
        
    Returns:
        Tuple of (assigned_plasmids, contig_assignments)
        where contig_assignments maps contig_id -> plasmid_id
    """
    if config is None:
        config = AnalysisConfig()
        
    logger = logging.getLogger(__name__)
    
    # Get database info for cluster tracking

    try:
        all_pling_counts = _get_pling_type_counts(plasmids) if config.USE_PLING_COMMUNITY_COUNTS else pd.Series()
    except Exception:
        all_pling_counts = pd.Series()
        
    # Sort plasmids by score (highest first)
    plasmids_to_process = sorted(
        plasmids.keys(),
        key=lambda p: -plasmids[p].score
    )
    
    # Track assignments
    plasmids_assigned = []
    contigs_assigned = {}
    clusters_assigned = {
        'pling_type': [],
        'mob_cluster': []
    }
    rep_types_assigned = []
    contigs_available = set(contigs.keys())
    
    logger.info(f"Starting assignment of {len(plasmids_to_process)} plasmids to {len(contigs_available)} contigs")
    
    while plasmids_to_process and contigs_available:
        # Select best candidate plasmid
        best_plasmid_id = _select_best_plasmid(
            plasmids_to_process, plasmids, all_pling_counts
        )
        
        if best_plasmid_id is None:
            logger.warning("No suitable plasmid found, stopping assignment")
            break
            
        plasmid = plasmids[best_plasmid_id]
        logger.debug(f"Processing plasmid {plasmid.id} (score: {plasmid.score:.3f})")
        
        # Get contigs that should be assigned to this plasmid
        candidate_contigs = plasmid.get_contigs_from_blast()
        
        # Remove contigs that are already assigned (except repetitive ones)
        contigs_to_remove = []
        for contig_id in plasmid.contigs:
            if (contig_id in contigs_assigned and 
                contigs[contig_id].contig_type != 'repetitive'):
                contigs_to_remove.append(contig_id)
            elif contig_id not in candidate_contigs:
                contigs_to_remove.append(contig_id)
            elif contig_id not in contigs:
                contigs_to_remove.append(contig_id)
                
        # Update plasmid by removing unavailable contigs
        if contigs_to_remove:
            plasmid.remove_contigs(contigs_to_remove)
            
        # Check if plasmid still meets assignment criteria
        if not _meets_assignment_criteria(plasmid, config):
            plasmids_to_process.remove(best_plasmid_id)
            logger.debug(f"Plasmid {plasmid.id} does not meet assignment criteria {plasmid.scores}")
            logger.debug(plasmid.scores)
            continue
            
        # Check for conflicts with already assigned clusters
        if _has_cluster_conflict(plasmid, clusters_assigned, rep_types_assigned):
            plasmids_to_process.remove(best_plasmid_id)
            logger.debug(f"Plasmid {plasmid.id} conflicts with assigned clusters")
            continue
            
        # Assign this plasmid
        plasmids_assigned.append(plasmid)
        clusters_assigned['pling_type'].append(plasmid.other.get('pling_type'))
        clusters_assigned['mob_cluster'].append(plasmid.cluster_id)
        
        if plasmid.rep_type:
            rep_types_assigned.append(plasmid.rep_type)
            
        logger.info(f"Assigned plasmid {plasmid.id} with {len(plasmid.contigs)} contigs")
        

        # Mark contigs as assigned
        for contig_id in plasmid.contigs:
            contigs[contig_id].assign(plasmid.id, plasmid.cluster_id)
            contigs_assigned[contig_id] = plasmid.id
            
            # Remove non-repetitive contigs from available pool
            if contigs[contig_id].contig_type != 'repetitive':
                contigs_available.discard(contig_id)
        plasmids_to_process.remove(plasmid.id)
        
        # Update remaining plasmids by removing assigned contigs
        _update_remaining_plasmids(
            plasmids_to_process, plasmids, contigs_assigned, contigs
        )
        
        # Filter remaining plasmids by updated criteria
        plasmids_to_process = _filter_remaining_plasmids(
            plasmids_to_process, plasmids, clusters_assigned, 
            rep_types_assigned, config
        )
        
        logger.debug(f"Remaining plasmids to process: {len(plasmids_to_process)}")
        
    logger.info(f"Assignment complete: {len(plasmids_assigned)} plasmids assigned")
    return plasmids_assigned, contigs_assigned


def _get_pling_type_counts(plasmids: Dict[str, Plasmid]) -> pd.Series:
    """Get counts of pling types in the plasmid database."""
    pling_types = [p.other.get('pling_type') for p in plasmids.values() 
                  if p.other.get('pling_type')]
    return pd.Series(pling_types).value_counts()


def _select_best_plasmid(plasmids_to_process: List[str], 
                        plasmids: Dict[str, Plasmid],
                        all_pling_counts: pd.Series) -> Optional[str]:
    """
    Select the best plasmid candidate for assignment.
    
    Prioritizes by pling type frequency and then by score.
    """
    if not plasmids_to_process:
        return None
        
    # Filter to high-scoring plasmids
    high_score_plasmids = [
        p for p in plasmids_to_process 
        if plasmids[p].score >= 0.80
    ]
    
    if not high_score_plasmids:
        # Fall back to all plasmids if none are high-scoring
        high_score_plasmids = plasmids_to_process
        
    # Try to select by pling type preference
    if not all_pling_counts.empty:
        try:
            # Get pling type counts for current candidates
            current_pling_counts = {}
            for pid in high_score_plasmids:
                pling_type = plasmids[pid].other.get('pling_type')
                if pling_type:
                    current_pling_counts[pling_type] = current_pling_counts.get(pling_type, 0) + 1
                    
            if current_pling_counts:
                # Calculate relative frequencies
                current_series = pd.Series(current_pling_counts)
                pling_ratios = current_series / all_pling_counts.reindex(current_series.index, fill_value=1)
                
                # Select plasmid with most frequent pling type
                best_pling_type = pling_ratios.idxmax()
                candidates = [
                    p for p in high_score_plasmids
                    if plasmids[p].other.get('pling_type') == best_pling_type
                ]
                
                if candidates:
                    return candidates[0]  # Already sorted by score
                    
        except Exception:
            # Fall back to score-based selection
            pass
            
    # Default: return highest scoring plasmid
    return high_score_plasmids[0]


def _meets_assignment_criteria(plasmid: Plasmid, config: AnalysisConfig) -> bool:
    """Check if plasmid meets criteria for assignment."""
    pass_mash_score = plasmid.mash_score >= config.MIN_MASH_SCORE
    return (
        plasmid.scores['blast_similarity'] >= config.MIN_BLAST_SIMILARITY and
        pass_mash_score
    )


def _has_cluster_conflict(plasmid: Plasmid, 
                         clusters_assigned: Dict[str, List],
                         rep_types_assigned: List[str]) -> bool:
    """Check if plasmid conflicts with already assigned clusters."""
    
    # Check pling type conflict
    if plasmid.other.get('pling_type') in clusters_assigned['pling_type']:
        return True
        
    # Check MOB cluster conflict  
    if plasmid.cluster_id in clusters_assigned['mob_cluster']:
        return True
        
    # Check replicon type conflict
    if plasmid.rep_type and plasmid.rep_type in rep_types_assigned:
        return True
        
    return False


def _update_remaining_plasmids(plasmids_to_process: List[str],
                              plasmids: Dict[str, Plasmid],
                              contigs_assigned: Dict[str, str],
                              contigs: Dict[str, Contig]) -> None:
    """Update remaining plasmids by removing assigned contigs."""
    
    for plasmid_id in plasmids_to_process:
        plasmid = plasmids[plasmid_id]
        
        # Find contigs to remove
        contigs_to_remove = []
        for contig_id in plasmid.contigs:
            if (contig_id in contigs_assigned and 
                contigs[contig_id].contig_type != 'repetitive'):
                contigs_to_remove.append(contig_id)
                
        # Remove conflicting contigs
        if contigs_to_remove:
            plasmid.remove_contigs(contigs_to_remove)


def _filter_remaining_plasmids(plasmids_to_process: List[str],
                              plasmids: Dict[str, Plasmid],
                              clusters_assigned: Dict[str, List],
                              rep_types_assigned: List[str],
                              config: AnalysisConfig) -> List[str]:
    """Filter remaining plasmids by assignment criteria."""
    
    filtered_plasmids = []
    
    for plasmid_id in plasmids_to_process:
        plasmid = plasmids[plasmid_id]
        
        # Skip if doesn't meet basic criteria
        if not _meets_assignment_criteria(plasmid, config):
            continue
            
        # Skip if conflicts with assigned clusters
        if _has_cluster_conflict(plasmid, clusters_assigned, rep_types_assigned):
            continue
            
        filtered_plasmids.append(plasmid_id)
        
    # Re-sort by score
    return sorted(filtered_plasmids, key=lambda p: -plasmids[p].score)


def calculate_assignment_confidence(plasmid: Plasmid) -> float:
    """
    Calculate confidence score for plasmid assignment.
    
    Args:
        plasmid: Plasmid object
        
    Returns:
        Confidence score between 0 and 1
    """
    factors = []
    
    # High BLAST similarity increases confidence
    if plasmid.scores['blast_similarity'] > 90:
        factors.append(0.3)
    elif plasmid.scores['blast_similarity'] > 80:
        factors.append(0.2)
    else:
        factors.append(0.1)
        
    # High Mash score increases confidence
    if plasmid.mash_score > 0.99:
        factors.append(0.25)
    elif plasmid.mash_score > 0.95:
        factors.append(0.15)
    else:
        factors.append(0.05)
        
    # Good synteny increases confidence
    if plasmid.scores['synteny_score'] > 0.8:
        factors.append(0.2)
    elif plasmid.scores['synteny_score'] > 0.5:
        factors.append(0.1)
    else:
        factors.append(0.05)
        
    # Length similarity increases confidence
    if plasmid.scores['length_score'] > 0.9:
        factors.append(0.15)
    elif plasmid.scores['length_score'] > 0.7:
        factors.append(0.1)
    else:
        factors.append(0.05)
        
    # Read data increases confidence if available
    if plasmid.read_cov_data is not None:
        if plasmid.scores['read_gaps_score'] > 0.8:
            factors.append(0.1)
        else:
            factors.append(0.05)
            
    return min(sum(factors), 1.0)


def resolve_assignment_conflicts(assignments: Dict[str, str],
                                plasmids: Dict[str, Plasmid],
                                contigs: Dict[str, Contig]) -> Dict[str, str]:
    """
    Resolve conflicts in plasmid assignments.
    
    Args:
        assignments: Dictionary of contig_id -> plasmid_id
        plasmids: Dictionary of plasmid objects
        contigs: Dictionary of contig objects
        
    Returns:
        Resolved assignments dictionary
    """
    logger = logging.getLogger(__name__)
    
    # Find contigs assigned to multiple plasmids
    conflicts = {}
    for contig_id, plasmid_id in assignments.items():
        if contig_id in conflicts:
            conflicts[contig_id].append(plasmid_id)
        else:
            conflicts[contig_id] = [plasmid_id]
            
    conflicts = {k: v for k, v in conflicts.items() if len(v) > 1}
    
    if not conflicts:
        return assignments
        
    logger.info(f"Resolving {len(conflicts)} assignment conflicts")
    
    resolved_assignments = assignments.copy()
    
    for contig_id, conflicting_plasmids in conflicts.items():
        # Choose plasmid with highest confidence
        best_plasmid = None
        best_confidence = 0
        
        for plasmid_id in conflicting_plasmids:
            confidence = calculate_assignment_confidence(plasmids[plasmid_id])
            if confidence > best_confidence:
                best_confidence = confidence
                best_plasmid = plasmid_id
                
        if best_plasmid:
            resolved_assignments[contig_id] = best_plasmid
            logger.debug(f"Resolved conflict for {contig_id}: assigned to {best_plasmid}")
            
    return resolved_assignments


def validate_assignments(assignments: Dict[str, str],
                        plasmids: Dict[str, Plasmid],
                        contigs: Dict[str, Contig]) -> List[str]:
    """
    Validate plasmid assignments and return warnings.
    
    Args:
        assignments: Dictionary of contig_id -> plasmid_id
        plasmids: Dictionary of plasmid objects  
        contigs: Dictionary of contig objects
        
    Returns:
        List of validation warning messages
    """
    warnings = []
    
    # Check for size mismatches
    plasmid_sizes = {}
    for contig_id, plasmid_id in assignments.items():
        if plasmid_id not in plasmid_sizes:
            plasmid_sizes[plasmid_id] = {
                'expected': plasmids[plasmid_id].length,
                'actual': 0,
                'contigs': []
            }
        plasmid_sizes[plasmid_id]['actual'] += contigs[contig_id].length
        plasmid_sizes[plasmid_id]['contigs'].append(contig_id)
        
    for plasmid_id, size_info in plasmid_sizes.items():
        size_diff = abs(size_info['expected'] - size_info['actual'])
        size_ratio = size_diff / size_info['expected']
        
        if size_ratio > 0.5:  # More than 50% size difference
            warnings.append(
                f"Large size mismatch for {plasmid_id}: "
                f"expected {size_info['expected']}bp, "
                f"got {size_info['actual']}bp from contigs {size_info['contigs']}"
            )
            
    # Check for missing essential genes (if data available)
    # This would require additional gene annotation data
    
    # Check for unusual copy numbers
    for contig_id, plasmid_id in assignments.items():
        contig = contigs[contig_id]
        if hasattr(contig, 'copy_num') and contig.copy_num:
            if contig.copy_num > 5:
                warnings.append(
                    f"High copy number ({contig.copy_num:.1f}x) for contig {contig_id} "
                    f"assigned to {plasmid_id}"
                )
            elif contig.copy_num < 0.5:
                warnings.append(
                    f"Low copy number ({contig.copy_num:.1f}x) for contig {contig_id} "
                    f"assigned to {plasmid_id}"
                )
                
    return warnings