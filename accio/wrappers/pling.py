import pandas as pd
import subprocess
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple


from ..config import AnalysisConfig
from .external_tools import ExternalToolError


class PlingRunner:
    """Wrapper for Pling operations"""

    def __init__(self, config: Optional[AnalysisConfig] = None):
        self.config = config or AnalysisConfig()

    def run_pling(self, input_fasta: str, output_dir: str, pling_args: str = None,
                  previous_pling_dir: str = None) -> None:
        """
        Run Pling
        """
        
        pling_command = ['add'] if previous_pling_dir else ['cluster', 'align']
        pling_cmd = ['pling'] + pling_command + [ input_fasta, output_dir, '--cores', str(self.config.THREADS)]

        if pling_args:
            pling_cmd += pling_args.strip('"').strip("'").split()
        if previous_pling_dir:
            pling_cmd += [previous_pling_dir]

        try:
            subprocess.check_call(pling_cmd)
        except subprocess.CalledProcessError as e:
            if previous_pling_dir:
                self.run_pling(input_fasta, output_dir, pling_args=pling_args)
            else:
                raise ExternalToolError(f"Pling failed: {str(e)}")
        self.output_dir = output_dir
    
    def run_plasnet(self, new_dcj_dists: str, dcj_threshold: float = None, previous_pling_dir: str = None):
        comm_pkl = Path(self.output_dir) / 'containment' / 'containment_communities' / 'objects' / 'communities.pkl'
        cmd = ['plasnet', 'type', '-d', str(dcj_threshold), 
               '--output-type', 'both', 
               str(comm_pkl),  new_dcj_dists, 
               str(Path(self.output_dir) / 'pling_dcj_length_subcommunities')]
        if previous_pling_dir:
            prev_pling_dist = Path(previous_pling_dir) / 'pling_dcj_length_subcommunities' / 'objects' / 'typing.tsv'
            cmd += ['--prev_typing', str(prev_pling_dist), '--reclustering_method', 'biased']
        try:
            subprocess.check_call(cmd)
            return Path(self.output_dir) / 'pling_dcj_length_subcommunities' 
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(f"Plasnet failed: {str(e)}")
    
    
    def get_dcj_dist_by_length(self, row: pd.Series, threshold: float =.05, min_len: float = 120000, dcj_threshold: int = None):
        min_length = min(row['length1'], row['length2'])
        dcj_dist = row['distance_y'] / (min_length / 1000)
        if min_length <= min_len and row['distance_y'] <= dcj_threshold and dcj_dist > threshold:
            dcj_dist = threshold - .01
        if min_length <= min_len and row['distance_y'] > dcj_threshold and dcj_dist < threshold:
            dcj_dist = threshold + .01
        return dcj_dist

    def edit_dcj(self, plasmid_info: pd.DataFrame, threshold: float = None, plasmid_length: float = None, min_dcj = None, 
                 previous_pling_dir: str = None):
        dcj_threshold = threshold or self.config.DCJ_THRESHOLD
        dcj_num = min_dcj or self.config.DCJ_THRESHOLD_NUM
        min_len = plasmid_length or self.config.DCJ_THRESHOLD_LEN
        new_dcj_dists = Path(self.output_dir) / 'pling_dcj_dists_by_length.tsv'
        plasmid_info = plasmid_info.set_index('plasmid_id') if 'plasmid_id' in plasmid_info.columns else plasmid_info
        pling_dcj_dists = pd.read_csv(Path(self.output_dir) / 'all_plasmids_distances.tsv', sep='\t')

        # Normalize column names in distance dataframe if needed
        col_map = {}
        for c in pling_dcj_dists.columns:
            c_clean = str(c).strip().lower().replace('-', '_').replace(' ', '_')
            if c_clean in ['plasmid1', 'query', 'p1', 'plasmid_1']:
                col_map[c] = 'plasmid_1'
            elif c_clean in ['plasmid2', 'ref', 'reference', 'p2', 'plasmid_2']:
                col_map[c] = 'plasmid_2'
            elif c_clean in ['dist', 'distance']:
                col_map[c] = 'distance'
        if col_map:
            pling_dcj_dists = pling_dcj_dists.rename(columns=col_map)

        # Detect plasmid size/length column name
        size_col = None
        for col_name in ['size', 'length', 'total_length', 'len']:
            if col_name in plasmid_info.columns:
                size_col = col_name
                break

        def _get_len(pid):
            if size_col and pid in plasmid_info.index:
                val = plasmid_info.loc[pid, size_col]
                if isinstance(val, pd.Series):
                    val = val.iloc[0]
                try:
                    return float(val) if pd.notna(val) else min_len
                except (ValueError, TypeError):
                    return min_len
            return min_len

        pling_dcj_dists['distance_y'] = pling_dcj_dists['distance']
        pling_dcj_dists['length1'] = pling_dcj_dists['plasmid_1'].map(_get_len)
        pling_dcj_dists['length2'] = pling_dcj_dists['plasmid_2'].map(_get_len)
        pling_dcj_dists = pling_dcj_dists.dropna(how='any')
        pling_dcj_dists['distance'] = pling_dcj_dists.apply(lambda r: self.get_dcj_dist_by_length(r, threshold=dcj_threshold, min_len=min_len, dcj_threshold=dcj_num), axis=1)
    
        pling_dcj_dists.to_csv(new_dcj_dists, columns=['plasmid_1', 'plasmid_2', 'distance'], sep='\t')

        new_dcj_types = self.run_plasnet(new_dcj_dists, dcj_threshold, previous_pling_dir=previous_pling_dir)
        return new_dcj_types    
    

