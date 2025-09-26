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

    def run_pling(self, input_fasta: str, output_dir: str, pling_args: str = None) -> None:
        """
        Run Pling
        """

        pling_cmd = ['pling', input_fasta, output_dir, 'align', '--cores', str(self.config.THREADS)]

        if pling_args:
            pling_cmd += pling_args.strip('"').strip("'").split()
        try:
            subprocess.check_call(pling_cmd)
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(f"Pling failed: {str(e)}")
        self.output_dir = output_dir
    
    def run_plasnet(self, new_dcj_dists: str, dcj_threshold: float = None):
        comm_pkl = Path(self.output_dir) / 'containment' / 'containment_communities' / 'objects' / 'communities.pkl'
        cmd = ['plasnet', 'type', '-d', str(dcj_threshold), 
               '--output-type', 'both', 
               str(comm_pkl),  new_dcj_dists, 
               str(Path(self.output_dir) / 'pling_dcj_length_subcommunities')]
        try:
            subprocess.check_call(cmd)
            return Path(self.output_dir) / 'pling_dcj_length_subcommunities' 
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(f"Plasnet failed: {str(e)}")
    
    
    def get_dcj_dist_by_length(self, row: pd.Series, threshold: float =.05, min_len: float = 120000, dcj_threshold: int = None):
        min_length = min(row['length1'], row['length2'])
        dcj_dist = row['distance_y'] / (min_length / 1000)
        if min_length < min_len and row['distance_y'] <= dcj_threshold and dcj_dist > threshold:
            dcj_dist = threshold - .01
        if min_length < min_len and row['distance_y'] > dcj_threshold and dcj_dist < threshold:
            dcj_dist = threshold + .01
        return dcj_dist

    def edit_dcj(self, plasmid_info: pd.DataFrame, threshold: float = None, plasmid_length: float = None, min_dcj = None):
        dcj_threshold = threshold or self.config.DCJ_THRESHOLD
        dcj_num = min_dcj or self.config.DCJ_THRESHOLD_NUM
        min_len = plasmid_length or self.config.DCJ_THRESHOLD_LEN
        new_dcj_dists = Path(self.output_dir) / 'pling_dcj_dists_by_length.tsv'
        pling_dcj_dists = pd.read_csv(Path(self.output_dir) / 'all_plasmids_distances.tsv', sep='\t')
        pling_dcj_dists['distance_y'] = pling_dcj_dists['distance']
        pling_dcj_dists['length1'] = pling_dcj_dists['plasmid_1'].map(lambda s : plasmid_info.loc[s, 'size'])
        pling_dcj_dists['length2'] = pling_dcj_dists['plasmid_2'].map(lambda s: plasmid_info.loc[s, 'size'])
        pling_dcj_dists['distance'] = pling_dcj_dists.apply(lambda r: self.get_dcj_dist_by_length(r, threshold=dcj_threshold, min_len=min_len, dcj_threshold=dcj_num), axis=1)
    
        pling_dcj_dists.to_csv(new_dcj_dists, columns=['plasmid_1', 'plasmid_2', 'distance'], sep='\t')

        new_dcj_types = self.run_plasnet(new_dcj_dists, dcj_threshold)
        return new_dcj_types    
    

