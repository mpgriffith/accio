"""Database builder for Accio plasmid databases."""

import os
import logging
import hashlib
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import scipy.cluster.hierarchy as sch

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..wrappers import (BLASTRunner, MashRunner, MOBRunner, PlingRunner, SkaniRunner)
from ..config import MOB_TYPER_COLS


class DatabaseBuilder: 
    """Build Accio plasmid databases from reference sequences."""
    
    def __init__(self, threads: int = 4, min_length: int = 10000, max_length: int = 1000000):
        """
        Initialize database builder.
        
        Args:
            threads: Number of threads for parallel processing
            min_length: Minimum sequence length to include
            max_length: Maximum sequence length to include
        """
        self.threads = threads
        self.min_length = min_length
        self.max_length = max_length
        self.logger = logging.getLogger(__name__)
        
        # Initialize tool runners
        self.blast_runner = BLASTRunner()
        self.mash_runner = MashRunner()
        self.mob_runner = MOBRunner()
        self.skani_runner = SkaniRunner()
        self.pling_runner = PlingRunner()

    def build_database(self, 
                      input_files: List[str],
                      output_dir: str,
                      metadata_file: Optional[str] = None,
                      update_existing: bool = False,
                      check_circular: bool = True,
                      no_cluster: bool=False,
                      deduplicate: bool = False,
                      validate_sequences: bool = False,
                      cluster_threshold: float = 0.95,
                      use_all: bool = False) -> Dict[str, Any]:
        """
        Build plasmid database from input sequences.
        
        Args:
            input_files: List of FASTA files containing plasmid sequences
            output_dir: Output directory for database files
            metadata_file: Optional CSV file with plasmid metadata
            update_existing: Whether to update existing database
            check_circular: Check for circular sequences
            deduplicate: Remove duplicate sequences
            validate_sequences: Validate sequence quality
            cluster_threshold: Mash distance threshold for clustering
            
        Returns:
            Dictionary with build statistics
        """
        self.logger.info("Starting database construction...")
        
        check_circular = True if use_all == False else check_circular

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        (output_path / 'fastas').mkdir(exist_ok=True)
        (output_path / 'mob').mkdir(exist_ok=True)

        # Load existing data if updating
        existing_sequences = {}
        existing_metadata = pd.DataFrame()
        if update_existing and os.path.exists(output_dir):
            existing_sequences, existing_metadata = self._load_existing_database(output_dir)
            self.logger.info(f"Loaded {len(existing_sequences)} existing sequences")
        
        # Process input sequences
        self.logger.info("Processing input sequences...")
        sequences, stats = self._process_input_sequences(
            input_files, validate_sequences, check_circular, use_all
        )
        self.logger.info(f"Processed {stats['total_input']} input sequences")
        self.logger.info(f"Filtered to {len(sequences)} valid sequences")
        
        # Merge with existing sequences
        if existing_sequences:
            sequences.update(existing_sequences)
            self.logger.info(f"Combined with existing: {len(sequences)} total unique sequences")

        # Remove duplicates if requested
        if deduplicate:
            sequences = self._deduplicate_sequences(sequences)
            self.logger.info(f"After deduplication: {len(sequences)} sequences")
            stats['duplicates_removed'] = stats['total_input'] - len(sequences)
        
        if not update_existing:
            self.download_plasmidfinder(output_dir)

        # Type sequences and filter for plasmids
        plasmids, mob_df, plasmidfinder_results, mob_typer_file = self._type_and_filter_sequences(
            sequences, output_path, use_all
        )

        # Cluster plasmids to find representatives
        self.logger.info("Clustering sequences...")
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta") as temp_plasmids_file:
            SeqIO.write(plasmids.values(), temp_plasmids_file, 'fasta')
            plasmids_fasta_path = temp_plasmids_file.name

        rep_plasmids, cluster_info = self.cluster_plasmids_skani(
            plasmids_fasta_path, output_path, sim_threshold=cluster_threshold, skip_cluster=no_cluster
        )
        
        # Process and merge metadata from all sources
        metadata_df = self._process_metadata(metadata_file, existing_metadata, sequences.keys())
        metadata_df = self._merge_cluster_info(metadata_df, cluster_info, mob_df, plasmidfinder_results)
        self.logger.info("\n" + metadata_df.head().to_string())
        metadata_df.to_csv(output_dir + '/plasmid_db_before_pling.csv')
        metadata_df = self._create_plasmid_communities(rep_plasmids, metadata_df, output_dir)

        metadata_df = self._run_mob_cluster(rep_plasmids, metadata_df, mob_typer_file, output_dir)
        
        # Write final database files
        self._write_database_files(rep_plasmids, metadata_df, output_dir)

        # Build BLAST and Mash databases
        self._build_blast_database(output_dir)
        self._build_mash_database(output_dir)
        
        # Create summary statistics
        final_stats = {
            'total_sequences': stats['total_input'],
            'final_sequences': len(rep_plasmids),
            'clusters_created': cluster_info['cluster'].nunique() if not cluster_info.empty else 0,
            'duplicates_removed': stats.get('duplicates_removed', 0),
            'circular_sequences': stats.get('circular_count', 0),
            'database_files': self._get_database_files(output_dir)
        }
        
        # Clean up temporary files
        os.remove(plasmids_fasta_path)

        self.logger.info("Database construction completed!")
        return final_stats
    
    def _load_existing_database(self, db_dir: str) -> Tuple[Dict[str, SeqRecord], pd.DataFrame]:
        """Load existing database sequences and metadata."""
        sequences = {}
        metadata = pd.DataFrame()
        
        # Load sequences
        fasta_path = os.path.join(db_dir, "plasmidDB.fasta")
        if os.path.exists(fasta_path):
            sequences = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
            
        # Load metadata
        metadata_path = os.path.join(db_dir, "plasmidDB_info.csv")
        if os.path.exists(metadata_path):
            metadata = pd.read_csv(metadata_path)
        
        return sequences, metadata
    
    def download_plasmidfinder(self, output_dir):
        plasmidfinder_dir = Path(output_dir) / "plasmidfinder_db"
        plasmidfinder_dir.mkdir(exist_ok=True)
        plasmidfinder_fasta = Path(output_dir) / 'PlasmidFinder.fasta'
        plasmidfinder_db_path = Path(output_dir) / 'PlasmidFinder'
        self.logger.info("Downloading PlasmidFinder database...")
        try:
            # Clone the PlasmidFinder repository
            clone_cmd = ['git', 'clone', 'https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git', str(plasmidfinder_dir)]
            clone_cmd = subprocess.run(clone_cmd, check=True, capture_output=True)
            fastas = list(plasmidfinder_dir.glob('*.fsa'))
            with plasmidfinder_fasta.open("w") as outfile:
                for file_path in fastas:
                    with file_path.open("r") as infile:
                        for line in infile:
                            outfile.write(line)
            self.blast_runner.run_makeblastdb(str(plasmidfinder_fasta), str(plasmidfinder_db_path))
            # Clean up
            shutil.rmtree(plasmidfinder_dir)
            self.logger.info("PlasmidFinder database downloaded and set up successfully.")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to download or set up PlasmidFinder database: {e.stderr.decode()}")
            raise

    def _type_and_filter_sequences(self, sequences: Dict[str, SeqRecord], output_path: Path, use_all: bool):
        """Run typing tools and filter sequences to identify plasmids."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta") as temp_fasta:
            SeqIO.write(sequences.values(), temp_fasta, 'fasta')
            all_seqs_file = temp_fasta.name

        # Run Mob Typer on all sequences
        mob_df, mob_typer_file = self._run_mob_typer(all_seqs_file, output_path)
        
        # Run PlasmidFinder on all sequences
        plasmidfinder_results, _ = self._run_plasmidfinder(all_seqs_file, str(output_path))

        # Filter for sequences that are likely plasmids
        plasmids = self._filter_plasmids(sequences, mob_df, plasmidfinder_results, use_all)
        return plasmids, mob_df, plasmidfinder_results, mob_typer_file
    
    def _run_plasmidfinder(self, sequence_file: str, output_dir: str) -> pd.DataFrame:
        """Run PlasmidFinder on sequences."""
        self.logger.info("Running PlasmidFinder...")
        
        results_file = str(Path(output_dir) / "plasmidfinder_results.tsv")
        plasmidfinder_results = self.blast_runner.run_blastn(
            query=sequence_file, 
            database=str(Path(output_dir) / 'PlasmidFinder'),
            output_file=results_file
        )
        
        return plasmidfinder_results, results_file
    
    def _run_mob_typer(self, sequence_file: str, output_dir: Path):
        """Run Mob Typer on sequences."""
        mob_runner = MOBRunner()
        results_file =  str(output_dir / 'mob_typer_results.tsv')
        results_file = mob_runner.run_mob_typer(sequence_file, str(output_dir / 'mob_typer_results.tsv'), multi=True)
        mob_df = pd.read_csv(results_file, sep='\t')
        return mob_df, results_file

    def _filter_plasmids(self, sequences: Dict[str, SeqRecord], mob_results: pd.DataFrame,
                       plasmidfinder_results: pd.DataFrame,
                       use_all: bool=False) -> Dict[str, SeqRecord]:
        """Filter plasmids based on MOB-typer and PlasmidFinder results."""
        keep_ids = set()
        
        for p_id, p_record in sequences.items():
            mob_plasmid = mob_results[mob_results['sample_id'] == p_id]
            pf_plasmid = plasmidfinder_results[plasmidfinder_results['qseqid'] == p_id]
            
            # Keep if `use_all` is true, or if it has any PlasmidFinder hits
            if use_all or not pf_plasmid.empty:
                keep_ids.add(p_id)
            # Keep if MOB-typer found a replicon type
            elif not mob_plasmid.empty and mob_plasmid.iloc[0].get('rep_type(s)', '-') != '-':
                keep_ids.add(p_id)
                    
        self.logger.info(f"Filtered to {len(keep_ids)} plasmids with replicon evidence.")
        return {p_id: sequences[p_id] for p_id in keep_ids}
        
    def _process_input_sequences(self, input_files: List[str], 
                               validate: bool, check_circular: bool,
                               use_all: bool = False) -> Tuple[Dict[str, SeqRecord], Dict[str, Any]]:
        """Process and filter input sequences."""
        sequences = {}
        stats = {
            'total_input': 0,
            'filtered_length': 0,
            'filtered_quality': 0,
            'circular_count': 0
        }
        
        for file_path in input_files:
            self.logger.info(f"Processing {file_path}")
            
            try:
                for record in SeqIO.parse(file_path, "fasta"):
                    stats['total_input'] += 1
                    

                    # Filter by length
                    if  not (self.min_length <= len(record.seq) <= self.max_length):
                        stats['filtered_length'] += 1
                        continue
                    
                    # Validate sequence quality
                    if validate and not self._validate_sequence(record):
                        stats['filtered_quality'] += 1
                        continue
                    
                    # Check for circularity
                    if not use_all and check_circular:
                        is_circular = self._check_circularity(record)
                        if is_circular:
                            record.annotations['topology'] = 'circular'
                            stats['circular_count'] += 1
                    
                    # Clean sequence ID
                    clean_id = self._clean_sequence_id(record.id)
                    record.id = clean_id
                    record.name = clean_id
                    
                    sequences[clean_id] = record
                    
            except Exception as e:
                self.logger.warning(f"Error processing {file_path}: {e}")
                continue
        
        return sequences, stats
    
    def _validate_sequence(self, record: SeqRecord) -> bool:
        """Validate sequence quality."""
        seq_str = str(record.seq).upper()
        
        # Check for excessive N's
        n_content = seq_str.count('N') / len(seq_str)
        if n_content > 0.05:  # More than 5% N's
            return False
            
        # Check for valid nucleotides
        valid_bases = set('ATCGRYSWKMBDHVN')
        invalid_bases = set(seq_str) - valid_bases
        if invalid_bases:
            return False
            
        return True
    
    def _check_circularity(self, record: SeqRecord) -> bool:
        """Check if a sequence is likely circular based on circular in description or terminal overlaps."""
        seq_str = str(record.seq)
        seq_len = len(seq_str)
        if 'circular=true' in record.description.lower():
            return True
        # Check for terminal overlaps (common in circular assemblies)
        overlap_sizes = [20, 50, 100, 200]
        
        for overlap_size in overlap_sizes:
            if overlap_size > seq_len // 4:  # Don't check overlaps larger than 1/4 of sequence
                continue
                
            start_seq = seq_str[:overlap_size]
            end_seq = seq_str[-overlap_size:]
            
            # Check for exact match
            if start_seq == end_seq:
                return True
                
            # Check for reverse complement match (indicates circular with wrong orientation)
            from Bio.Seq import Seq
            start_rc = str(Seq(start_seq).reverse_complement())
            if start_rc == end_seq:
                return True
        
        return False
    
    def _clean_sequence_id(self, seq_id: str) -> str:
        """Clean and standardize sequence IDs."""
        # Remove common prefixes and problematic characters
        clean_id = seq_id.split()[0]  # Take first part before whitespace
        clean_id = clean_id.replace('|', '_')
        clean_id = clean_id.replace(':', '_')
        clean_id = clean_id.replace(';', '_')
        clean_id = clean_id.replace(',', '_')
        
        return clean_id
    
    def _deduplicate_sequences(self, sequences: Dict[str, SeqRecord]) -> Dict[str, SeqRecord]:
        """Remove duplicate sequences based on sequence hash."""
        unique_sequences = {}
        seen_hashes = set()
        
        for seq_id, record in sequences.items():
            seq_hash = hashlib.md5(str(record.seq).encode()).hexdigest()
            
            if seq_hash not in seen_hashes:
                unique_sequences[seq_id] = record
                seen_hashes.add(seq_hash)
            else:
                self.logger.debug(f"Removing duplicate sequence: {seq_id}")
        
        return unique_sequences
    
    def _process_metadata(self, metadata_file: Optional[str], 
                         existing_metadata: pd.DataFrame,
                         sequence_ids: List[str]) -> pd.DataFrame:
        """Process and validate metadata."""
        metadata_df = pd.DataFrame()
        
        # Load new metadata if provided
        if metadata_file and os.path.exists(metadata_file):
            try:
                metadata_df = pd.read_csv(metadata_file)
                self.logger.info(f"Loaded metadata for {len(metadata_df)} sequences")
            except Exception as e:
                self.logger.warning(f"Error loading metadata: {e}")
                metadata_df = pd.DataFrame()
        
        # Combine with existing metadata
        if not existing_metadata.empty:
            metadata_df = pd.concat([existing_metadata, metadata_df], ignore_index=True)
            metadata_df = metadata_df.drop_duplicates(subset=['plasmid_id'], keep='last')
        
        # Create default metadata for sequences without metadata
        all_seq_ids = set(sequence_ids)
        existing_ids = set(metadata_df.get('plasmid_id', []))
        missing_ids = all_seq_ids - existing_ids
        
        if missing_ids:
            default_metadata = pd.DataFrame({
                'plasmid_id': list(missing_ids),
                'cluster_id': ['unknown'] * len(missing_ids),
                'rep_types': [None] * len(missing_ids),
                'length': [0] * len(missing_ids),  # Will be updated later
                'source': ['user_provided'] * len(missing_ids)
            })
            metadata_df = pd.concat([metadata_df, default_metadata], ignore_index=True)
        
        # Ensure required columns exist
        required_columns = ['plasmid_id', 'cluster_id', 'rep_type', 'length']
        for col in required_columns:
            if col not in metadata_df.columns:
                metadata_df[col] = 'unknown' if col != 'length' else 0
        
        return metadata_df
    
    def _cluster_sequences(self, sequences: Dict[str, SeqRecord], 
                          threshold: float = 0.95) -> Dict[str, str]:
        #TODO: cluster using skani instead
        """Cluster sequences using Mash distances."""
        if len(sequences) < 2:
            # Single sequence or no sequences
            return {list(sequences.keys())[0]: 'cluster_1'} if sequences else {}
        
        cluster_info = {}
        
        try:
            # Create temporary files
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
                SeqIO.write(sequences.values(), temp_fasta, 'fasta')
                temp_fasta_path = temp_fasta.name
            
            # Create Mash sketch
            sketch_path = self.mash_runner.run_sketch(temp_fasta_path, temp_fasta_path)
            
            # Calculate distances
            mash_dist_cmd = ['mash', 'dist', sketch_path, sketch_path]
            result = subprocess.run(mash_dist_cmd, check=True, capture_output=True, text=True)
            
            # Parse distance matrix
            distances = [] # This part is complex, should be a helper
            for line in result.stdout.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 5:
                    seq1, seq2, distance = parts[0], parts[1], float(parts[2])
                    distances.append((seq1, seq2, distance))
            
            # Simple clustering based on distance threshold
            cluster_assignments = self._simple_clustering(distances, threshold)
            
            # Map sequence IDs to cluster IDs
            seq_id_map = {os.path.basename(seq_file): seq_id 
                         for seq_id, seq_file in enumerate(sequences.keys())}
            
            for seq_file, cluster_id in cluster_assignments.items():
                seq_id = seq_id_map.get(seq_file, seq_file)
                cluster_info[seq_id] = f'cluster_{cluster_id}'
            
            # Clean up temporary files
            for temp_file in [temp_fasta_path, sketch_path]:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
        
        except subprocess.CalledProcessError as e:
            self.logger.warning(f"Mash clustering failed: {e}")
            # Fall back to single cluster per sequence
            for i, seq_id in enumerate(sequences.keys()):
                cluster_info[seq_id] = f'cluster_{i+1}'
        
        return cluster_info
    
    def cluster_plasmids_skani(self, plasmid_file: str, output_dir: Path, sim_threshold: float = 99.95,
                              cov_threshold: float = 98.0, skip_cluster: bool = False) -> Tuple[List[SeqRecord], pd.DataFrame]:
        """
        Cluster plasmids using skani.
        
        Args:
            plasmid_file: Input plasmid FASTA file
            sim_threshold: ANI threshold for clustering
            cov_threshold: Alignment fraction threshold
            skip_cluster: Skip clustering (treat each as separate cluster)
            
        Returns:
            Tuple of (representative sequences, cluster data)
        """
        self.logger.info(f"Clustering plasmids with skani (threshold={sim_threshold}%)")
        from scipy.spatial.distance import squareform
        # Load sequences
        clus_seqs = list(SeqIO.parse(plasmid_file, 'fasta'))
        
        if skip_cluster:
            clus_data = pd.DataFrame({
                'Plasmid': [seq.id for seq in clus_seqs],
                'cluster': [str(i) for i in range(len(clus_seqs))],
                'len': [len(seq) for seq in clus_seqs],
                'pct_ident': [100.0] * len(clus_seqs),
                'Representative': [seq.id for seq in clus_seqs]
            })
            return clus_seqs, clus_data
            
        # Run skani
        skani_output = str(output_dir / 'skani_results.tsv')
        self.skani_runner.run_ska_dist(plasmid_file, skani_output)
        skani_df = pd.read_csv(skani_output, sep='\t')
        # Calculate combined distance and other metrics
        thres = 100.0 - sim_threshold
        skani_df['ref'] = skani_df['Ref_name'].map(lambda s: s.split()[0])
        skani_df['query'] = skani_df['Query_name'].map(lambda s: s.split()[0])
        skani_df['dANI'] = 100.0 - skani_df['ANI']
        skani_df['minAF'] = skani_df.apply(lambda r: min(r['Align_fraction_ref'], r['Align_fraction_query']), axis=1)
        skani_df['dAF'] = 100.0 - skani_df['minAF']
        af_mult = (100.0 - cov_threshold) / thres
        skani_df['combined_dist'] = (0.5 * (skani_df['dANI'])) + (0.5 * (skani_df['dAF']) / af_mult)
        skani_df.loc[skani_df['ref'] == skani_df['query'], 'combined_dist'] = 0
        
        # Create distance matrix
        all_ids = sorted(set(list(skani_df['ref']) + list(skani_df['query'])))
        skani_dm = skani_df[skani_df['minAF'] >= cov_threshold].pivot(
            index='query', columns='ref', values='combined_dist'
        ).fillna(100)
        
        skani_dm = skani_dm.reindex(index=all_ids, columns=all_ids)
        
        import itertools
        # Fill symmetric matrix
        for i1, i2 in itertools.combinations_with_replacement(all_ids, 2):
            if i1 == i2:
                skani_dm.at[i1, i2] = 0
            else:
                d1 = skani_dm.at[i1, i2] if i1 in skani_dm.index and i2 in skani_dm.columns else np.nan
                d2 = skani_dm.at[i2, i1] if i2 in skani_dm.index and i1 in skani_dm.columns else np.nan
                
                if not pd.isna(d1) and not pd.isna(d2):
                    max_dist = max(d1, d2)
                    skani_dm.at[i1, i2] = max_dist
                    skani_dm.at[i2, i1] = max_dist
                    
        skani_dm = skani_dm.fillna(100)
        
        # Perform clustering
        try:
            linkage_matrix = sch.linkage(squareform(skani_dm), method='average')
            clusters = sch.fcluster(linkage_matrix, t=thres, criterion='distance')
            
            clus = pd.DataFrame(clusters, index=skani_dm.index, columns=['cluster'])
            
            # Add sequence lengths
            for seq in clus_seqs:
                if seq.id in clus.index:
                    clus.at[seq.id, 'len'] = len(seq)
                    
            clus['Plasmid'] = clus.index
            
            # Select representatives (longest in each cluster)
            clus = clus.sort_values(['cluster', 'len', 'Plasmid'], ascending=[True, False, True])
            
            for _, df in clus.groupby('cluster'):
                rep = df.iloc[0]['Plasmid']
                clus.loc[df.index, 'Representative'] = rep
                
            # Get representative sequences
            rep_seq_ids = list(clus['Representative'].unique())
            rep_seqs = [s for s in clus_seqs if s.id in rep_seq_ids]
            
            self.logger.info(f"Clustered into {len(rep_seqs)} representative sequences")
            return rep_seqs, clus
            
        except Exception as e:
            self.logger.error(f"Clustering failed: {e}")
            return clus_seqs, pd.DataFrame()

    def _simple_clustering(self, distances: List[Tuple[str, str, float]],
                          threshold: float) -> Dict[str, int]:
        """Simple clustering based on distance threshold."""
        # Create adjacency list for sequences within threshold
        adjacency = {}
        all_seqs = set()
        
        # Ensure distances is a list of tuples/lists with at least 3 elements
        valid_distances = [d for d in distances if isinstance(d, (list, tuple)) and len(d) >= 3]

        for item in valid_distances:
            seq1, seq2, dist, *_ = item  # Unpack safely
            all_seqs.add(seq1)
            all_seqs.add(seq2)
            
            if dist <= (1.0 - threshold):  # Convert similarity to distance
                if seq1 not in adjacency:
                    adjacency[seq1] = []
                if seq2 not in adjacency:
                    adjacency[seq2] = []
                adjacency[seq1].append(seq2)
                adjacency[seq2].append(seq1)
        
        # Find connected components (clusters)
        visited = set()
        clusters = {}
        cluster_id = 1
        
        for seq in all_seqs:
            if seq not in visited:
                # BFS to find connected component
                cluster_members = []
                queue = [seq]
                
                while queue:
                    current = queue.pop(0)
                    if current not in visited:
                        visited.add(current)
                        cluster_members.append(current)
                        
                        # Add neighbors
                        if current in adjacency:
                            for neighbor in adjacency[current]:
                                if neighbor not in visited:
                                    queue.append(neighbor)
                
                # Assign cluster ID to all members
                for member in cluster_members:
                    clusters[member] = cluster_id
                cluster_id += 1
        
        return clusters
    
    def _combine_rep_types(self, rep_types: List[str]):
        rep_types_found = set()
        for types in rep_types:
            if type(types) != str:
                continue
            types = types.strip()
            clean_types = [t.strip().split('_')[0].split('(')[0].split('-')[0] for t in types.split(',')]
            rep_types_found.update(clean_types)
        if '-' in rep_types_found:
            rep_types_found.remove('-')
        self.logger.debug(f"Found rep types: {rep_types_found}")
        return ','.join(rep_types_found)

    def _merge_cluster_info(self, metadata_df: pd.DataFrame, 
                           cluster_info: pd.DataFrame,
                           mob_df: pd.DataFrame,
                           plasmidfinder_results: pd.DataFrame) -> pd.DataFrame:
        """Merge clustering information with metadata."""
        metadata_df = metadata_df.set_index('plasmid_id')

        if not cluster_info.empty:
            cluster_info = cluster_info.set_index('Plasmid')
            metadata_df['primary_cluster_id'] = cluster_info['cluster']
            metadata_df['representative'] = cluster_info['Representative']
            metadata_df = metadata_df.loc[cluster_info['Representative'].unique()]
            metadata_df['plasmid_id'] = metadata_df.index

        mob_df['plasmid_id'] = mob_df['sample_id'].map(lambda s: s.split()[0])
        mob_df = mob_df.set_index('plasmid_id')
        metadata_df = metadata_df.merge(mob_df, left_index=True, right_index=True, suffixes=['_x', ''], how='left')
        pf_summary = plasmidfinder_results.groupby('qseqid')['sseqid'].apply(lambda x: ','.join(x.unique())).reset_index()

        metadata_df['plasmidfinder_reps'] = metadata_df.index.map(lambda s: pf_summary.loc[s, 'sseqid'] if s in pf_summary.index else None)
        self.logger.debug("\n" + metadata_df.head().to_string())

        metadata_df['rep_types'] = metadata_df.apply(lambda r: self._combine_rep_types([r['rep_type(s)'], r['plasmidfinder_reps']]), axis=1)
        metadata_df = metadata_df.drop(columns=[c for c in metadata_df.columns if '_x' in c])
        metadata_df['sample_id'] = metadata_df.index

        return metadata_df
    
    def _run_mob_cluster(self, sequences, metadata_df, mob_typer_file, output_dir):
        """Run Mob Cluster on representative sequences."""
        tax_df = pd.DataFrame(index=metadata_df.index, columns=['id','organism'])
        tax_df['id'] = metadata_df.index
        tax_df['organism'] = 'Bacteria'
        tax_df.to_csv(str(Path(output_dir) / 'taxonomy.tsv'), sep='\t', index=False)
        self.logger.info("Running Mob Cluster...")
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta") as temp_fasta:
            for seq in sequences:
                seq.name = seq.id
                seq.description = seq.id
            SeqIO.write(sequences, temp_fasta, 'fasta')
            plasmids_file = temp_fasta.name
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".tsv") as temp_tsv:
            metadata_df.to_csv(temp_tsv, sep='\t', index=False, columns=MOB_TYPER_COLS)
            cluster_mob_typer = temp_tsv.name
        self.mob_runner.run_mob_cluster(plasmids_file, str(Path(output_dir) / 'mob_cluster'), cluster_mob_typer, str(Path(output_dir) / 'taxonomy.tsv'))
        mob_clusters = pd.read_csv(str(Path(output_dir) / 'mob_cluster' / 'clusters.txt' ), sep='\t')
        metadata_df.update(mob_clusters.set_index('sample_id'))
        metadata_df['cluster_id'] = metadata_df['secondary_cluster_id']
        return metadata_df

    def _write_database_files(self, sequences: Dict[str, SeqRecord], 
                             metadata_df: pd.DataFrame, output_dir: str) -> None:
        """Write database files to output directory."""
        os.makedirs(output_dir, exist_ok=True)
        
        # Write FASTA file
        fasta_path = os.path.join(output_dir, "plasmidDB.fasta")
        with open(fasta_path, 'w') as f:
            SeqIO.write(sequences, f, 'fasta')
        
        # Update sequence lengths in metadata
        for record in sequences:
            seq_id = record.id
            mask = metadata_df['plasmid_id'] == seq_id
            if mask.any():
                metadata_df.loc[mask, 'length'] = len(record.seq)
        
        # Write metadata file
        metadata_path = os.path.join(output_dir, "plasmidDB_info.csv")
        metadata_df.to_csv(metadata_path)
        
        # Write database statistics
        stats_path = os.path.join(output_dir, "database_stats.txt")
        with open(stats_path, 'w') as f:
            f.write(f"Accio Plasmid Database Statistics\n")
            f.write(f"=================================\n\n")
            f.write(f"Total sequences: {len(sequences)}\n")
            f.write(f"Total clusters: {len(set(metadata_df['cluster_id']))}\n")
            f.write(f"Average length: {metadata_df['length'].mean():.1f} bp\n")
            f.write(f"Length range: {metadata_df['length'].min()} - {metadata_df['length'].max()} bp\n")
            
            # Cluster size distribution
            cluster_sizes = metadata_df['cluster_id'].value_counts()
            f.write(f"\nCluster size distribution:\n")
            f.write(f"  1 sequence: {sum(cluster_sizes == 1)} clusters\n")
            f.write(f"  2-5 sequences: {sum((cluster_sizes >= 2) & (cluster_sizes <= 5))} clusters\n")
            f.write(f"  6-10 sequences: {sum((cluster_sizes >= 6) & (cluster_sizes <= 10))} clusters\n")
            f.write(f"  >10 sequences: {sum(cluster_sizes > 10)} clusters\n")
    
    def _create_plasmid_communities(self, sequences: Dict[str, SeqRecord], 
                                    metadata_df: pd.DataFrame, output_dir: str) -> pd.DataFrame:
        pling_runner = PlingRunner()
        pling_dir = str(Path(output_dir) / 'pling')
        fasta_dir = Path(output_dir) / 'fastas'
        fasta_dir.mkdir(exist_ok=True)
        for seq in sequences:
            SeqIO.write(seq, str(fasta_dir / f'{seq.id}.fasta'), 'fasta')
        fasta_path_file = Path(output_dir) / 'pling_fasta_paths.txt'
        with open(fasta_path_file, 'w') as f:
            f.write('\n'.join([str(fasta_dir / f'{seq.id}.fasta') for seq in sequences]))
        pling_runner.run_pling(str(fasta_path_file), pling_dir, pling_args='--containment_distance 0.03')
        new_dcj_types_dir = pling_runner.edit_dcj(metadata_df)
        
        # Assuming edit_dcj returns the path to the new typing results directory
        typing_file = Path(new_dcj_types_dir) / 'objects' / 'typing.tsv'
        if typing_file.exists():
            pling_types = pd.read_csv(typing_file, sep='\t', index_col=0, header=0, names=['pling_type'])
            hub_plasmids = pd.read_csv(str(Path(pling_dir) / 'dcj_thresh_4_graph' /'objects'/ 'hub_plasmids.csv'), index_col=0)
            communities = pd.read_csv(str(Path(pling_dir) / 'containment' /'containment_communities' /'objects'/ 'communities.tsv'), sep='\t', index_col=0, header=0, names=['community'])

            metadata_df = metadata_df.join(pling_types, how='left')
            metadata_df = metadata_df.join(communities, how='left')
            for i, h in enumerate(hub_plasmids.index):
                metadata_df.loc[h, 'pling_type'] = communities.loc[h, 'community'] + "_hub{i}"

        return metadata_df

    def _build_blast_database(self, output_dir: str) -> None:
        """Build BLAST database."""
        fasta_path = os.path.join(output_dir, "plasmidDB.fasta")
        try:
            self.blast_runner.run_makeblastdb(
                input_fasta=fasta_path,
                output_file=os.path.join(output_dir, 'plasmidDB'),
                title='Accio Plasmid Database'
            )
            self.logger.info("BLAST database created successfully")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to create BLAST database: {e}")
            raise
    
    def _build_mash_database(self, output_dir: str) -> None:
        """Build Mash database."""
        fasta_path = os.path.join(output_dir, "plasmidDB.fasta")
        
        self.mash_runner.run_sketch(fasta_path, fasta_path)
        self.logger.info('Mash database created successfully')
        
    
    def _get_database_files(self, output_dir: str) -> List[str]:
        """Get list of created database files."""
        base_files = [
            "plasmidDB.fasta",
            "plasmidDB_info.csv",
            "plasmidDB.fasta.msh",
            "database_stats.txt"
        ]
        
        # BLAST database files
        blast_extensions = ['.nhr', '.nin', '.nsq']
        for ext in blast_extensions:
            base_files.append(f"plasmidDB{ext}")
        
        # Filter to existing files
        existing_files = []
        for filename in base_files:
            filepath = os.path.join(output_dir, filename)
            if os.path.exists(filepath):
                existing_files.append(filename)
        
        return existing_files