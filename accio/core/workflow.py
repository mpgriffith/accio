"""Main workflow orchestrator for plasmid analysis."""

import logging
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Any
from Bio import SeqIO

from ..config import AnalysisConfig, MAX_PLASMID_LENGTH, FilePatterns
from ..wrappers.external_tools import (
    BLASTRunner, MashRunner, NucmerRunner, PLASMeRunner,
    ReadMappingRunner, MOBRunner, ExternalToolError
)
from .contig import Contig
from .plasmid import Plasmid
from .assignment import assign_plasmids
from ..utils.filters import filter_overlaps


class PlasmidAnalysisWorkflow:
    """Main workflow orchestrator for plasmid analysis."""
    
    def __init__(self, config: AnalysisConfig):
        """
        Initialize the workflow.
        
        Args:
            config: Analysis configuration
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Initialize tool runners
        self.blast_runner = BLASTRunner(config)
        self.mash_runner = MashRunner(config)
        self.nucmer_runner = NucmerRunner(config)
        self.plasme_runner = PLASMeRunner(config)
        self.mapping_runner = ReadMappingRunner(config)
        self.mob_runner = MOBRunner(config)
        
    def run(self, fasta_file: str, reads: Optional[List[str]], 
            db_dir: str, output_dir: str, 
            sample_name: Optional[str] = None,
            use_minimap: bool = False,
            keep_intermediate: bool = False,
            output_format: Optional[str] = None) -> Dict[str, Any]:
        """
        Run the complete plasmid analysis workflow.
        
        Args:
            fasta_file: Path to assembly FASTA file
            reads: List of read files (optional)
            db_dir: Database directory path
            output_dir: Output directory path
            use_minimap: Whether to use minimap2 for read mapping
            keep_intermediate: Whether to keep intermediate files
            
        Returns:
            Dictionary containing analysis results
        """
        self.logger.info("Starting plasmid analysis workflow")
        
        # Extract sample name
        sample_name = sample_name if sample_name else Path(fasta_file).stem
        
        # Load assembly and create contigs
        self.logger.info("Loading assembly contigs")
        contigs = self._load_contigs(fasta_file)
        
        # Load plasmid database information
        self.logger.info("Loading plasmid database")
        plasmids = self._load_plasmids(db_dir)
        
        # Run sequence analysis
        self.logger.info("Running sequence similarity analysis")
        blast_results = self._run_sequence_analysis(
            fasta_file, db_dir, output_dir, sample_name, keep_intermediate
        )
        if self.logger.isEnabledFor(logging.DEBUG):
            for name, df in blast_results.items():
                df.to_csv(Path(output_dir) / f'{sample_name}_{name}.blast.csv')
        # Run read-based analysis if reads provided
        coverage_data = None
        if reads:
            self.logger.info("Running read-based analysis")
            coverage_data = self._run_read_analysis(
                fasta_file, reads, db_dir, output_dir, sample_name, use_minimap
            )
        
        # Process contigs with analysis results
        self.logger.info("Processing contig analysis results")
        self._process_contigs(contigs, blast_results, coverage_data)
        
        # Process plasmids with analysis results
        self.logger.info("Processing plasmid analysis results")
        plasmid_matches = self._process_plasmids(contigs, plasmids, blast_results)
        
        # Assign contigs to plasmids
        self.logger.info("Assigning contigs to plasmids")
        assignment_results = self._assign_plasmids(contigs, plasmids)
        
        # Handle novel plasmid contigs
        self.logger.info("Processing novel plasmid contigs")
        novel_results = self._process_novel_plasmids(
            contigs, output_dir, sample_name
        )
        
        # Generate output files
        self.logger.info("Generating output files")
        output_files = self._generate_outputs(
            assignment_results, novel_results, contigs, 
            output_dir, sample_name
        )
        if not self.logger.isEnabledFor(logging.DEBUG):
            
            self._clean_up_files(output_dir, output_files.values())

        return {
            'contigs': contigs,
            'plasmids_assigned': assignment_results['plasmids_assigned'],
            'novel_plasmid_contigs': novel_results.get('contigs', []),
            'output_files': output_files,
            'summary': self._create_summary(assignment_results, novel_results)
        }
        
    def _load_contigs(self, fasta_file: str) -> Dict[str, Contig]:
        """Load contigs from FASTA file."""
        contigs = {}
        
        for seq in SeqIO.parse(fasta_file, 'fasta'):
            # Only consider contigs that could be plasmids
            if len(seq) < 1.5 * MAX_PLASMID_LENGTH:
                is_circular = 'circular=true' in seq.description.lower()
                contigs[seq.id] = Contig(seq.id, seq, length=len(seq), circular=is_circular)
                
        self.logger.info(f"Loaded {len(contigs)} potential plasmid contigs")
        return contigs
        
    def _load_plasmids(self, db_dir: str) -> Dict[str, Plasmid]:
        """Load plasmid database information."""
        info_file = Path(db_dir) / 'plasmidDB_info.csv'
        mash_db = Path(db_dir) / 'plasmidDB.fasta.msh'
        
        try:
            plasmid_info = pd.read_csv(info_file)
            plasmids = {}
            
            for idx, row in plasmid_info.iterrows():
                plasmids[row['sample_id']] = Plasmid(
                    id=row['sample_id'],
                    pling_type=row['pling_type'],
                    cluster_id=row['primary_cluster_id'],
                    length=row['size'],
                    plasmid_mash_db=mash_db,
                    other={'pling_type': row['pling_type']},
                    rep_type=row.get('rep_types', None)
                )
                
            self.logger.info(f"Loaded {len(plasmids)} reference plasmids")
            return plasmids
            
        except Exception as e:
            raise RuntimeError(f"Failed to load plasmid database: {e}")
            
    def _run_sequence_analysis(self, fasta_file: str, db_dir: str, 
                             output_dir: str, sample_name: str,
                             keep_intermediate: bool) -> Dict[str, pd.DataFrame]:
        """Run sequence-based analysis tools."""
        results = {}
        self.intermediate_files = []
        
        output_dir_path = Path(output_dir)
        # Define database paths
        db_dir_path = Path(db_dir)
        plasmid_db = str(db_dir_path / 'plasmidDB')
        plasmidfinder_db = str(db_dir_path / 'PlasmidFinder')
        rep_db = str(db_dir_path / 'repetitive.dna.fas')
        mash_db = str(db_dir_path / 'plasmidDB.fasta.msh')
        
        try:
            # BLAST against plasmid database
            self.logger.info("Running BLAST against plasmid database")
            output_file = output_dir_path / FilePatterns.BLAST_OUTPUT.format(sample=sample_name)
            blast_out_path = str(output_file) if keep_intermediate else None
            results['blast'] = self.blast_runner.run_blastn(
                fasta_file, plasmid_db, blast_out_path, sample_name=sample_name
            )
            if blast_out_path: self.intermediate_files.append(blast_out_path)
            
            # BLAST against PlasmidFinder
            self.logger.info("Running BLAST against PlasmidFinder database")
            output_file = output_dir_path / FilePatterns.PLASMIDFINDER_OUTPUT.format(sample=sample_name)
            pf_out_path = str(output_file) if keep_intermediate else None
            plasmidfinder_results = self.blast_runner.run_blastn(
                fasta_file, plasmidfinder_db, pf_out_path, sample_name=sample_name
            )
            if pf_out_path: self.intermediate_files.append(pf_out_path)
            
            # Process PlasmidFinder results
            if not plasmidfinder_results.empty:
                plasmidfinder_results['sim_score'] = (
                    plasmidfinder_results['pident'] * plasmidfinder_results['length'] / plasmidfinder_results['slen']
                )
                plasmidfinder_results = filter_overlaps(
                    plasmidfinder_results, 
                    sort_cols=['sim_score', 'bitscore'], 
                    sort_ascending=[False, False]
                )
                results['plasmidfinder'] = plasmidfinder_results[plasmidfinder_results['sim_score'] > 80]
            else:
                results['plasmidfinder'] = pd.DataFrame()
                
            # BLAST against repetitive elements
            self.logger.info("Running BLAST against repetitive elements")
            output_file = output_dir_path / FilePatterns.REP_BLAST_OUTPUT.format(sample=sample_name)
            rep_out_path = str(output_file) if keep_intermediate else None
            results['repetitive'] = self.blast_runner.run_blastn(
                fasta_file, rep_db, rep_out_path, sample_name=sample_name
            )
            if rep_out_path: self.intermediate_files.append(rep_out_path)
            
            # Mash screening
            self.logger.info("Running Mash screening")
            output_file = output_dir_path / FilePatterns.MASH_OUTPUT.format(sample=sample_name)
            mash_out_path = str(output_file) if keep_intermediate else None
            results['mash'] = self.mash_runner.run_screen(
                mash_db, fasta_file, mash_out_path
            )
            if mash_out_path: self.intermediate_files.append(mash_out_path)
            
            # Nucmer alignment
            self.logger.info("Running Nucmer alignment")
            prefix = output_dir_path / FilePatterns.NUCMER_PREFIX.format(sample=sample_name)
            plasmid_fasta = db_dir_path / 'plasmidDB.fasta'
            coords_df, snps_df, nucmer_files = self.nucmer_runner.run_nucmer(str(plasmid_fasta), fasta_file, str(prefix))
            results['nucmer_coords'] = coords_df
            results['nucmer_snps'] = snps_df
            if keep_intermediate:
                self.intermediate_files.extend(nucmer_files.values())
            
            # PLASMe analysis
            self.logger.info("Running PLASMe analysis")
            try:
                plasme_results, plasme_files = self.plasme_runner.run_plasme(fasta_file, str(Path(output_dir) / 'plasme'))
                results['plasme'] = plasme_results
                if keep_intermediate:
                    self.intermediate_files.extend(plasme_files.values())
            except ExternalToolError as e:
                self.logger.warning(f"PLASMe analysis failed: {e}")
                results['plasme'] = pd.DataFrame()
                
        except ExternalToolError as e:
            self.logger.error(f"Sequence analysis failed: {e}")
            raise
            
        return results
        
    def _run_read_analysis(self, fasta_file: str, reads: List[str], 
                          db_dir: str, output_dir: str, sample_name: str,
                          use_minimap: bool) -> Dict[str, Any]:
        """Run read-based analysis."""
        results = {}
        output_dir_path = Path(output_dir)
        
        try:
            # Map reads to assembly for coverage calculation
            self.logger.info("Mapping reads to assembly")
            assembly_bam = self.mapping_runner.map_reads(
                fasta_file, reads, 
                str(output_dir_path / f"{sample_name}_self"),
                use_minimap
            )
            
            # Calculate coverage statistics
            cov_file = output_dir_path / FilePatterns.CONTIG_COV_OUTPUT.format(sample=sample_name)
            coverage_df, cov_path = self.mapping_runner.calculate_coverage(assembly_bam, str(cov_file))
            
            # Calculate median coverage and copy numbers
            median_cov = self._get_median_coverage(coverage_df)
            coverage_df['copy_num'] = coverage_df['meandepth'] / median_cov
            
            results['coverage'] = coverage_df
            results['median_coverage'] = median_cov
            
            # Map reads to plasmid database
            self.logger.info("Mapping reads to plasmid database")
            plasmid_db_fasta = Path(db_dir) / 'plasmidDB.fasta'
            plasmid_bam = self.mapping_runner.map_reads(
                str(plasmid_db_fasta), reads,
                str(output_dir_path / f"{sample_name}_plasmidDB"),
                use_minimap
            )
            
            # Mash screening of reads
            # self.logger.info("Mash screening of reads")
            # mash_db = os.path.join(db_dir, 'plasmidDB.fasta.msh')
            
            # # Combine reads temporarily
            # combined_reads = os.path.join(output_dir, 'combined_reads.fastq')
            # with open(combined_reads, 'w') as outfile:
            #     for read_file in reads:
            #         with open(read_file, 'r') as infile:
            #             outfile.write(infile.read())
                        
            # read_mash_results = self.mash_runner.run_screen(mash_db, combined_reads)
            # os.remove(combined_reads)
            
            # results['read_mash'] = read_mash_results
            results['plasmid_bam'] = plasmid_bam
            
        except ExternalToolError as e:
            self.logger.error(f"Read analysis failed: {e}")
            raise
            
        return results
        
    def _get_median_coverage(self, coverage_df: pd.DataFrame, 
                           top_n: int = 10) -> float:
        """Calculate median coverage from top contigs."""
        if coverage_df.empty:
            return 30.0  # Default fallback
            
        # Use top N longest contigs for median calculation
        top_contigs = coverage_df.nlargest(top_n, 'length')
        
        if top_contigs.empty:
            return 30.0
            
        # Calculate length-weighted median
        total_length = top_contigs['length'].sum()
        cumulative_length = 0
        
        for _, row in top_contigs.iterrows():
            cumulative_length += row['length']
            if cumulative_length > total_length / 2:
                return row['meandepth']
                
        return top_contigs['meandepth'].median()
        
    def _process_contigs(self, contigs: Dict[str, Contig], 
                        blast_results: Dict[str, pd.DataFrame],
                        coverage_data: Optional[Dict[str, Any]]) -> None:
        """Process contigs with analysis results."""
        
        for contig_id, contig in contigs.items():
            # Add BLAST data
            contig_blast = blast_results['blast'][blast_results['blast']['qseqid'] == contig_id]
            contig.add_blast_data(contig_blast, self.config)
            
            # Add PLASMe data
            if not blast_results['plasme'].empty:
                plasme_data = blast_results['plasme'][blast_results['plasme']['contig'] == contig_id]
                contig.add_plasme_data(plasme_data)
                
            # Add replicon type data
            if not blast_results['plasmidfinder'].empty:
                rep_data = blast_results['plasmidfinder'][blast_results['plasmidfinder']['qseqid'] == contig_id]
                contig.add_rep_type(rep_data)
                
            # Add repetitive element data
            if not blast_results['repetitive'].empty:
                rep_elem_data = blast_results['repetitive'][blast_results['repetitive']['qseqid'] == contig_id]
                contig.add_repetitive_data(rep_elem_data, self.config)
                
            # Add Nucmer data
            coords_data = blast_results['nucmer_coords'][blast_results['nucmer_coords']['qseqid'] == contig_id]
            snps_data = blast_results['nucmer_snps'][blast_results['nucmer_snps']['query'] == contig_id]
            contig.add_nucmer_data(coords_data, snps_data)
            
            # Add coverage data if available
            if coverage_data:
                cov_df = coverage_data['coverage']
                contig_cov = cov_df[cov_df['contig'] == contig_id]
                if not contig_cov.empty:
                    contig.add_copy_num(contig_cov.iloc[0]['copy_num'], self.config)
                    
    def _process_plasmids(self, contigs: Dict[str, Contig], 
                         plasmids: Dict[str, Plasmid],
                         blast_results: Dict[str, pd.DataFrame]) -> List[str]:
        """Process plasmids and filter by replicon types."""
        
        # Get replicon types found in sample
        plasmidfinder_results = blast_results['plasmidfinder']
        if not plasmidfinder_results.empty:
            found_rep_types = [
                p.split('_')[0].split('(')[0].split('-')[0] 
                for p in plasmidfinder_results['sseqid']
            ]
        else:
            found_rep_types = []
            
        # Load plasmid database info for filtering
        plasmid_matches = set()
        filtered_plasmids = {}
        
        for pid, plasmid in plasmids.items():
            # Check if plasmid's replicon types are present in sample
            plasmid_rep_types = [rep_type.split('_')[0].split('(')[0].split('-')[0] for rep_type in plasmid.rep_type.split(',')] if type(plasmid.rep_type) == str else []

            if plasmid_rep_types:
                has_matching_rep = any(
                    rep_type.split('_')[0].split('(')[0] in found_rep_types
                    for rep_type in plasmid_rep_types
                )
            else:
                has_matching_rep = True  # Allow plasmids without specific rep types
                
            if has_matching_rep:
                filtered_plasmids[pid] = plasmid
                
                # Check if any contigs match this plasmid
                blast_data = blast_results['blast']
                if not blast_data.empty:
                    plasmid_hits = blast_data[blast_data['sseqid'] == pid]
                    if not plasmid_hits.empty:
                        plasmid_matches.add(pid)
                        
        # Update plasmids dict
        plasmids.clear()
        plasmids.update(filtered_plasmids)
        
        self.logger.info(f"Filtered to {len(filtered_plasmids)} plasmids with matching replicon types")
        self.logger.info(f"Found {len(plasmid_matches)} plasmids with BLAST matches")
        
        return list(plasmid_matches)
        
    def _assign_plasmids(self, contigs: Dict[str, Contig], 
                        plasmids: Dict[str, Plasmid]) -> Dict[str, Any]:
        """Assign contigs to plasmids using the assignment algorithm."""
        
        # Separate circular contigs for priority processing
        circular_contigs = {cid: c for cid, c in contigs.items() if c.circular}
        linear_contigs = {cid: c for cid, c in contigs.items() if not c.circular}
        
        all_assigned_plasmids = []
        all_assigned_contigs = {}
        
        # Process circular contigs first (higher confidence)
        if circular_contigs:
            self.logger.info(f"Processing {len(circular_contigs)} circular contigs")
            
            for contig_id, contig in circular_contigs.items():
                contig.contig_type = 'circular'
                
                # Find matching plasmids for this contig
                contig_plasmid_matches = {}
                
                for plasmid_id, plasmid in plasmids.items():
                    if not contig.blast_data.empty and plasmid_id in contig.blast_data['sseqid'].values:
                        # Create a copy of plasmid for this analysis
                        test_plasmid = Plasmid(
                            plasmid.id, plasmid.cluster_id, plasmid.pling_type, plasmid.length,
                            plasmid.plasmid_db_mash, plasmid.other, plasmid.rep_type
                        )
                        test_plasmid.add_contigs(contig)
                        test_plasmid.calc_scores(config=self.config)
                        
                        contig_plasmid_matches[plasmid_id] = test_plasmid
                        
                # Assign best matching plasmid
                if contig_plasmid_matches:
                    self.logger.debug(f'Circular contig {contig_id} matching to plasmids: {contig_plasmid_matches}')
                    assigned_plasmids, assigned_contigs = assign_plasmids(
                        contig_plasmid_matches, {contig_id: contig}
                    )
                    
                    all_assigned_plasmids.extend(assigned_plasmids)
                    all_assigned_contigs.update(assigned_contigs)
                    
        # Process linear contigs
        if linear_contigs:
            self.logger.info(f"Processing {len(linear_contigs)} linear contigs")
            
            # Add contigs to matching plasmids
            for plasmid_id, plasmid in plasmids.items():
                if plasmid_id in [p.id for p in all_assigned_plasmids]:
                    continue  # Skip already assigned plasmids
                    
                matching_contigs = []
                for contig_id, contig in linear_contigs.items():
                    if (contig_id not in all_assigned_contigs and
                        not contig.blast_data.empty and
                        plasmid_id in contig.blast_data['sseqid'].values):
                        matching_contigs.append(contig)
                        
                if matching_contigs:
                    plasmid.add_contigs(matching_contigs)
                    plasmid.calc_scores(config=self.config)
                    
            # Assign linear contigs to plasmids
            remaining_contigs = {
                cid: c for cid, c in linear_contigs.items() 
                if cid not in all_assigned_contigs
            }
            
            remaining_plasmids = {
                pid: p for pid, p in plasmids.items()
                if pid not in [pl.id for pl in all_assigned_plasmids]
            }
            
            assigned_plasmids, assigned_contigs = assign_plasmids(
                remaining_plasmids, remaining_contigs
            )
            
            all_assigned_plasmids.extend(assigned_plasmids)
            all_assigned_contigs.update(assigned_contigs)
            
        return {
            'plasmids_assigned': all_assigned_plasmids,
            'contigs_assigned': all_assigned_contigs
        }
        
    def _run_and_parse_mob_recon(self, novel_fasta: Path, output_dir: str):
        """
        Run MOB-recon on a FASTA file and parse the output reports.

        Args:
            novel_fasta: Path to the FASTA file with novel contigs.
            output_dir: The main output directory for the workflow.

        Returns:
            A tuple containing two DataFrames: (mob_typer_results, mob_contig_report).
            Returns empty DataFrames if the analysis fails or files are not found.
        """
        mob_output_dir = Path(output_dir) / 'mob_novel_plasmids'
        mob_results = pd.DataFrame()
        mob_contigs = pd.DataFrame()

        try:
            mob_files = self.mob_runner.run_mob_recon(str(novel_fasta), str(mob_output_dir))
            self.logger.info("MOB-recon analysis completed for novel plasmids.")

            # Define report file paths
            mob_report_file = Path(mob_files["mobtyper_results"])
            mob_contig_file = Path(mob_files["contig_report"])

            # Safely parse results
            if mob_report_file.exists() and mob_report_file.stat().st_size > 0 and mob_contig_file.exists() and mob_contig_file.stat().st_size > 0:
                mob_results = pd.read_csv(mob_report_file, sep='\t', on_bad_lines='skip')
                mob_contigs = pd.read_csv(mob_contig_file, sep='\t')
            else:
                self.logger.warning("MOB-recon did not produce the expected report files.")

        except ExternalToolError as e:
            self.logger.warning(f"MOB-recon processing failed: {e}")
        except Exception as e:
            self.logger.error(f"An unexpected error occurred during MOB-recon parsing: {e}")

        return mob_results, mob_contigs

    def _process_novel_plasmids(self, contigs: Dict[str, Contig],
                               output_dir: str, sample_name: str) -> Dict[str, Any]:
        """Process contigs that appear to be novel plasmids."""
        novel_contigs = [c for c in contigs.values() if c.type == 'plasmid' and not c.assigned]

        if novel_contigs:
            self.logger.info(f"Found {len(novel_contigs)} novel plasmid contigs")
            novel_fasta = Path(output_dir) / FilePatterns.NOVEL_PLASMIDS_OUTPUT.format(sample=sample_name)

            try:
                novel_seqs = [c.seq for c in novel_contigs]
                SeqIO.write(novel_seqs, novel_fasta, 'fasta')
                mob_results, mob_contigs = self._run_and_parse_mob_recon(novel_fasta, output_dir)
            except (IOError, OSError) as e:
                self.logger.error(f"Failed to write novel plasmids FASTA file: {e}")
                return {'contigs': novel_contigs, 'mob_results': pd.DataFrame(), 'contig_results': pd.DataFrame()}
        else:
            mob_results, mob_contigs = pd.DataFrame(), pd.DataFrame()

        return {'contigs': novel_contigs, 'mob_results': mob_results, 'contig_results': mob_contigs}
        
    def _generate_outputs(self, assignment_results: Dict[str, Any],
                         novel_results: Dict[str, Any],
                         contigs: Dict[str, Contig],
                         output_dir: str, sample_name: str) -> Dict[str, str]:
        """Generate final output files."""
        
        output_files = {}
        output_dir_path = Path(output_dir)
        
        # Generate plasmid summary
        plasmid_data = []
        
        # Add assigned plasmids
        for plasmid in assignment_results['plasmids_assigned']:
            summary = plasmid.get_summary()
            plasmid_data.append(summary)
            

        # Add novel plasmid entries from MOB-recon results
        novel_contigs = novel_results.get('contigs', [])
        mob_results = novel_results.get('mob_results', pd.DataFrame())
        mob_contigs = novel_results.get('contig_results', pd.DataFrame())

        if not mob_results.empty:
            # Group by plasmid ID from MOB-recon
            for plasmid_id, group in mob_contigs.groupby('sample_id'):
                contig_ids = group['contig_id'].tolist()

                # Find corresponding contig objects
                plasmid_contigs = [c for c in novel_contigs if c.id in contig_ids]

                if not plasmid_contigs:
                    continue

                rep_types = set()
                for contig in plasmid_contigs:
                    if contig.rep_types:
                        rep_types.update(contig.rep_types)

                # Get data from the first row of the group, assuming it's representative
                rep_row = group.iloc[0]

                novel_entry = {
                    'plasmid_id': plasmid_id,
                    'cluster_id': 'novel',
                    'pling_type': 'novel',
                    'length': sum(c.length for c in plasmid_contigs),
                    'num_contigs': len(plasmid_contigs),
                    'contig_ids': [c.id for c in plasmid_contigs],
                    'rep_types': ','.join(sorted(list(rep_types))),
                    'mob_plasmid_id': rep_row.get('plasmid_id'),
                    'mob_cluster_id': rep_row.get('primary_cluster_id'),
                    'mob_rep_types': rep_row.get('rep_type(s)'),
                    'mob_mash_nearest_neighbor': rep_row.get('mash_nearest_neighbor'),
                    'mob_mash_dist': rep_row.get('mash_dist'),
                }
                plasmid_data.append(novel_entry)
        novel_unclass_contigs = list(set(novel_contigs) - set(mob_contigs.contig_id.tolist())) if not mob_contigs.empty else set(novel_contigs)
        if novel_unclass_contigs:
            # Fallback if MOB-recon fails but novel contigs are found
            novel_entry = {
                'plasmid_id': 'novel_plasmid_1',
                'cluster_id': 'novel',
                'pling_type': 'novel',
                'length': sum(c.length for c in novel_unclass_contigs),
                'num_contigs': len(novel_unclass_contigs),
                'contig_ids': [c.id for c in novel_unclass_contigs],
                'rep_types': ','.join(set(
                    rep for c in novel_unclass_contigs if c.rep_types
                    for rep in c.rep_types
                ))
            }
            plasmid_data.append(novel_entry)
            

        # Write plasmid summary

        if True: #plasmid_data:
            plasmid_file = output_dir_path / FilePatterns.PLASMID_DATA_OUTPUT.format(sample=sample_name)

            pd.DataFrame(plasmid_data).to_csv(str(plasmid_file), index=False)
            output_files['plasmids'] = str(plasmid_file)
        
            

        # Generate contig summary
        contig_data = []
        for contig in contigs.values():
            summary = contig.get_summary()
            summary['novel_plasmid_contig'] = (
                contig.type == 'plasmid' and not contig.assigned
            )
            contig_data.append(summary)
            
        contig_file = output_dir_path / FilePatterns.CONTIG_DATA_OUTPUT.format(sample=sample_name)
        pd.DataFrame(contig_data).to_csv(str(contig_file), index=False)
        output_files['contigs'] = str(contig_file)
        
        if plasmid_data:
            for row in plasmid_data:
                seqs = []
                for contig_id in row['contig_ids']:
                    seqs.append(contigs[contig_id].seq)
                fasta_output = output_dir_path / f"{sample_name}_plasmid_{row['pling_type']}_contigs.fasta"
                SeqIO.write(seqs, str(fasta_output), 'fasta')
                output_files[f"plasmid_{row['pling_type']}_seq"] = str(fasta_output)

        return output_files

    def _clean_up_files(self, output_dir: str, output_files_to_keep: List[str]):
        """Clean up temporary files."""
        all_output_files = Path(output_dir).glob('*')
        for f in all_output_files:
            if str(f) not in output_files_to_keep:
                if f.is_dir():
                    for fd in f.glob('*'):
                        fd.unlink(missing_ok=True)
                    f.rmdir()
                else:
                    f.unlink(missing_ok=True)

    def _create_summary(self, assignment_results: Dict[str, Any],
                       novel_results: Dict[str, Any]) -> Dict[str, Any]:
        """Create analysis summary."""
        
        assigned_plasmids = assignment_results['plasmids_assigned']
        novel_contigs = novel_results.get('contigs', [])
        
        return {
            'num_plasmids_assigned': len(assigned_plasmids),
            'num_novel_plasmid_contigs': len(novel_contigs),
            'assigned_plasmid_ids': [p.id for p in assigned_plasmids],
            'novel_contig_ids': [c.id for c in novel_contigs],
            'total_assigned_length': sum(p.length for p in assigned_plasmids),
            'total_novel_length': sum(c.length for c in novel_contigs)
        }