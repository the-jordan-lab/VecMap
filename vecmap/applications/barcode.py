"""
Barcode and UMI Processing for Single-Cell Genomics
===================================================

Ultra-fast exact matching for:
- Cell barcode identification
- UMI deduplication  
- Sample demultiplexing (hashtags, CMOs)
- Feature barcoding (antibodies, CRISPR guides)
"""

import numpy as np
from typing import List, Dict, Tuple, Set, Optional
from collections import defaultdict, Counter
from ..core.mapper import vecmap


class BarcodeProcessor:
    """
    High-performance barcode processing for single-cell assays.
    
    Optimized for 10x Genomics, Parse Biosciences, and similar platforms.
    """
    
    def __init__(self, 
                 barcode_whitelist: Optional[Set[str]] = None,
                 barcode_length: int = 16,
                 umi_length: int = 10,
                 max_hamming_dist: int = 1):
        """
        Initialize barcode processor.
        
        Args:
            barcode_whitelist: Set of valid barcodes (e.g., 10x whitelist)
            barcode_length: Expected barcode length
            umi_length: Expected UMI length
            max_hamming_dist: Maximum Hamming distance for barcode correction
        """
        self.barcode_whitelist = barcode_whitelist
        self.barcode_length = barcode_length
        self.umi_length = umi_length
        self.max_hamming_dist = max_hamming_dist
        
        # Build reference for exact matching if whitelist provided
        if barcode_whitelist:
            self._build_barcode_reference()
    
    def _build_barcode_reference(self):
        """Build concatenated reference for VecMap matching."""
        self.barcode_to_position = {}
        self.position_to_barcode = {}
        
        # Sort for consistent ordering
        if not self.barcode_whitelist:
            return
        sorted_barcodes = sorted(self.barcode_whitelist)
        
        # Build reference with spacers
        self.barcode_reference = ""
        position = 0
        
        for barcode in sorted_barcodes:
            self.barcode_to_position[barcode] = position
            self.position_to_barcode[position] = barcode
            self.barcode_reference += barcode + "N" * 10
            position += len(barcode) + 10
    
    def extract_barcodes(self, reads: List[Tuple[str, str]], 
                        barcode_start: int = 0) -> Dict[str, str]:
        """
        Extract barcodes from reads.
        
        Args:
            reads: List of (sequence, read_id) tuples
            barcode_start: Start position of barcode in read
            
        Returns:
            Dict mapping read_id to extracted barcode
        """
        barcode_map = {}
        
        for seq, read_id in reads:
            if len(seq) >= barcode_start + self.barcode_length:
                barcode = seq[barcode_start:barcode_start + self.barcode_length]
                barcode_map[read_id] = barcode
        
        return barcode_map
    
    def correct_barcodes(self, barcode_map: Dict[str, str]) -> Dict[str, str]:
        """
        Correct barcodes to whitelist using exact matching.
        
        This is much faster than traditional Hamming distance calculation
        for large whitelists.
        """
        if not self.barcode_whitelist:
            return barcode_map
        
        corrected = {}
        
        # Prepare reads for VecMap
        barcode_reads = [(barcode, read_id) 
                        for read_id, barcode in barcode_map.items()]
        
        # Use VecMap for exact matching
        alignments = vecmap(self.barcode_reference, barcode_reads, self.barcode_length)
        
        # Process results
        for pos, mismatches, read_id in alignments:
            if mismatches <= self.max_hamming_dist:
                # Find which barcode this position corresponds to
                barcode_pos = (pos // (self.barcode_length + 10)) * (self.barcode_length + 10)
                if barcode_pos in self.position_to_barcode:
                    corrected[read_id] = self.position_to_barcode[barcode_pos]
        
        return corrected
    
    def extract_umis(self, reads: List[Tuple[str, str]], 
                    umi_start: Optional[int] = None) -> Dict[str, str]:
        """
        Extract UMIs from reads.
        
        Args:
            reads: List of (sequence, read_id) tuples
            umi_start: Start position of UMI (defaults to after barcode)
            
        Returns:
            Dict mapping read_id to UMI
        """
        if umi_start is None:
            umi_start = self.barcode_length
        
        umi_map = {}
        
        for seq, read_id in reads:
            if len(seq) >= umi_start + self.umi_length:
                umi = seq[umi_start:umi_start + self.umi_length]
                umi_map[read_id] = umi
        
        return umi_map
    
    def deduplicate_umis(self, 
                        barcode_umi_gene_tuples: List[Tuple[str, str, str]]) -> Dict[str, Dict[str, int]]:
        """
        Deduplicate UMIs per barcode per gene.
        
        Args:
            barcode_umi_gene_tuples: List of (barcode, umi, gene) tuples
            
        Returns:
            Dict mapping barcode -> gene -> unique UMI count
        """
        # Group by barcode and gene
        barcode_gene_umis = defaultdict(lambda: defaultdict(set))
        
        for barcode, umi, gene in barcode_umi_gene_tuples:
            barcode_gene_umis[barcode][gene].add(umi)
        
        # Count unique UMIs
        counts = {}
        for barcode, gene_umis in barcode_gene_umis.items():
            counts[barcode] = {
                gene: len(umis) for gene, umis in gene_umis.items()
            }
        
        return counts


class HashtagDemultiplexer:
    """
    Demultiplex samples using hashtag oligonucleotides (HTOs) or CMOs.
    
    Common in:
    - Cell hashing (CITE-seq)
    - Multiplexed scRNA-seq
    - Sample pooling experiments
    """
    
    def __init__(self, hashtag_sequences: Dict[str, str]):
        """
        Initialize with hashtag sequences.
        
        Args:
            hashtag_sequences: Dict mapping sample names to hashtag sequences
        """
        self.hashtag_sequences = hashtag_sequences
        self.hashtag_length = len(next(iter(hashtag_sequences.values())))
        
        # Build reference for matching
        self._build_hashtag_reference()
    
    def _build_hashtag_reference(self):
        """Build reference for VecMap matching."""
        self.hashtag_positions = {}
        self.reference = ""
        
        position = 0
        for sample_name, hashtag_seq in self.hashtag_sequences.items():
            self.hashtag_positions[position] = sample_name
            self.reference += hashtag_seq + "N" * 20  # Spacer
            position += len(hashtag_seq) + 20
    
    def demultiplex_cells(self, 
                         hashtag_reads: List[Tuple[str, str]], 
                         cell_barcodes: Dict[str, str]) -> Dict[str, str]:
        """
        Assign cells to samples based on hashtag reads.
        
        Args:
            hashtag_reads: List of (sequence, read_id) from hashtag library
            cell_barcodes: Dict mapping read_id to cell barcode
            
        Returns:
            Dict mapping cell barcode to sample name
        """
        # Detect hashtags in reads
        alignments = vecmap(self.reference, hashtag_reads, self.hashtag_length)
        
        # Count hashtags per cell
        cell_hashtag_counts = defaultdict(lambda: defaultdict(int))
        
        for pos, mismatches, read_id in alignments:
            if mismatches == 0 and read_id in cell_barcodes:
                cell_barcode = cell_barcodes[read_id]
                
                # Find which hashtag this is
                hashtag_pos = (pos // (self.hashtag_length + 20)) * (self.hashtag_length + 20)
                if hashtag_pos in self.hashtag_positions:
                    sample_name = self.hashtag_positions[hashtag_pos]
                    cell_hashtag_counts[cell_barcode][sample_name] += 1
        
        # Assign cells to samples (simple maximum)
        cell_assignments = {}
        for cell_barcode, hashtag_counts in cell_hashtag_counts.items():
            if hashtag_counts:
                # Assign to sample with most reads
                assigned_sample = max(hashtag_counts.items(), key=lambda x: x[1])[0]
                cell_assignments[cell_barcode] = assigned_sample
        
        return cell_assignments
    
    def identify_doublets(self, 
                         cell_hashtag_counts: Dict[str, Dict[str, int]], 
                         min_ratio: float = 0.2) -> Set[str]:
        """
        Identify potential doublets based on multiple hashtags.
        
        Args:
            cell_hashtag_counts: Dict mapping cell -> sample -> count
            min_ratio: Minimum ratio of second hashtag to call doublet
            
        Returns:
            Set of cell barcodes identified as doublets
        """
        doublets = set()
        
        for cell_barcode, hashtag_counts in cell_hashtag_counts.items():
            if len(hashtag_counts) >= 2:
                sorted_counts = sorted(hashtag_counts.values(), reverse=True)
                if sorted_counts[1] / sorted_counts[0] >= min_ratio:
                    doublets.add(cell_barcode)
        
        return doublets


class FeatureBarcodeDetector:
    """
    Detect feature barcodes (antibodies, guide RNAs, etc.) in CITE-seq and similar assays.
    """
    
    def __init__(self, feature_barcodes: Dict[str, str]):
        """
        Initialize with feature barcode sequences.
        
        Args:
            feature_barcodes: Dict mapping feature names to barcode sequences
        """
        self.feature_barcodes = feature_barcodes
        self.barcode_length = len(next(iter(feature_barcodes.values())))
        
        # Build reference
        self._build_feature_reference()
    
    def _build_feature_reference(self):
        """Build reference for VecMap matching."""
        self.feature_positions = {}
        self.reference = ""
        
        position = 0
        for feature_name, barcode_seq in self.feature_barcodes.items():
            self.feature_positions[position] = feature_name
            self.reference += barcode_seq + "N" * 15
            position += len(barcode_seq) + 15
    
    def detect_features(self, 
                       feature_reads: List[Tuple[str, str]], 
                       allow_mismatches: int = 1) -> Dict[str, List[str]]:
        """
        Detect feature barcodes in reads.
        
        Args:
            feature_reads: List of (sequence, read_id) tuples
            allow_mismatches: Maximum mismatches allowed
            
        Returns:
            Dict mapping read_id to list of detected features
        """
        # Use VecMap for detection
        alignments = vecmap(self.reference, feature_reads, self.barcode_length)
        
        # Process results
        read_features = defaultdict(list)
        
        for pos, mismatches, read_id in alignments:
            if mismatches <= allow_mismatches:
                # Find which feature this is
                feature_pos = (pos // (self.barcode_length + 15)) * (self.barcode_length + 15)
                if feature_pos in self.feature_positions:
                    feature_name = self.feature_positions[feature_pos]
                    read_features[read_id].append(feature_name)
        
        return dict(read_features)


def process_10x_data(r1_fastq: str, r2_fastq: str,
                    barcode_whitelist: Set[str],
                    feature_reference: Optional[Dict[str, str]] = None,
                    max_reads: Optional[int] = None) -> Dict[str, Dict[str, int]]:
    """Simple 10x Genomics data processing pipeline.

    This utility parses paired FASTQ files, extracts barcodes and UMIs
    from ``r1_fastq`` and optionally detects feature barcodes in ``r2_fastq``.

    Args:
        r1_fastq: Path to Read 1 FASTQ containing cell barcodes and UMIs.
        r2_fastq: Path to Read 2 FASTQ containing cDNA or feature reads.
        barcode_whitelist: Set of valid cell barcodes.
        feature_reference: Optional feature barcode reference mapping feature
            names to sequences.
        max_reads: Optional limit on reads processed.

    Returns:
        Dict mapping cell barcodes to feature/gene counts.
    """
    processor = BarcodeProcessor(
        barcode_whitelist=barcode_whitelist,
        barcode_length=16,
        umi_length=10,
    )

    r1_reads: List[Tuple[str, str]] = []
    r2_reads: List[Tuple[str, str]] = []

    with open(r1_fastq, "r") as f1, open(r2_fastq, "r") as f2:
        read_count = 0
        while True:
            h1 = f1.readline().strip()
            if not h1:
                break
            if not h1.startswith('@'):
                raise ValueError(f"Invalid FASTQ header format in R1 at read {read_count + 1}")
            s1 = f1.readline().strip()
            if not s1:
                raise ValueError(f"Incomplete FASTQ record in R1 at read {read_count + 1}")
            f1.readline()
            f1.readline()
            h2 = f2.readline().strip()
            if not h2:
                raise ValueError(f"R2 file has fewer reads than R1 at read {read_count + 1}")
            if not h2.startswith('@'):
                raise ValueError(f"Invalid FASTQ header format in R2 at read {read_count + 1}")
            s2 = f2.readline().strip()
            if not s2:
                raise ValueError(f"Incomplete FASTQ record in R2 at read {read_count + 1}")
            f2.readline()
            f2.readline()

            read_id = h1[1:]
            r1_reads.append((s1, read_id))
            r2_reads.append((s2, read_id))

            read_count += 1
            if max_reads and read_count >= max_reads:
                break
    barcode_map = processor.extract_barcodes(r1_reads)
    corrected = processor.correct_barcodes(barcode_map)
    umi_map = processor.extract_umis(r1_reads)

    if feature_reference:
        detector = FeatureBarcodeDetector(feature_reference)
        feature_hits = detector.detect_features(r2_reads)
    else:
        feature_hits = {read_id: [] for _, read_id in r1_reads}

    tuples = []
    for _, read_id in r1_reads:
        barcode = corrected.get(read_id)
        umi = umi_map.get(read_id)
        if not (barcode and umi):
            continue
        feats = feature_hits.get(read_id) or ["gene"]
        for feat in feats:
            tuples.append((barcode, umi, feat))

    return processor.deduplicate_umis(tuples)