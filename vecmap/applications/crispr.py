"""
CRISPR Guide Detection for Single-Cell Screens
==============================================

Optimized exact matching for Perturb-seq, CROP-seq, and similar assays.

Key advantages over traditional aligners:
1. No tolerance for mismatches (critical for guide assignment)
2. Handles millions of reads in seconds
3. Direct integration with single-cell Python ecosystem
4. Minimal memory footprint
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Set
from collections import defaultdict
from ..core.mapper import vecmap

class CRISPRGuideDetector:
    """
    High-performance CRISPR guide detection for single-cell screens.
    
    Designed for:
    - Perturb-seq / CROP-seq experiments
    - Guide RNA detection in pooled screens
    - Barcode-guide association in single cells
    """
    
    def __init__(self, guide_library: Dict[str, str], guide_length: int = 20):
        """
        Initialize detector with guide library.
        
        Args:
            guide_library: Dict mapping guide names to sequences
            guide_length: Expected length of guide sequences (default: 20)
        """
        self.guide_library = guide_library
        self.guide_length = guide_length
        
        # Build reference from all guides
        self.guide_positions = {}
        self.reference = ""
        
        position = 0
        for guide_name, guide_seq in guide_library.items():
            if len(guide_seq) != guide_length:
                raise ValueError(f"Guide {guide_name} has length {len(guide_seq)}, expected {guide_length}")
            
            self.guide_positions[position] = guide_name
            self.reference += guide_seq + "N" * 10  # Add spacer
            position += guide_length + 10
    
    def detect_guides(self, reads: List[Tuple[str, str]], 
                     allow_reverse_complement: bool = True) -> Dict[str, List[str]]:
        """
        Detect guide RNAs in reads using exact matching.
        
        Args:
            reads: List of (read_sequence, read_id) tuples
            allow_reverse_complement: Also search reverse complement
            
        Returns:
            Dict mapping read_id to list of detected guide names
        """
        results = defaultdict(list)
        
        # Forward strand detection
        alignments = vecmap(self.reference, reads, self.guide_length)
        
        for (pos, mismatch_count, read_id) in alignments:
            if mismatch_count == 0:  # Exact match only
                # Find which guide this position corresponds to
                guide_position = (pos // (self.guide_length + 10)) * (self.guide_length + 10)
                if guide_position in self.guide_positions:
                    results[read_id].append(self.guide_positions[guide_position])
        
        # Reverse complement detection
        if allow_reverse_complement:
            rc_reads = [(self._reverse_complement(seq), read_id) 
                       for seq, read_id in reads]
            
            rc_alignments = vecmap(self.reference, rc_reads, self.guide_length)
            
            for (pos, mismatch_count, read_id) in rc_alignments:
                if mismatch_count == 0:
                    guide_position = (pos // (self.guide_length + 10)) * (self.guide_length + 10)
                    if guide_position in self.guide_positions:
                        guide_name = self.guide_positions[guide_position]
                        if guide_name not in results[read_id]:
                            results[read_id].append(guide_name + "_RC")
        
        return dict(results)
    
    def detect_guides_with_context(self, reads: List[Tuple[str, str]], 
                                 upstream_context: str = "", 
                                 downstream_context: str = "") -> Dict[str, List[str]]:
        """
        Detect guides with sequence context (e.g., constant regions).
        
        Useful for:
        - Validating guide context in expression vectors
        - Filtering false positives
        - Detecting truncated guides
        """
        context_reference = ""
        context_positions = {}
        position = 0
        
        for guide_name, guide_seq in self.guide_library.items():
            full_seq = upstream_context + guide_seq + downstream_context
            context_positions[position + len(upstream_context)] = guide_name
            context_reference += full_seq + "N" * 10
            position += len(full_seq) + 10
        
        # Search for full context
        results = defaultdict(list)
        
        search_len = len(upstream_context) + self.guide_length + len(downstream_context)
        alignments = vecmap(context_reference, reads, search_len)
        
        for (pos, mismatch_count, read_id) in alignments:
            if mismatch_count == 0:
                # Check if this is a valid guide position
                for guide_pos, guide_name in context_positions.items():
                    if abs(pos - guide_pos) < 5:  # Allow small position variation
                        results[read_id].append(guide_name)
                        break
        
        return dict(results)
    
    def _reverse_complement(self, seq: str) -> str:
        """Compute reverse complement of DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement.get(base, 'N') for base in seq[::-1])
    
    def summarize_detection(self, detection_results: Dict[str, List[str]]) -> Dict[str, int]:
        """
        Summarize guide detection results.
        
        Returns:
            Dict mapping guide names to read counts
        """
        guide_counts = defaultdict(int)
        
        for read_guides in detection_results.values():
            for guide in read_guides:
                guide_counts[guide] += 1
        
        return dict(guide_counts)


class BarcodeGuideMatcher:
    """
    Match cell barcodes to guide RNAs for single-cell CRISPR screens.
    
    Handles the common workflow:
    1. Extract cell barcode from Read 1
    2. Extract guide sequence from Read 2
    3. Build barcode->guide mapping
    """
    
    def __init__(self, barcode_length: int = 16, umi_length: int = 10,
                 guide_detector: Optional[CRISPRGuideDetector] = None):
        """
        Initialize barcode-guide matcher.
        
        Args:
            barcode_length: Length of cell barcode (default: 16 for 10x)
            umi_length: Length of UMI (default: 10 for 10x)
            guide_detector: Pre-configured guide detector
        """
        self.barcode_length = barcode_length
        self.umi_length = umi_length
        self.guide_detector = guide_detector
    
    def process_read_pairs(self, 
                          read1_data: List[Tuple[str, str]], 
                          read2_data: List[Tuple[str, str]]) -> Dict[str, Dict[str, int]]:
        """
        Process paired-end reads to extract barcode-guide associations.
        
        Args:
            read1_data: List of (sequence, read_id) from Read 1 (barcodes)
            read2_data: List of (sequence, read_id) from Read 2 (guides)
            
        Returns:
            Dict mapping cell barcodes to guide counts
        """
        # Extract barcodes from Read 1
        barcode_map = {}
        for seq, read_id in read1_data:
            if len(seq) >= self.barcode_length + self.umi_length:
                barcode = seq[:self.barcode_length]
                umi = seq[self.barcode_length:self.barcode_length + self.umi_length]
                barcode_map[read_id] = (barcode, umi)
        
        # Detect guides in Read 2
        if self.guide_detector:
            guide_results = self.guide_detector.detect_guides(read2_data)
        else:
            # If no detector provided, extract fixed position
            guide_results = {}
            for seq, read_id in read2_data:
                if len(seq) >= 20:
                    guide_results[read_id] = [seq[:20]]  # First 20bp as guide
        
        # Build barcode->guide mapping
        barcode_guide_counts = defaultdict(lambda: defaultdict(int))
        
        for read_id, guides in guide_results.items():
            if read_id in barcode_map:
                barcode, umi = barcode_map[read_id]
                for guide in guides:
                    # Count unique UMIs per guide per barcode
                    barcode_guide_counts[barcode][guide] += 1
        
        return {bc: dict(guides) for bc, guides in barcode_guide_counts.items()}
    
    def filter_barcodes(self, 
                       barcode_guide_counts: Dict[str, Dict[str, int]], 
                       whitelist: Optional[Set[str]] = None,
                       min_umi_count: int = 3) -> Dict[str, Dict[str, int]]:
        """
        Filter results to valid cell barcodes.
        
        Args:
            barcode_guide_counts: Raw barcode->guide mapping
            whitelist: Set of valid cell barcodes (e.g., from 10x)
            min_umi_count: Minimum UMIs to call a guide
            
        Returns:
            Filtered barcode->guide mapping
        """
        filtered = {}
        
        for barcode, guide_counts in barcode_guide_counts.items():
            # Check whitelist
            if whitelist and barcode not in whitelist:
                continue
            
            # Filter guides by UMI count
            filtered_guides = {
                guide: count for guide, count in guide_counts.items()
                if count >= min_umi_count
            }
            
            if filtered_guides:
                filtered[barcode] = filtered_guides
        
        return filtered


def detect_crispr_guides(fastq_file: str, guide_library: Dict[str, str], 
                        max_reads: Optional[int] = None) -> Dict[str, int]:
    """
    Convenience function for guide detection from FASTQ.
    
    Args:
        fastq_file: Path to FASTQ file
        guide_library: Dict mapping guide names to sequences
        max_reads: Maximum reads to process (None for all)
        
    Returns:
        Dict mapping guide names to read counts
    """
    detector = CRISPRGuideDetector(guide_library)
    
    # Parse FASTQ and collect reads
    reads = []
    with open(fastq_file, 'r') as f:
        count = 0
        while True:
            header = f.readline().strip()
            if not header:
                break
            
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            
            read_id = header.split()[0][1:]  # Remove @
            reads.append((seq, read_id))
            
            count += 1
            if max_reads and count >= max_reads:
                break
    
    # Detect guides
    results = detector.detect_guides(reads)
    
    # Summarize
    return detector.summarize_detection(results) 