#!/usr/bin/env python3
"""
VecMap Production Version
========================
Enhanced with paired-end support, indel handling, and splice-aware alignment
"""

import numpy as np
from collections import defaultdict
import time
import random
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass

@dataclass
class AlignmentResult:
    """Store alignment results with full CIGAR support"""
    query_id: str
    ref_pos: int
    mapq: int
    cigar: str
    mismatches: int
    is_paired: bool = False
    mate_pos: Optional[int] = None
    insert_size: Optional[int] = None
    is_reverse: bool = False
    is_splice_junction: bool = False
    
class VecMapProduction:
    """Production-ready VecMap with full feature support"""
    
    def __init__(self, reference: str, 
                 seed_len: int = 20,
                 seed_step: int = 10,
                 max_mismatches: int = 5,
                 max_indel_size: int = 10,
                 splice_motifs: Optional[List[str]] = None):
        """
        Initialize VecMap with configurable parameters
        
        Args:
            reference: Reference sequence
            seed_len: Length of k-mer seeds
            seed_step: Step size for seed extraction
            max_mismatches: Maximum allowed mismatches
            max_indel_size: Maximum indel size to consider
            splice_motifs: Known splice site motifs (e.g., ['GT-AG', 'GC-AG'])
        """
        self.reference = reference
        self.ref_len = len(reference)
        self.seed_len = seed_len
        self.seed_step = seed_step
        self.max_mismatches = max_mismatches
        self.max_indel_size = max_indel_size
        self.splice_motifs = splice_motifs or ['GT-AG', 'GC-AG', 'AT-AC']
        
        # Build indices
        self.seed_index = self._build_seed_index()
        self.ref_array = np.array(list(reference))
        
        # Splice site index for RNA-seq
        self.splice_sites = self._find_splice_sites() if splice_motifs else {}
        
    def _build_seed_index(self) -> Dict[str, List[int]]:
        """Build k-mer index with configurable seed length"""
        index = defaultdict(list)
        for i in range(self.ref_len - self.seed_len + 1):
            seed = self.reference[i:i + self.seed_len]
            index[seed].append(i)
        return dict(index)
    
    def _find_splice_sites(self) -> Dict[int, List[Tuple[int, str]]]:
        """Find potential splice sites based on motifs"""
        splice_sites = defaultdict(list)
        
        for motif in self.splice_motifs:
            donor, acceptor = motif.split('-')
            
            # Find donor sites
            for i in range(self.ref_len - len(donor)):
                if self.reference[i:i+len(donor)] == donor:
                    # Look for nearby acceptor sites
                    for j in range(i + 50, min(i + 100000, self.ref_len - len(acceptor))):
                        if self.reference[j:j+len(acceptor)] == acceptor:
                            splice_sites[i].append((j, motif))
                            
        return dict(splice_sites)
    
    def _vectorized_scoring_with_indels(self, read_seq: str, candidates: List[int]) -> Tuple[np.ndarray, List[str]]:
        """
        Enhanced scoring that considers indels
        Returns scores and CIGAR strings
        """
        if not candidates:
            return np.array([]), []
            
        read_len = len(read_seq)
        read_arr = np.array(list(read_seq))
        scores = []
        cigars = []
        
        for start_pos in candidates:
            if start_pos < 0 or start_pos + read_len > self.ref_len:
                scores.append(float('inf'))
                cigars.append('')
                continue
                
            # Try exact match first
            ref_segment = self.ref_array[start_pos:start_pos + read_len]
            exact_mismatches = (ref_segment != read_arr).sum()
            
            if exact_mismatches <= self.max_mismatches:
                scores.append(exact_mismatches)
                cigars.append(f'{read_len}M')
                continue
            
            # Try with indels using banded alignment
            score, cigar = self._banded_alignment(read_seq, start_pos)
            scores.append(score)
            cigars.append(cigar)
            
        return np.array(scores), cigars
    
    def _banded_alignment(self, read: str, ref_start: int) -> Tuple[int, str]:
        """
        Perform banded alignment allowing indels
        Returns alignment score and CIGAR string
        """
        read_len = len(read)
        ref_end = min(ref_start + read_len + self.max_indel_size, self.ref_len)
        ref_segment = self.reference[ref_start:ref_end]
        
        # Simple banded DP for demonstration
        # In production, use more efficient algorithms
        m, n = len(read), len(ref_segment)
        dp = np.full((m + 1, n + 1), float('inf'))
        dp[0, :] = np.arange(n + 1)
        dp[:, 0] = np.arange(m + 1)
        
        for i in range(1, m + 1):
            for j in range(max(1, i - self.max_indel_size), 
                           min(n + 1, i + self.max_indel_size + 1)):
                if j < n + 1:
                    match = dp[i-1, j-1] + (0 if read[i-1] == ref_segment[j-1] else 1)
                    insert = dp[i-1, j] + 1
                    delete = dp[i, j-1] + 1 if j > 0 else float('inf')
                    dp[i, j] = min(match, insert, delete)
        
        # Traceback for CIGAR
        cigar = self._traceback_cigar(dp, read, ref_segment)
        score = dp[m, n]
        
        return int(score), cigar
    
    def _traceback_cigar(self, dp: np.ndarray, read: str, ref: str) -> str:
        """Generate CIGAR string from DP matrix"""
        m, n = len(read), len(ref)
        i, j = m, n
        cigar_ops = []
        
        while i > 0 and j > 0:
            if i > 0 and j > 0 and dp[i, j] == dp[i-1, j-1] + (0 if read[i-1] == ref[j-1] else 1):
                cigar_ops.append('M')
                i -= 1
                j -= 1
            elif i > 0 and dp[i, j] == dp[i-1, j] + 1:
                cigar_ops.append('I')
                i -= 1
            else:
                cigar_ops.append('D')
                j -= 1
                
        # Compress CIGAR
        cigar_ops.reverse()
        cigar = []
        if cigar_ops:
            count = 1
            current = cigar_ops[0]
            for op in cigar_ops[1:]:
                if op == current:
                    count += 1
                else:
                    cigar.append(f'{count}{current}')
                    current = op
                    count = 1
            cigar.append(f'{count}{current}')
            
        return ''.join(cigar)
    
    def _find_splice_alignments(self, read: str, seed_hits: Dict[int, List[int]]) -> List[AlignmentResult]:
        """Find splice-aware alignments for RNA-seq reads"""
        splice_alignments = []
        read_len = len(read)
        
        # Look for split alignments
        for pos1, hits1 in seed_hits.items():
            if pos1 > read_len // 2:
                continue
                
            for pos2, hits2 in seed_hits.items():
                if pos2 <= pos1 + self.seed_len:
                    continue
                    
                # Check if this could be a splice junction
                gap_in_read = pos2 - (pos1 + self.seed_len)
                
                for hit1 in hits1:
                    for hit2 in hits2:
                        gap_in_ref = hit2 - (hit1 + self.seed_len)
                        
                        # Large gap in reference suggests intron
                        if 50 <= gap_in_ref <= 100000 and gap_in_read < 10:
                            # Verify splice motif
                            donor_pos = hit1 + pos1 + self.seed_len
                            acceptor_pos = hit2 + pos2 - 2
                            
                            if donor_pos in self.splice_sites:
                                for acceptor, motif in self.splice_sites[donor_pos]:
                                    if abs(acceptor - acceptor_pos) < 5:
                                        # Found valid splice junction
                                        cigar = f'{pos1 + self.seed_len}M{gap_in_ref}N{read_len - pos2}M'
                                        result = AlignmentResult(
                                            query_id='',
                                            ref_pos=hit1,
                                            mapq=60,
                                            cigar=cigar,
                                            mismatches=0,
                                            is_splice_junction=True
                                        )
                                        splice_alignments.append(result)
                                        
        return splice_alignments
    
    def align_single(self, read_seq: str, read_id: str = '') -> List[AlignmentResult]:
        """Align a single read with full feature support"""
        results = []
        
        # Extract seeds at multiple positions
        seed_hits = defaultdict(list)
        for offset in range(0, len(read_seq) - self.seed_len + 1, self.seed_step):
            seed = read_seq[offset:offset + self.seed_len]
            if seed in self.seed_index:
                for hit in self.seed_index[seed]:
                    seed_hits[offset].append(hit - offset)
        
        # Check for splice junctions if RNA-seq mode
        if self.splice_motifs:
            splice_results = self._find_splice_alignments(read_seq, seed_hits)
            results.extend(splice_results)
        
        # Collect unique candidate positions
        candidates = set()
        for hits in seed_hits.values():
            candidates.update(hits)
        candidates = sorted(candidates)
        
        # Score candidates with indel support
        scores, cigars = self._vectorized_scoring_with_indels(read_seq, candidates)
        
        if len(scores) > 0:
            best_idx = int(np.argmin(scores))  # Cast to int
            best_score = scores[best_idx]
            
            if best_score <= self.max_mismatches:
                # Calculate mapping quality
                mapq = self._calculate_mapq(scores, best_idx)
                
                result = AlignmentResult(
                    query_id=read_id,
                    ref_pos=candidates[best_idx],
                    mapq=mapq,
                    cigar=cigars[best_idx],
                    mismatches=int(best_score)
                )
                results.append(result)
                
        return results
    
    def align_paired(self, read1: str, read2: str, 
                    read_id: str = '',
                    min_insert: int = 100,
                    max_insert: int = 1000) -> Tuple[Optional[AlignmentResult], Optional[AlignmentResult]]:
        """
        Align paired-end reads considering insert size constraints
        """
        # Align both reads independently first
        results1 = self.align_single(read1, read_id + '/1')
        results2 = self.align_single(read2, read_id + '/2')
        
        if not results1 or not results2:
            # Try to rescue unmapped mate
            if results1 and not results2:
                results2 = self._rescue_mate(read2, results1[0].ref_pos, max_insert)
            elif results2 and not results1:
                results1 = self._rescue_mate(read1, results2[0].ref_pos, max_insert)
            else:
                return None, None
        
        # Find best pair considering insert size
        best_pair = None
        best_score = float('inf')
        
        for r1 in results1:
            for r2 in results2:
                # Check orientation and insert size
                if r1.ref_pos < r2.ref_pos:
                    insert_size = r2.ref_pos + len(read2) - r1.ref_pos
                    if min_insert <= insert_size <= max_insert:
                        pair_score = r1.mismatches + r2.mismatches
                        if pair_score < best_score:
                            best_score = pair_score
                            r1.is_paired = True
                            r2.is_paired = True
                            r1.mate_pos = r2.ref_pos
                            r2.mate_pos = r1.ref_pos
                            r1.insert_size = insert_size
                            r2.insert_size = insert_size
                            best_pair = (r1, r2)
                            
        return best_pair if best_pair else (results1[0] if results1 else None,
                                           results2[0] if results2 else None)
    
    def _rescue_mate(self, unmapped_read: str, mate_pos: int, max_insert: int) -> List[AlignmentResult]:
        """Try to align unmapped read near its mate's position"""
        results = []
        search_start = max(0, mate_pos - max_insert)
        search_end = min(self.ref_len, mate_pos + max_insert)
        
        # More intensive search in mate's vicinity
        candidates = list(range(search_start, search_end, 50))
        scores, cigars = self._vectorized_scoring_with_indels(unmapped_read, candidates)
        
        if len(scores) > 0:
            best_idx = int(np.argmin(scores))  # Cast to int
            if scores[best_idx] <= self.max_mismatches * 1.5:  # More lenient for rescue
                result = AlignmentResult(
                    query_id='',
                    ref_pos=candidates[best_idx],
                    mapq=0,  # Lower quality for rescued
                    cigar=cigars[best_idx],
                    mismatches=int(scores[best_idx])
                )
                results.append(result)
                
        return results
    
    def _calculate_mapq(self, scores: np.ndarray, best_idx: int) -> int:
        """Calculate mapping quality score"""
        if len(scores) == 1:
            return 60
            
        sorted_scores = np.sort(scores)
        best_score = sorted_scores[0]
        second_best = sorted_scores[1] if len(sorted_scores) > 1 else float('inf')
        
        # Simple MAPQ calculation
        if second_best == float('inf'):
            mapq = 60
        else:
            mapq = min(60, int(6 * (second_best - best_score)))
            
        return max(0, mapq)

def align_reads(reference: str, 
                reads: List[Tuple[str, str]],
                paired_reads: Optional[List[Tuple[str, str, str]]] = None,
                **kwargs) -> List[AlignmentResult]:
    """
    Main alignment function with full feature support
    
    Args:
        reference: Reference sequence
        reads: List of (read_seq, read_id) tuples for single-end
        paired_reads: List of (read1, read2, read_id) for paired-end
        **kwargs: Additional parameters for VecMapProduction
        
    Returns:
        List of AlignmentResult objects
    """
    aligner = VecMapProduction(reference, **kwargs)
    results = []
    
    # Process single-end reads
    if reads:
        for read_seq, read_id in reads:
            alignments = aligner.align_single(read_seq, read_id)
            results.extend(alignments)
    
    # Process paired-end reads
    if paired_reads:
        for read1, read2, read_id in paired_reads:
            r1, r2 = aligner.align_paired(read1, read2, read_id)
            if r1:
                results.append(r1)
            if r2:
                results.append(r2)
                
    return results

# Backward compatibility wrapper
def vecmap(ref, reads, read_len, **kwargs):
    """Wrapper for backward compatibility with original vecmap interface"""
    single_reads = [(read, f'read_{i}') for i, (read, _) in enumerate(reads)]
    results = align_reads(ref, single_reads, **kwargs)
    
    # Convert to old format
    mappings = []
    for i, (read, true_pos) in enumerate(reads):
        if i < len(results):
            r = results[i]
            mappings.append((r.ref_pos, r.mismatches, true_pos))
        else:
            mappings.append((-1, -1, true_pos))
            
    return mappings 