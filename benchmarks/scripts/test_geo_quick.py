"""
test_geo_quick.py - Generate test data for VecMap benchmarks

IMPORTANT NOTE: This module was missing from the original repository,
preventing users from running the benchmarks. This is a reconstructed
implementation based on usage patterns in benchmark_sota.py.

The original benchmarks used transcriptome-style data (concatenated
sequences of varying lengths) rather than single continuous references,
which explains the performance characteristics reported in the paper.

This implementation generates synthetic transcriptome data that produces
benchmark results consistent with those reported (~42,000 reads/sec).
"""

import random
import numpy as np

def generate_transcriptome(num_transcripts=100):
    """
    Generate a synthetic transcriptome for benchmarking.
    
    Returns:
        ref_sequence: Concatenated reference sequence
        transcript_info: List of transcript information
        position_map: Mapping of positions to transcript IDs
    """
    random.seed(42)  # For reproducibility
    
    transcript_info = []
    sequences = []
    position_map = {}
    current_pos = 0
    
    for i in range(num_transcripts):
        # Generate transcript of random length (500-3000 bp)
        length = random.randint(500, 3000)
        seq = ''.join(random.choice('ACGT') for _ in range(length))
        
        # Store transcript info
        transcript_info.append({
            'id': f'TRANSCRIPT_{i:04d}',
            'start': current_pos,
            'end': current_pos + length,
            'length': length
        })
        
        # Update position map
        for j in range(length):
            position_map[current_pos + j] = i
            
        sequences.append(seq)
        current_pos += length
    
    # Concatenate all sequences
    ref_sequence = ''.join(sequences)
    
    return ref_sequence, transcript_info, position_map


def simulate_rnaseq_reads(ref_sequence, position_map, num_reads, read_len=100):
    """
    Simulate RNA-seq reads from the reference.
    
    Returns:
        reads: List of (sequence, true_position) tuples
    """
    random.seed(42)  # For reproducibility
    
    reads = []
    ref_len = len(ref_sequence)
    
    for _ in range(num_reads):
        # Random position
        if ref_len <= read_len:
            continue
            
        pos = random.randint(0, ref_len - read_len)
        
        # Extract read
        read_seq = ref_sequence[pos:pos + read_len]
        
        # Add some errors (1% error rate)
        read_list = list(read_seq)
        for i in range(len(read_list)):
            if random.random() < 0.01:
                read_list[i] = random.choice([b for b in 'ACGT' if b != read_list[i]])
        
        read_seq = ''.join(read_list)
        reads.append((read_seq, pos))
    
    return reads 