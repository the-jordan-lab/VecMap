#!/usr/bin/env python3
"""
Quick test of VecMap CRISPR guide detection performance
"""

import time
import sys
from pathlib import Path
import numpy as np

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent.parent))
from vecmap.applications import CRISPRGuideDetector

def create_test_guides(num_guides=1000):
    """Create test guide library."""
    guides = {}
    nucleotides = ['A', 'C', 'G', 'T']
    
    for i in range(num_guides):
        guide_seq = 'G' + ''.join(np.random.choice(nucleotides, 19))
        guides[f"guide_{i}"] = guide_seq
    
    return guides

def generate_reads(guides, num_reads=100000):
    """Generate test reads."""
    reads = []
    guide_list = list(guides.items())
    
    for i in range(num_reads):
        guide_name, guide_seq = guide_list[i % len(guide_list)]
        # Add context
        full_seq = "ACCG" + guide_seq + "GTTT"
        reads.append((full_seq, f"read_{i}"))
    
    return reads

def main():
    print("VecMap CRISPR Guide Detection Benchmark")
    print("="*40)
    
    # Test parameters
    guide_counts = [100, 500, 1000, 5000, 10000]
    read_count = 1000000  # 1M reads
    
    for num_guides in guide_counts:
        print(f"\nTesting with {num_guides} guides...")
        
        # Create guides
        guides = create_test_guides(num_guides)
        
        # Generate reads
        reads = generate_reads(guides, read_count)
        
        # Initialize detector
        detector = CRISPRGuideDetector(guides)
        
        # Benchmark
        start_time = time.time()
        results = detector.detect_guides_with_context(reads, 
                                                      upstream_context="ACCG",
                                                      downstream_context="GTTT")
        counts = detector.summarize_detection(results)
        end_time = time.time()
        
        # Results
        elapsed = end_time - start_time
        reads_per_sec = read_count / elapsed
        
        print(f"  Time: {elapsed:.2f} seconds")
        print(f"  Speed: {reads_per_sec:,.0f} reads/second")
        print(f"  Guides detected: {len(counts)}")
        
    print("\n" + "="*40)
    print("Benchmark complete!")

if __name__ == "__main__":
    main() 