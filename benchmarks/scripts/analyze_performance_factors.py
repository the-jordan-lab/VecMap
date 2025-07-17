#!/usr/bin/env python3
"""Analyze what makes VecMap fast or slow"""

from vecmap.applications.crispr import CRISPRGuideDetector
import time
import random
import numpy as np

print("VecMap Performance Factor Analysis")
print("="*60)

# Test different factors that might affect performance

# Factor 1: Number of guides in library
print("\n1. GUIDE LIBRARY SIZE IMPACT")
print("-"*40)

for num_guides in [10, 50, 100, 500, 1000, 5000]:
    guides = {}
    for i in range(num_guides):
        guides[f"guide_{i}"] = ''.join(random.choice('ACGT') for _ in range(20))
    
    # Generate 10k reads
    reads = []
    for i in range(10000):
        guide_name = random.choice(list(guides.keys()))
        full_seq = "ACCG" + guides[guide_name] + "GTTT"
        reads.append((full_seq, f"read_{i}"))
    
    detector = CRISPRGuideDetector(guides)
    
    # Time it
    start = time.time()
    _ = detector.detect_guides_with_context(reads, "ACCG", "GTTT")
    elapsed = time.time() - start
    
    speed = len(reads) / elapsed
    print(f"  {num_guides:5d} guides: {speed:8,.0f} reads/sec")

# Factor 2: Read count impact
print("\n2. READ COUNT IMPACT (100 guides)")
print("-"*40)

guides = {}
for i in range(100):
    guides[f"guide_{i}"] = ''.join(random.choice('ACGT') for _ in range(20))

detector = CRISPRGuideDetector(guides)

for num_reads in [1000, 5000, 10000, 50000, 100000]:
    reads = []
    for i in range(num_reads):
        guide_name = random.choice(list(guides.keys()))
        full_seq = "ACCG" + guides[guide_name] + "GTTT"
        reads.append((full_seq, f"read_{i}"))
    
    start = time.time()
    _ = detector.detect_guides_with_context(reads, "ACCG", "GTTT")
    elapsed = time.time() - start
    
    speed = len(reads) / elapsed
    print(f"  {num_reads:6d} reads: {speed:8,.0f} reads/sec")

# Factor 3: K-mer diversity (what fraction of guides share k-mers)
print("\n3. K-MER DIVERSITY IMPACT")
print("-"*40)

def create_guides_with_overlap(num_guides, overlap_fraction):
    """Create guides where some share k-mers"""
    guides = {}
    
    if overlap_fraction == 0:
        # All unique
        for i in range(num_guides):
            guides[f"guide_{i}"] = ''.join(random.choice('ACGT') for _ in range(20))
    else:
        # Create some base sequences
        num_bases = int(num_guides * (1 - overlap_fraction))
        base_seqs = []
        for i in range(num_bases):
            base_seqs.append(''.join(random.choice('ACGT') for _ in range(20)))
        
        # Create guides with overlaps
        for i in range(num_guides):
            if i < num_bases:
                guides[f"guide_{i}"] = base_seqs[i]
            else:
                # Create variant of existing guide
                base_idx = random.randint(0, num_bases - 1)
                seq = list(base_seqs[base_idx])
                # Change 1-2 positions
                for _ in range(random.randint(1, 2)):
                    pos = random.randint(0, 19)
                    seq[pos] = random.choice('ACGT')
                guides[f"guide_{i}"] = ''.join(seq)
    
    return guides

for overlap in [0.0, 0.25, 0.5, 0.75]:
    guides = create_guides_with_overlap(500, overlap)
    
    # Generate reads
    reads = []
    for i in range(10000):
        guide_name = random.choice(list(guides.keys()))
        full_seq = "ACCG" + guides[guide_name] + "GTTT"
        reads.append((full_seq, f"read_{i}"))
    
    detector = CRISPRGuideDetector(guides)
    
    start = time.time()
    _ = detector.detect_guides_with_context(reads, "ACCG", "GTTT")
    elapsed = time.time() - start
    
    speed = len(reads) / elapsed
    print(f"  {overlap*100:3.0f}% overlap: {speed:8,.0f} reads/sec")

print("\n" + "="*60)
print("KEY INSIGHTS:")
print("- Smaller guide libraries = faster performance")
print("- Read count has minimal impact on speed per read")
print("- K-mer diversity matters for general VecMap but less for CRISPR")
print("  (since we search for exact 20bp guides)") 