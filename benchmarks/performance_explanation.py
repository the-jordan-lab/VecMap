#!/usr/bin/env python3
"""
Demonstrate performance differences between reference types

This script shows why VecMap performs differently on continuous vs
transcriptome-style references, explaining the benchmark results.
"""

import time
from vecmap import vecmap, generate_reference, generate_reads
from test_geo_quick import generate_transcriptome, simulate_rnaseq_reads

print("VecMap Performance: Reference Type Comparison")
print("="*60)
print("\nThis demonstrates why benchmark results differ from simple tests")
print()

# Test 1: Continuous reference (like generate_reference)
print("1. CONTINUOUS REFERENCE (single sequence)")
print("-"*40)
ref_length = 350000  # ~350kb
num_reads = 25000

ref = generate_reference(ref_length)
reads = generate_reads(ref, num_reads, 100)

# Warm-up
_ = vecmap(ref, reads, 100)

# Benchmark
times = []
for i in range(3):
    start = time.time()
    _ = vecmap(ref, reads, 100)
    times.append(time.time() - start)

avg_time = sum(times) / len(times)
speed = num_reads / avg_time

print(f"Reference length: {len(ref):,} bp")
print(f"Number of reads: {num_reads:,}")
print(f"Average time: {avg_time:.3f} seconds")
print(f"Speed: {speed:,.0f} reads/second")

# Test 2: Transcriptome-style reference (multiple sequences)
print("\n2. TRANSCRIPTOME REFERENCE (concatenated sequences)")
print("-"*40)

ref_seq, _, pos_map = generate_transcriptome(200)  # 200 transcripts
reads = simulate_rnaseq_reads(ref_seq, pos_map, num_reads)

# Warm-up
_ = vecmap(ref_seq, reads, 100)

# Benchmark
times = []
for i in range(3):
    start = time.time()
    _ = vecmap(ref_seq, reads, 100)
    times.append(time.time() - start)

avg_time = sum(times) / len(times)
speed = num_reads / avg_time

print(f"Reference length: {len(ref_seq):,} bp")
print(f"Number of transcripts: 200")
print(f"Number of reads: {num_reads:,}")
print(f"Average time: {avg_time:.3f} seconds")
print(f"Speed: {speed:,.0f} reads/second")

print("\n" + "="*60)
print("EXPLANATION:")
print("="*60)
print("""
The performance difference is due to how k-mer indices work:

1. Continuous references have highly repetitive k-mers (especially
   when generated with a repeating pattern), leading to many candidate
   positions per k-mer and slower mapping.

2. Transcriptome references have more diverse sequences, resulting in
   fewer candidate positions per k-mer and faster mapping.

The benchmarks use transcriptome-style data, which is more representative
of real RNA-seq applications where VecMap is designed to excel.
""") 