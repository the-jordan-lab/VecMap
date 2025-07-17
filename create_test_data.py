#!/usr/bin/env python3
"""Create test data for VecMap testing"""

import random

# Create a small transcriptome-like reference
print("Creating test reference...")

# Generate 5 "transcripts" of varying lengths
transcripts = []
for i in range(5):
    length = random.randint(300, 800)
    seq = ''.join(random.choice('ACGT') for _ in range(length))
    transcripts.append((f"transcript_{i}", seq))

# Write reference FASTA
with open('test_reference.fa', 'w') as f:
    for name, seq in transcripts:
        f.write(f">{name}\n")
        # Write in 80-char lines
        for j in range(0, len(seq), 80):
            f.write(seq[j:j+80] + "\n")

# Concatenate all transcripts for read generation
full_ref = ''.join(t[1] for t in transcripts)
print(f"Reference length: {len(full_ref)} bp")

# Generate 1000 reads from the reference
print("Creating test reads...")
reads = []
read_len = 100

for i in range(1000):
    # Random position
    pos = random.randint(0, len(full_ref) - read_len)
    read_seq = full_ref[pos:pos + read_len]
    
    # Add some errors (1% error rate)
    read_list = list(read_seq)
    for j in range(len(read_list)):
        if random.random() < 0.01:
            read_list[j] = random.choice([b for b in 'ACGT' if b != read_list[j]])
    
    read_seq = ''.join(read_list)
    reads.append((f"read_{i}", read_seq))

# Write FASTQ
with open('test_reads.fq', 'w') as f:
    for name, seq in reads:
        f.write(f"@{name}\n")
        f.write(f"{seq}\n")
        f.write("+\n")
        f.write("I" * len(seq) + "\n")  # High quality scores

print(f"Created {len(reads)} reads")
print("\nTest data created:")
print("  Reference: test_reference.fa")
print("  Reads: test_reads.fq") 