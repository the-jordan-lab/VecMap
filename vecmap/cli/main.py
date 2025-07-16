#!/usr/bin/env python3
"""
VecMap Command-Line Interface
============================

Simple CLI for running VecMap alignment.
"""

import argparse
import sys
import time
from typing import List, Tuple, Optional
from ..core.mapper import vecmap


def parse_fasta(filename: str) -> Tuple[str, str]:
    """Parse a simple FASTA file (single sequence)."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    header = lines[0].strip()
    sequence = ''.join(line.strip() for line in lines[1:])
    
    return header, sequence


def parse_fastq(filename: str, max_reads: Optional[int] = None) -> List[Tuple[str, str]]:
    """Parse FASTQ file and return (sequence, read_id) tuples."""
    reads = []
    
    with open(filename, 'r') as f:
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
    
    return reads


def main():
    """Main entry point for VecMap CLI."""
    parser = argparse.ArgumentParser(
        description="VecMap: Vectorized sequence alignment for exact matching",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Align RNA-seq reads to transcriptome
  vecmap -r transcriptome.fa -q reads.fq -o alignments.txt
  
  # Process subset of reads
  vecmap -r reference.fa -q reads.fq -n 10000 -o test.txt
  
  # Specify k-mer size
  vecmap -r reference.fa -q reads.fq -k 20 -o alignments.txt
        """
    )
    
    parser.add_argument('-r', '--reference', required=True,
                        help='Reference sequence (FASTA format)')
    parser.add_argument('-q', '--query', required=True,
                        help='Query sequences (FASTQ format)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output alignment file')
    parser.add_argument('-k', '--kmer', type=int, default=20,
                        help='Seed length for indexing (default: 20)')
    parser.add_argument('-n', '--max-reads', type=int, default=None,
                        help='Maximum number of reads to process')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')
    
    args = parser.parse_args()
    
    # Load reference
    if args.verbose:
        print(f"Loading reference from {args.reference}...")
    
    ref_header, reference = parse_fasta(args.reference)
    
    if args.verbose:
        print(f"Reference loaded: {len(reference):,} bp")
    
    # Load reads
    if args.verbose:
        print(f"Loading reads from {args.query}...")
    
    reads = parse_fastq(args.query, args.max_reads)

    if args.verbose:
        print(f"Loaded {len(reads):,} reads")

    # Determine read length
    read_len = len(reads[0][0]) if reads else 0
    if any(len(seq) != read_len for seq, _ in reads):
        print("Warning: reads have varying lengths; using first read length")

    # Run alignment
    if args.verbose:
        print(f"Running VecMap alignment (seed_len={args.kmer}, read_len={read_len})...")

    start_time = time.time()
    alignments = vecmap(reference, reads, read_len, seed_len=args.kmer)
    elapsed_time = time.time() - start_time
    
    if args.verbose:
        print(f"Alignment complete in {elapsed_time:.2f} seconds")
        print(f"Speed: {len(reads)/elapsed_time:,.0f} reads/second")
    
    # Write output
    with open(args.output, 'w') as f:
        f.write(f"# VecMap alignment results\n")
        f.write(f"# Reference: {ref_header}\n")
        f.write(f"# Query: {args.query}\n")
        f.write(f"# Seed length: {args.kmer}\n")
        f.write(f"# Read length: {read_len}\n")
        f.write(f"# Reads processed: {len(reads)}\n")
        f.write(f"# Time: {elapsed_time:.2f} seconds\n")
        f.write(f"# Speed: {len(reads)/elapsed_time:,.0f} reads/second\n")
        f.write("#\n")
        f.write("# Format: read_id\tref_position\tmismatches\n")
        
        for pos, mismatches, read_id in alignments:
            f.write(f"{read_id}\t{pos}\t{mismatches}\n")
    
    if args.verbose:
        print(f"Results written to {args.output}")
        print(f"Aligned {len(alignments)} reads ({len(alignments)/len(reads)*100:.1f}%)")


if __name__ == "__main__":
    main() 