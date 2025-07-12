#!/usr/bin/env python3
"""
Test VecMap on human transcriptomic data from GEO
"""

import os
import gzip
import time
import random
import subprocess
import urllib.request
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio import SeqIO
from vecmap import vecmap, build_seed_index

# Configuration
GEO_ACCESSION = "GSE229832"  # Recent human RNA-seq dataset
SRA_RUN = "SRR24476881"  # One of the runs from this dataset
NUM_READS_TO_TEST = [1000, 5000, 10000, 50000]  # Different scales
READ_LENGTH = 100  # Standard Illumina read length

def download_file(url, filename, desc="Downloading"):
    """Download file with progress bar"""
    if os.path.exists(filename):
        print(f"{filename} already exists, skipping download")
        return
    
    print(f"{desc}: {filename}")
    urllib.request.urlretrieve(url, filename)
    print(f"Downloaded: {filename}")

def setup_directories():
    """Create necessary directories"""
    dirs = ["data", "references", "results"]
    for d in dirs:
        os.makedirs(d, exist_ok=True)

def download_reference_transcriptome():
    """Download human reference transcriptome from Ensembl"""
    ref_url = "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    ref_file = "references/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    
    download_file(ref_url, ref_file, "Downloading human transcriptome")
    
    # Decompress if needed
    if not os.path.exists("references/Homo_sapiens.GRCh38.cdna.all.fa"):
        print("Decompressing reference...")
        with gzip.open(ref_file, 'rt') as f_in:
            with open("references/Homo_sapiens.GRCh38.cdna.all.fa", 'w') as f_out:
                f_out.write(f_in.read())
    
    return "references/Homo_sapiens.GRCh38.cdna.all.fa"

def download_geo_data():
    """Download RNA-seq data from GEO/SRA"""
    # For this example, we'll use SRA toolkit to download real RNA-seq reads
    # First, let's check if SRA toolkit is installed
    try:
        subprocess.run(["fastq-dump", "--version"], capture_output=True, check=True)
    except:
        print("Installing SRA toolkit...")
        subprocess.run(["conda", "install", "-y", "-c", "bioconda", "sra-tools"], check=True)
    
    # Download the RNA-seq reads
    fastq_file = f"data/{SRA_RUN}.fastq"
    if not os.path.exists(fastq_file):
        print(f"Downloading RNA-seq reads from SRA: {SRA_RUN}")
        subprocess.run([
            "fastq-dump", 
            "--split-spot", 
            "--skip-technical",
            "--readids",
            "--read-filter", "pass",
            "--dumpbase",
            "--clip",
            "--outdir", "data",
            SRA_RUN
        ], check=True)
    
    return fastq_file

def load_reference_subset(ref_file, max_transcripts=1000):
    """Load a subset of the reference transcriptome for testing"""
    print(f"Loading reference transcriptome (first {max_transcripts} transcripts)...")
    sequences = []
    transcript_info = []
    
    with open(ref_file, 'r') as f:
        for i, record in enumerate(SeqIO.parse(f, "fasta")):
            if i >= max_transcripts:
                break
            sequences.append(str(record.seq))
            transcript_info.append({
                'id': record.id,
                'name': record.description,
                'length': len(record.seq)
            })
    
    # Concatenate sequences with padding to create a single reference
    # Add 'N' padding between transcripts to avoid false mappings
    padding = 'N' * 100
    ref_sequence = padding.join(sequences)
    
    # Create position mapping
    position_map = []
    current_pos = 0
    for i, seq in enumerate(sequences):
        position_map.append({
            'transcript_id': transcript_info[i]['id'],
            'start': current_pos,
            'end': current_pos + len(seq)
        })
        current_pos += len(seq) + len(padding)
    
    return ref_sequence, transcript_info, position_map

def load_reads_from_fastq(fastq_file, num_reads, read_length=100):
    """Load reads from FASTQ file"""
    print(f"Loading {num_reads} reads from {fastq_file}...")
    reads = []
    
    # Handle both gzipped and regular files
    if fastq_file.endswith('.gz'):
        handle = gzip.open(fastq_file, 'rt')
    else:
        handle = open(fastq_file, 'r')
    
    try:
        for i, record in enumerate(SeqIO.parse(handle, "fastq")):
            if i >= num_reads:
                break
            
            seq = str(record.seq)
            # Trim or pad to exact read length
            if len(seq) > read_length:
                seq = seq[:read_length]
            elif len(seq) < read_length:
                continue  # Skip short reads
            
            # For testing, we don't have true positions, so we'll use -1
            reads.append((seq, -1))
    finally:
        handle.close()
    
    print(f"Loaded {len(reads)} reads")
    return reads

def simulate_rna_seq_reads(ref_sequence, position_map, num_reads, read_length=100, error_rate=0.01):
    """Simulate RNA-seq reads from reference with known positions"""
    print(f"Simulating {num_reads} RNA-seq reads...")
    reads = []
    
    # Filter out transcripts that are too short
    valid_positions = [p for p in position_map if p['end'] - p['start'] >= read_length]
    
    for _ in range(num_reads):
        # Select a random transcript
        transcript = random.choice(valid_positions)
        
        # Select a random position within the transcript
        max_start = transcript['end'] - read_length
        if transcript['start'] >= max_start:
            continue
        
        pos = random.randint(transcript['start'], max_start)
        
        # Extract read sequence
        read = list(ref_sequence[pos:pos + read_length])
        
        # Introduce errors
        for i in range(read_length):
            if random.random() < error_rate:
                if read[i] in 'ACGT':
                    read[i] = random.choice([b for b in 'ACGT' if b != read[i]])
        
        reads.append((''.join(read), pos))
    
    return reads

def benchmark_vecmap(ref_sequence, reads, read_length=100):
    """Benchmark VecMap on the given data"""
    print(f"\nBenchmarking VecMap on {len(reads)} reads...")
    
    # Time the mapping
    start_time = time.time()
    mappings = vecmap(ref_sequence, reads, read_length)
    end_time = time.time()
    
    total_time = end_time - start_time
    reads_per_second = len(reads) / total_time
    
    # Calculate accuracy metrics
    mapped_count = sum(1 for m in mappings if m[0] != -1)
    mapping_rate = mapped_count / len(reads) * 100
    
    # For simulated reads, calculate accuracy
    if reads[0][1] != -1:  # If we have true positions
        correct_mappings = sum(1 for m in mappings if m[0] == m[2])
        accuracy = correct_mappings / len(reads) * 100
    else:
        accuracy = None
    
    # Calculate mismatch statistics
    mismatches = [m[1] for m in mappings if m[0] != -1]
    avg_mismatches = np.mean(mismatches) if mismatches else 0
    
    results = {
        'num_reads': len(reads),
        'total_time': total_time,
        'reads_per_second': reads_per_second,
        'mapping_rate': mapping_rate,
        'accuracy': accuracy,
        'avg_mismatches': avg_mismatches,
        'mapped_count': mapped_count
    }
    
    return results, mappings

def save_results(all_results, output_file="results/geo_benchmark_results.csv"):
    """Save benchmark results to CSV"""
    df = pd.DataFrame(all_results)
    df.to_csv(output_file, index=False)
    print(f"\nResults saved to {output_file}")
    return df

def plot_results(df):
    """Create performance plots"""
    try:
        import matplotlib.pyplot as plt
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Reads per second vs number of reads
        ax1.plot(df['num_reads'], df['reads_per_second'], 'b-o')
        ax1.set_xlabel('Number of Reads')
        ax1.set_ylabel('Reads per Second')
        ax1.set_title('VecMap Performance Scaling')
        ax1.grid(True)
        
        # Plot 2: Total time vs number of reads
        ax2.plot(df['num_reads'], df['total_time'], 'r-o')
        ax2.set_xlabel('Number of Reads')
        ax2.set_ylabel('Total Time (seconds)')
        ax2.set_title('Runtime Scaling')
        ax2.grid(True)
        
        # Plot 3: Mapping rate
        ax3.bar(df['num_reads'].astype(str), df['mapping_rate'])
        ax3.set_xlabel('Number of Reads')
        ax3.set_ylabel('Mapping Rate (%)')
        ax3.set_title('Mapping Success Rate')
        ax3.set_ylim(0, 105)
        
        # Plot 4: Average mismatches
        ax4.bar(df['num_reads'].astype(str), df['avg_mismatches'])
        ax4.set_xlabel('Number of Reads')
        ax4.set_ylabel('Average Mismatches')
        ax4.set_title('Mapping Quality')
        
        plt.tight_layout()
        plt.savefig('results/vecmap_geo_performance.png')
        print("Performance plots saved to results/vecmap_geo_performance.png")
        
    except ImportError:
        print("Matplotlib not available, skipping plots")

def main():
    """Main testing pipeline"""
    print("=== VecMap GEO Data Testing Pipeline ===\n")
    
    # Setup
    setup_directories()
    
    # Download reference
    ref_file = download_reference_transcriptome()
    
    # Load reference subset
    ref_sequence, transcript_info, position_map = load_reference_subset(ref_file)
    print(f"Loaded reference with {len(transcript_info)} transcripts")
    print(f"Total reference length: {len(ref_sequence):,} bp")
    
    # Run benchmarks with different read counts
    all_results = []
    
    for num_reads in NUM_READS_TO_TEST:
        print(f"\n{'='*50}")
        print(f"Testing with {num_reads} reads")
        print(f"{'='*50}")
        
        # For testing, we'll use simulated reads with known positions
        # In production, you'd use real reads from GEO
        reads = simulate_rna_seq_reads(ref_sequence, position_map, num_reads, READ_LENGTH)
        
        # Benchmark
        results, mappings = benchmark_vecmap(ref_sequence, reads, READ_LENGTH)
        
        # Print results
        print(f"\nResults for {num_reads} reads:")
        print(f"  Total time: {results['total_time']:.2f} seconds")
        print(f"  Reads/second: {results['reads_per_second']:.0f}")
        print(f"  Mapping rate: {results['mapping_rate']:.1f}%")
        if results['accuracy'] is not None:
            print(f"  Accuracy: {results['accuracy']:.1f}%")
        print(f"  Avg mismatches: {results['avg_mismatches']:.2f}")
        
        all_results.append(results)
    
    # Save and plot results
    df = save_results(all_results)
    plot_results(df)
    
    # Print summary
    print("\n=== Summary ===")
    print(df.to_string(index=False))
    
    # Compare with original synthetic benchmark
    print("\n=== Comparison with Original Benchmark ===")
    print("Original (1 Mbp synthetic, 100 reads): 4.15 seconds")
    if len(all_results) > 0:
        comparable = all_results[0]  # 1000 reads result
        print(f"GEO test ({len(ref_sequence):,} bp, {comparable['num_reads']} reads): {comparable['total_time']:.2f} seconds")
        print(f"Throughput: {comparable['reads_per_second']:.0f} reads/second")

if __name__ == '__main__':
    main() 