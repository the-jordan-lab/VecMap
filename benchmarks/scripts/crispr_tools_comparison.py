#!/usr/bin/env python3
"""
Ultimate CRISPR Tools Benchmark: VecMap vs CRISPResso2 vs MAGeCK
================================================================

This benchmark compares guide detection performance across:
- VecMap: Our vectorized exact matching approach
- CRISPResso2: Standard tool for CRISPR analysis
- MAGeCK: Gold standard for CRISPR screen analysis

We focus on the guide counting/detection aspect for fair comparison.
"""

import os
import time
import subprocess
import tempfile
import pandas as pd
import numpy as np
from pathlib import Path
import gzip
import shutil
from typing import Dict, List, Tuple
import argparse
import json

# Import VecMap
import sys
sys.path.append(str(Path(__file__).parent.parent.parent))
from vecmap.applications import CRISPRGuideDetector


class CRISPRBenchmark:
    """Benchmark framework for CRISPR guide detection tools."""
    
    def __init__(self, guides: Dict[str, str], output_dir: str):
        self.guides = guides
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.results = []
        
    def generate_test_data(self, num_reads: int, guide_distribution: str = "uniform") -> Tuple[str, Dict[str, int]]:
        """Generate synthetic CRISPR screening data with known ground truth."""
        print(f"\nGenerating {num_reads:,} test reads...")
        
        # Create ground truth counts
        ground_truth = {}
        guide_names = list(self.guides.keys())
        
        if guide_distribution == "uniform":
            # Uniform distribution
            counts_per_guide = num_reads // len(guide_names)
            ground_truth = {name: counts_per_guide for name in guide_names}
            # Add remaining reads to random guides
            remaining = num_reads - (counts_per_guide * len(guide_names))
            for i in range(remaining):
                ground_truth[guide_names[i % len(guide_names)]] += 1
                
        elif guide_distribution == "power_law":
            # Power law distribution (more realistic)
            # Top 20% of guides get 80% of reads
            top_guides = int(0.2 * len(guide_names))
            top_reads = int(0.8 * num_reads)
            bottom_reads = num_reads - top_reads
            
            # Distribute top reads
            for i in range(top_guides):
                ground_truth[guide_names[i]] = top_reads // top_guides
                
            # Distribute bottom reads
            for i in range(top_guides, len(guide_names)):
                ground_truth[guide_names[i]] = bottom_reads // (len(guide_names) - top_guides)
                
        # Generate FASTQ file
        fastq_path = self.output_dir / f"test_reads_{num_reads}.fq"
        with open(fastq_path, 'w') as f:
            read_id = 0
            for guide_name, count in ground_truth.items():
                guide_seq = self.guides[guide_name]
                for _ in range(count):
                    # Add some context around the guide
                    full_seq = "ACCG" + guide_seq + "GTTT"  # Common CRISPR contexts
                    qual = "I" * len(full_seq)  # High quality
                    
                    f.write(f"@read_{read_id}\n")
                    f.write(f"{full_seq}\n")
                    f.write(f"+\n")
                    f.write(f"{qual}\n")
                    read_id += 1
                    
        return str(fastq_path), ground_truth
    
    def benchmark_vecmap(self, fastq_path: str, num_reads: int) -> Dict:
        """Benchmark VecMap CRISPR guide detection."""
        print("\nBenchmarking VecMap...")
        
        # Read sequences from FASTQ
        reads = []
        with open(fastq_path, 'r') as f:
            lines = f.readlines()
            for i in range(0, len(lines), 4):
                seq = lines[i+1].strip()
                read_id = lines[i].strip()[1:]
                reads.append((seq, read_id))
        
        # Initialize detector
        detector = CRISPRGuideDetector(self.guides)
        
        # Measure performance
        start_time = time.time()
        start_memory = self._get_memory_usage()
        
        # Detect guides with context
        results = detector.detect_guides_with_context(reads, 
                                                      upstream_context="ACCG",
                                                      downstream_context="GTTT")
        counts = detector.summarize_detection(results)
        
        end_time = time.time()
        end_memory = self._get_memory_usage()
        
        # Calculate metrics
        elapsed_time = end_time - start_time
        reads_per_second = num_reads / elapsed_time
        memory_used = max(0, end_memory - start_memory)
        
        return {
            'tool': 'VecMap',
            'num_reads': num_reads,
            'elapsed_time': elapsed_time,
            'reads_per_second': reads_per_second,
            'memory_mb': memory_used,
            'guide_counts': counts
        }
    
    def benchmark_crispresso2(self, fastq_path: str, num_reads: int) -> Dict:
        """Benchmark CRISPResso2 for guide counting."""
        print("\nBenchmarking CRISPResso2...")
        
        # Create amplicon sequence (representative)
        # CRISPResso2 needs an amplicon sequence to search within
        amplicon = "ACCG" + "N" * 20 + "GTTT"  # Template with guide region
        
        # Create guide file
        guide_file = self.output_dir / "guides_crispresso.txt"
        with open(guide_file, 'w') as f:
            for name, seq in self.guides.items():
                f.write(f"{name}\t{seq}\n")
        
        # Run CRISPResso2
        output_folder = self.output_dir / "crispresso_output"
        cmd = [
            "CRISPResso",
            "-r1", fastq_path,
            "-a", amplicon,
            "-g", str(guide_file),
            "-o", str(output_folder),
            "--name", "benchmark",
            "--no_rerun",
            "--suppress_report"
        ]
        
        start_time = time.time()
        start_memory = self._get_memory_usage()
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"CRISPResso2 error: {result.stderr}")
                return self._create_error_result('CRISPResso2', num_reads, "Tool execution failed")
        except FileNotFoundError:
            print("CRISPResso2 not found. Please install with: pip install CRISPResso2")
            return self._create_error_result('CRISPResso2', num_reads, "Tool not installed")
        
        end_time = time.time()
        end_memory = self._get_memory_usage()
        
        # Parse results (if needed - CRISPResso2 output is complex)
        elapsed_time = end_time - start_time
        reads_per_second = num_reads / elapsed_time
        memory_used = max(0, end_memory - start_memory)
        
        return {
            'tool': 'CRISPResso2',
            'num_reads': num_reads,
            'elapsed_time': elapsed_time,
            'reads_per_second': reads_per_second,
            'memory_mb': memory_used,
            'guide_counts': {}  # Would need to parse CRISPResso2 output files
        }
    
    def benchmark_mageck(self, fastq_path: str, num_reads: int) -> Dict:
        """Benchmark MAGeCK count module."""
        print("\nBenchmarking MAGeCK...")
        
        # Create library file for MAGeCK
        library_file = self.output_dir / "library_mageck.txt"
        with open(library_file, 'w') as f:
            f.write("sgRNA\tGene\tSequence\n")
            for name, seq in self.guides.items():
                gene = name.split('_')[0]  # Extract gene name
                f.write(f"{name}\t{gene}\t{seq}\n")
        
        # MAGeCK count command
        output_prefix = self.output_dir / "mageck_count"
        cmd = [
            "mageck", "count",
            "-l", str(library_file),
            "--sample-label", "test",
            "--fastq", fastq_path,
            "-n", str(output_prefix),
            "--norm-method", "none"  # No normalization for benchmark
        ]
        
        start_time = time.time()
        start_memory = self._get_memory_usage()
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"MAGeCK error: {result.stderr}")
                return self._create_error_result('MAGeCK', num_reads, "Tool execution failed")
        except FileNotFoundError:
            print("MAGeCK not found. Please install with: pip install mageck")
            return self._create_error_result('MAGeCK', num_reads, "Tool not installed")
        
        end_time = time.time()
        end_memory = self._get_memory_usage()
        
        # Parse MAGeCK count results
        count_file = f"{output_prefix}.count.txt"
        guide_counts = {}
        if os.path.exists(count_file):
            df = pd.read_csv(count_file, sep='\t')
            if 'test' in df.columns:
                guide_counts = dict(zip(df['sgRNA'], df['test']))
        
        elapsed_time = end_time - start_time
        reads_per_second = num_reads / elapsed_time
        memory_used = max(0, end_memory - start_memory)
        
        return {
            'tool': 'MAGeCK',
            'num_reads': num_reads,
            'elapsed_time': elapsed_time,
            'reads_per_second': reads_per_second,
            'memory_mb': memory_used,
            'guide_counts': guide_counts
        }
    
    def _get_memory_usage(self) -> float:
        """Get current memory usage in MB."""
        try:
            import psutil
            process = psutil.Process(os.getpid())
            return process.memory_info().rss / 1024 / 1024
        except ImportError:
            return 0.0
    
    def _create_error_result(self, tool: str, num_reads: int, error: str) -> Dict:
        """Create result dict for failed benchmark."""
        return {
            'tool': tool,
            'num_reads': num_reads,
            'elapsed_time': None,
            'reads_per_second': 0,
            'memory_mb': 0,
            'guide_counts': {},
            'error': error
        }
    
    def compare_accuracy(self, result: Dict, ground_truth: Dict[str, int]) -> float:
        """Compare detected counts against ground truth."""
        if 'error' in result:
            return 0.0
            
        detected = result.get('guide_counts', {})
        if not detected:
            return 0.0
            
        # Calculate correlation between detected and true counts
        guides = sorted(ground_truth.keys())
        true_counts = [ground_truth[g] for g in guides]
        detected_counts = [detected.get(g, 0) for g in guides]
        
        if sum(detected_counts) == 0:
            return 0.0
            
        correlation = np.corrcoef(true_counts, detected_counts)[0, 1]
        return correlation if not np.isnan(correlation) else 0.0
    
    def run_benchmark(self, read_counts: List[int], guide_distribution: str = "uniform"):
        """Run complete benchmark across all tools and datasets."""
        
        for num_reads in read_counts:
            print(f"\n{'='*60}")
            print(f"Testing with {num_reads:,} reads ({guide_distribution} distribution)")
            print(f"{'='*60}")
            
            # Generate test data
            fastq_path, ground_truth = self.generate_test_data(num_reads, guide_distribution)
            
            # Benchmark each tool
            for benchmark_func in [self.benchmark_vecmap, self.benchmark_crispresso2, self.benchmark_mageck]:
                result = benchmark_func(fastq_path, num_reads)
                
                # Calculate accuracy
                if 'guide_counts' in result and result['guide_counts']:
                    accuracy = self.compare_accuracy(result, ground_truth)
                    result['accuracy'] = accuracy
                else:
                    result['accuracy'] = 0.0
                
                self.results.append(result)
                
                # Print summary
                if result.get('elapsed_time'):
                    print(f"\n{result['tool']} Results:")
                    print(f"  Time: {result['elapsed_time']:.2f} seconds")
                    print(f"  Speed: {result['reads_per_second']:,.0f} reads/second")
                    print(f"  Memory: {result['memory_mb']:.1f} MB")
                    print(f"  Accuracy: {result['accuracy']:.3f}")
            
            # Clean up test file
            os.remove(fastq_path)
    
    def save_results(self):
        """Save benchmark results to CSV."""
        df = pd.DataFrame(self.results)
        output_file = self.output_dir / "crispr_tools_benchmark_results.csv"
        df.to_csv(output_file, index=False)
        print(f"\nResults saved to: {output_file}")
        
        # Create summary
        print("\n" + "="*60)
        print("BENCHMARK SUMMARY")
        print("="*60)
        
        # Average performance by tool
        for tool in df['tool'].unique():
            tool_data = df[df['tool'] == tool]
            # Filter out rows with null elapsed_time
            elapsed_times = tool_data['elapsed_time'].tolist()
            valid_indices = [i for i, t in enumerate(elapsed_times) if t is not None and not pd.isna(t)]
            if len(valid_indices) > 0:
                valid_data = tool_data.iloc[valid_indices]
                avg_speed = valid_data['reads_per_second'].mean()
                avg_memory = valid_data['memory_mb'].mean()
                avg_accuracy = valid_data['accuracy'].mean()
                
                print(f"\n{tool}:")
                print(f"  Average speed: {avg_speed:,.0f} reads/second")
                print(f"  Average memory: {avg_memory:.1f} MB")
                print(f"  Average accuracy: {avg_accuracy:.3f}")
        
        return df


def create_realistic_guide_library(num_guides: int = 1000) -> Dict[str, str]:
    """Create a realistic CRISPR guide library."""
    # Common gene targets in CRISPR screens
    gene_targets = ['TP53', 'KRAS', 'EGFR', 'MYC', 'BCL2', 'AKT1', 'PTEN', 'RB1', 
                    'CDKN2A', 'PIK3CA', 'BRAF', 'NRAS', 'JAK2', 'FLT3', 'IDH1',
                    'VHL', 'NOTCH1', 'SMAD4', 'APC', 'MLH1'] * 50  # Repeat to get enough
    
    guides = {}
    nucleotides = ['A', 'C', 'G', 'T']
    
    for i in range(num_guides):
        gene = gene_targets[i % len(gene_targets)]
        guide_num = (i // len(gene_targets)) + 1
        
        # Generate random 20bp guide
        guide_seq = ''.join(np.random.choice(nucleotides, 20))
        
        # Ensure guide starts with G (for U6 promoter compatibility)
        guide_seq = 'G' + guide_seq[1:]
        
        guide_name = f"{gene}_sg{guide_num}"
        guides[guide_name] = guide_seq
    
    return guides


def main():
    parser = argparse.ArgumentParser(description='Benchmark CRISPR guide detection tools')
    parser.add_argument('--num-guides', type=int, default=1000,
                        help='Number of guides in library (default: 1000)')
    parser.add_argument('--read-counts', type=str, default='10000,50000,100000,500000',
                        help='Comma-separated read counts to test (default: 10000,50000,100000,500000)')
    parser.add_argument('--distribution', type=str, default='power_law',
                        choices=['uniform', 'power_law'],
                        help='Guide count distribution (default: power_law)')
    parser.add_argument('--output-dir', type=str, default='benchmarks/results/crispr_comparison',
                        help='Output directory for results')
    
    args = parser.parse_args()
    
    # Parse read counts
    read_counts = [int(x) for x in args.read_counts.split(',')]
    
    print(f"CRISPR Tools Benchmark")
    print(f"=====================")
    print(f"Guide library size: {args.num_guides}")
    print(f"Read counts: {read_counts}")
    print(f"Distribution: {args.distribution}")
    
    # Create guide library
    guides = create_realistic_guide_library(args.num_guides)
    
    # Run benchmark
    benchmark = CRISPRBenchmark(guides, args.output_dir)
    benchmark.run_benchmark(read_counts, args.distribution)
    results_df = benchmark.save_results()
    
    # Calculate speedup factors
    print("\n" + "="*60)
    print("VECMAP PERFORMANCE ADVANTAGE")
    print("="*60)
    
    vecmap_data = results_df[results_df['tool'] == 'VecMap']
    
    for tool in ['CRISPResso2', 'MAGeCK']:
        tool_data = results_df[results_df['tool'] == tool]
        # Filter out rows with null elapsed_time
        elapsed_times = tool_data['elapsed_time'].tolist()
        valid_indices = [i for i, t in enumerate(elapsed_times) if t is not None and not pd.isna(t)]
        if len(valid_indices) > 0:
            valid_tool_data = tool_data.iloc[valid_indices]
            vecmap_avg = vecmap_data['reads_per_second'].mean()
            tool_avg = valid_tool_data['reads_per_second'].mean()
            
            if tool_avg > 0:
                speedup = vecmap_avg / tool_avg
                print(f"\nVecMap vs {tool}: {speedup:.1f}x faster")


if __name__ == "__main__":
    main() 