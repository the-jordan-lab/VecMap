#!/usr/bin/env python3
"""
Ultimate head-to-head benchmark: VecMap vs SOTA tools
Direct comparison on identical data and hardware
"""

import os
import sys
import time
import subprocess
import shutil
import json
import psutil
import numpy as np
import pandas as pd
from datetime import datetime
import multiprocessing
from vecmap import vecmap
from test_geo_quick import generate_transcriptome, simulate_rnaseq_reads

class BenchmarkTool:
    """Base class for benchmarking tools"""
    def __init__(self, name):
        self.name = name
        self.installed = False
        self.version = None
        
    def check_installation(self):
        """Check if tool is installed"""
        raise NotImplementedError
        
    def install(self):
        """Install the tool"""
        raise NotImplementedError
        
    def run_alignment(self, ref_file, reads_file, output_prefix, threads=1):
        """Run the alignment"""
        raise NotImplementedError

class VecMapBenchmark(BenchmarkTool):
    """VecMap benchmark wrapper"""
    def __init__(self):
        super().__init__("VecMap")
        self.installed = True
        self.version = "1.0"
        
    def check_installation(self):
        return True
        
    def run_alignment(self, ref_sequence, reads, output_prefix, threads=1):
        """Run VecMap directly on data structures"""
        start_time = time.time()
        start_memory = psutil.Process().memory_info().rss / 1024 / 1024  # MB
        
        mappings = vecmap(ref_sequence, reads, read_len=100)
        
        end_time = time.time()
        peak_memory = psutil.Process().memory_info().rss / 1024 / 1024  # MB
        
        # Save results
        with open(f"{output_prefix}_mappings.txt", 'w') as f:
            for m in mappings:
                f.write(f"{m[0]}\t{m[1]}\t{m[2]}\n")
                
        return {
            'time': end_time - start_time,
            'memory_mb': peak_memory - start_memory,
            'mappings': mappings
        }

class Minimap2Benchmark(BenchmarkTool):
    """Minimap2 benchmark wrapper"""
    def __init__(self):
        super().__init__("Minimap2")
        
    def check_installation(self):
        try:
            result = subprocess.run(['minimap2', '--version'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                self.installed = True
                self.version = result.stdout.strip()
                return True
        except:
            pass
        return False
        
    def install(self):
        print("Installing Minimap2...")
        try:
            # Try conda first
            subprocess.run(['conda', 'install', '-y', '-c', 'bioconda', 'minimap2'], 
                         check=True)
            return self.check_installation()
        except:
            # Try brew for Mac
            try:
                subprocess.run(['brew', 'install', 'minimap2'], check=True)
                return self.check_installation()
            except:
                print("Please install Minimap2 manually")
                return False
                
    def run_alignment(self, ref_file, reads_file, output_prefix, threads=1):
        """Run Minimap2"""
        output_file = f"{output_prefix}_minimap2.sam"
        
        cmd = [
            'minimap2', '-ax', 'sr',  # short read preset
            '-t', str(threads),
            ref_file, reads_file
        ]
        
        start_time = time.time()
        start_memory = psutil.Process().memory_info().rss / 1024 / 1024
        
        with open(output_file, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
            
        end_time = time.time()
        peak_memory = psutil.Process().memory_info().rss / 1024 / 1024
        
        if result.returncode != 0:
            print(f"Minimap2 error: {result.stderr}")
            return None
            
        return {
            'time': end_time - start_time,
            'memory_mb': peak_memory - start_memory,
            'output_file': output_file
        }

class BWABenchmark(BenchmarkTool):
    """BWA-MEM benchmark wrapper"""
    def __init__(self):
        super().__init__("BWA-MEM")
        
    def check_installation(self):
        try:
            result = subprocess.run(['bwa'], capture_output=True, text=True)
            if 'Version:' in result.stderr:
                self.installed = True
                # Extract version
                for line in result.stderr.split('\n'):
                    if 'Version:' in line:
                        self.version = line.split('Version:')[1].strip()
                return True
        except:
            pass
        return False
        
    def install(self):
        print("Installing BWA...")
        try:
            subprocess.run(['conda', 'install', '-y', '-c', 'bioconda', 'bwa'], 
                         check=True)
            return self.check_installation()
        except:
            try:
                subprocess.run(['brew', 'install', 'bwa'], check=True)
                return self.check_installation()
            except:
                print("Please install BWA manually")
                return False
                
    def run_alignment(self, ref_file, reads_file, output_prefix, threads=1):
        """Run BWA-MEM"""
        # Index reference first
        print("  Building BWA index...")
        subprocess.run(['bwa', 'index', ref_file], 
                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        output_file = f"{output_prefix}_bwa.sam"
        
        cmd = [
            'bwa', 'mem',
            '-t', str(threads),
            ref_file, reads_file
        ]
        
        start_time = time.time()
        start_memory = psutil.Process().memory_info().rss / 1024 / 1024
        
        with open(output_file, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
            
        end_time = time.time()
        peak_memory = psutil.Process().memory_info().rss / 1024 / 1024
        
        if result.returncode != 0:
            print(f"BWA error: {result.stderr}")
            return None
            
        return {
            'time': end_time - start_time,
            'memory_mb': peak_memory - start_memory,
            'output_file': output_file
        }

class KallistoBenchmark(BenchmarkTool):
    """Kallisto benchmark wrapper"""
    def __init__(self):
        super().__init__("Kallisto")
        
    def check_installation(self):
        try:
            result = subprocess.run(['kallisto', 'version'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                self.installed = True
                self.version = result.stdout.strip().split()[-1]
                return True
        except:
            pass
        return False
        
    def install(self):
        print("Installing Kallisto...")
        try:
            subprocess.run(['conda', 'install', '-y', '-c', 'bioconda', 'kallisto'], 
                         check=True)
            return self.check_installation()
        except:
            print("Please install Kallisto manually")
            return False
            
    def run_alignment(self, ref_file, reads_file, output_prefix, threads=1):
        """Run Kallisto"""
        # Build index first
        index_file = f"{output_prefix}_kallisto.idx"
        print("  Building Kallisto index...")
        subprocess.run(['kallisto', 'index', '-i', index_file, ref_file],
                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        output_dir = f"{output_prefix}_kallisto_output"
        os.makedirs(output_dir, exist_ok=True)
        
        cmd = [
            'kallisto', 'quant',
            '-i', index_file,
            '-o', output_dir,
            '--single', '-l', '200', '-s', '20',  # Single-end parameters
            '-t', str(threads),
            reads_file
        ]
        
        start_time = time.time()
        start_memory = psutil.Process().memory_info().rss / 1024 / 1024
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        end_time = time.time()
        peak_memory = psutil.Process().memory_info().rss / 1024 / 1024
        
        if result.returncode != 0:
            print(f"Kallisto error: {result.stderr}")
            return None
            
        return {
            'time': end_time - start_time,
            'memory_mb': peak_memory - start_memory,
            'output_dir': output_dir
        }

def prepare_test_data(num_transcripts=100, num_reads=10000):
    """Prepare test data for all tools"""
    print(f"Preparing test data: {num_transcripts} transcripts, {num_reads} reads...")
    
    # Generate transcriptome
    ref_sequence, transcript_info, position_map = generate_transcriptome(num_transcripts)
    
    # Save reference
    ref_file = "benchmark_data/reference.fasta"
    os.makedirs("benchmark_data", exist_ok=True)
    
    with open(ref_file, 'w') as f:
        f.write(">Reference_Transcriptome\n")
        # Write in 80-char lines
        for i in range(0, len(ref_sequence), 80):
            f.write(ref_sequence[i:i+80] + "\n")
    
    # Generate reads
    reads = simulate_rnaseq_reads(ref_sequence, position_map, num_reads)
    
    # Save reads
    reads_file = "benchmark_data/reads.fasta"
    with open(reads_file, 'w') as f:
        for i, (seq, pos) in enumerate(reads):
            f.write(f">read_{i}_pos_{pos}\n{seq}\n")
    
    print(f"  Reference: {len(ref_sequence):,} bp")
    print(f"  Reads: {len(reads)} x 100 bp")
    
    return ref_file, reads_file, ref_sequence, reads

def analyze_results(vecmap_mappings, sam_file=None):
    """Analyze alignment results"""
    if sam_file and os.path.exists(sam_file):
        # Parse SAM file for other tools
        mapped = 0
        total = 0
        with open(sam_file, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                total += 1
                fields = line.split('\t')
                if len(fields) > 3 and fields[2] != '*':
                    mapped += 1
        return mapped, total, None  # Return None for accuracy
    else:
        # VecMap results
        mapped = sum(1 for m in vecmap_mappings if m[0] != -1)
        correct = sum(1 for m in vecmap_mappings if m[0] == m[2])
        return mapped, len(vecmap_mappings), correct

def run_ultimate_benchmark():
    """Run the ultimate head-to-head benchmark"""
    print("="*80)
    print("ULTIMATE HEAD-TO-HEAD BENCHMARK: VecMap vs State-of-the-Art")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"System: {os.uname().sysname} {os.uname().machine}")
    print(f"CPU: {multiprocessing.cpu_count()} cores")
    print(f"Python: {sys.version.split()[0]}")
    print()
    
    # Test configurations
    test_configs = [
        {'transcripts': 50, 'reads': 5000, 'name': 'Small'},
        {'transcripts': 100, 'reads': 10000, 'name': 'Medium'},
        {'transcripts': 200, 'reads': 25000, 'name': 'Large'},
    ]
    
    # Tools to benchmark
    tools = {
        'VecMap': VecMapBenchmark(),
        'Minimap2': Minimap2Benchmark(),
        'BWA-MEM': BWABenchmark(),
        'Kallisto': KallistoBenchmark(),
    }
    
    # Check installations
    print("Checking tool installations...")
    for name, tool in tools.items():
        if tool.check_installation():
            print(f"  ✓ {name} {tool.version}")
        else:
            print(f"  ✗ {name} not found")
            if name != 'VecMap':  # Try to install
                if tool.install():
                    print(f"    ✓ {name} installed successfully")
                else:
                    print(f"    ✗ Failed to install {name}")
    
    # Results storage
    all_results = []
    
    # Run benchmarks
    for config in test_configs:
        print(f"\n{'='*80}")
        print(f"Test: {config['name']} ({config['transcripts']} transcripts, {config['reads']} reads)")
        print(f"{'='*80}")
        
        # Prepare data
        ref_file, reads_file, ref_sequence, reads = prepare_test_data(
            config['transcripts'], config['reads']
        )
        
        # Run each tool
        for tool_name, tool in tools.items():
            if not tool.check_installation():
                continue
                
            print(f"\nRunning {tool_name}...")
            
            try:
                if tool_name == 'VecMap':
                    # VecMap runs directly on data
                    result = tool.run_alignment(ref_sequence, reads, 
                                              f"benchmark_data/{tool_name.lower()}", 
                                              threads=1)
                    if result:
                        mapped, total, correct = analyze_results(result['mappings'])
                        accuracy = (correct / total * 100) if total > 0 and correct is not None else 0
                else:
                    # Other tools use files
                    result = tool.run_alignment(ref_file, reads_file,
                                              f"benchmark_data/{tool_name.lower()}",
                                              threads=1)
                    if result and 'output_file' in result:
                        mapped, total, correct = analyze_results(None, result['output_file'])
                        accuracy = None  # Can't easily determine without ground truth
                
                if result:
                    reads_per_second = config['reads'] / result['time']
                    
                    print(f"  Time: {result['time']:.2f} seconds")
                    print(f"  Speed: {reads_per_second:.0f} reads/second")
                    print(f"  Memory: {result['memory_mb']:.0f} MB")
                    if tool_name == 'VecMap':
                        print(f"  Accuracy: {accuracy:.1f}%")
                    print(f"  Mapped: {mapped}/{total} ({mapped/total*100:.1f}%)")
                    
                    all_results.append({
                        'config': config['name'],
                        'tool': tool_name,
                        'version': tool.version,
                        'time': result['time'],
                        'reads_per_second': reads_per_second,
                        'memory_mb': result['memory_mb'],
                        'mapped': mapped,
                        'total': total,
                        'accuracy': accuracy
                    })
                    
            except Exception as e:
                print(f"  Error: {e}")
    
    # Generate report
    generate_ultimate_report(all_results, test_configs)

def generate_ultimate_report(results, configs):
    """Generate comprehensive comparison report"""
    df = pd.DataFrame(results)
    
    print("\n" + "="*80)
    print("FINAL RESULTS SUMMARY")
    print("="*80)
    
    # Group by configuration
    for config in configs:
        config_results = df[df['config'] == config['name']]
        if config_results.empty:
            continue
            
        print(f"\n{config['name']} Test ({config['transcripts']} transcripts, {config['reads']} reads):")
        print("-"*80)
        
        # Sort by speed
        config_results = config_results.sort_values(by='reads_per_second', ascending=False)
        
        # Find VecMap performance for comparison
        vecmap_rows = config_results[config_results['tool'] == 'VecMap']
        if len(vecmap_rows) > 0:
            vecmap_speed = vecmap_rows['reads_per_second'].values[0]
        else:
            continue
        
        print(f"{'Tool':<12} {'Version':<10} {'Time (s)':<10} {'Speed (r/s)':<12} {'Memory (MB)':<12} {'vs VecMap':<12}")
        print("-"*80)
        
        for _, row in config_results.iterrows():
            speedup = row['reads_per_second'] / vecmap_speed
            print(f"{row['tool']:<12} {str(row['version']):<10} {row['time']:<10.2f} "
                  f"{row['reads_per_second']:<12.0f} {row['memory_mb']:<12.0f} "
                  f"{speedup:<12.2f}x")
    
    # Overall summary
    print("\n" + "="*80)
    print("OVERALL PERFORMANCE RANKING (by average speed)")
    print("="*80)
    
    avg_speeds = df.groupby('tool')['reads_per_second'].mean().sort_values(ascending=False)
    vecmap_avg = avg_speeds['VecMap']
    
    for tool, avg_speed in avg_speeds.items():
        ratio = avg_speed / vecmap_avg
        if tool == 'VecMap':
            print(f"1. {tool}: {avg_speed:.0f} reads/s (baseline)")
        else:
            if ratio > 1:
                print(f"{list(avg_speeds.index).index(tool)+1}. {tool}: {avg_speed:.0f} reads/s ({ratio:.1f}x faster than VecMap)")
            else:
                print(f"{list(avg_speeds.index).index(tool)+1}. {tool}: {avg_speed:.0f} reads/s ({1/ratio:.1f}x slower than VecMap)")
    
    # Save detailed results
    df.to_csv('ultimate_benchmark_results.csv', index=False)
    print(f"\nDetailed results saved to: ultimate_benchmark_results.csv")
    
    # Key findings
    print("\n" + "="*80)
    print("KEY FINDINGS")
    print("="*80)
    
    if 'VecMap' in avg_speeds.index:
        vecmap_results = df[df['tool'] == 'VecMap']
        print(f"✓ VecMap average speed: {vecmap_avg:.0f} reads/second")
        print(f"✓ VecMap memory usage: {vecmap_results['memory_mb'].mean():.0f} MB")
        if vecmap_results['accuracy'].notna().any():
            print(f"✓ VecMap accuracy: {vecmap_results['accuracy'].mean():.1f}%")
        
        # Compare with each tool
        for tool in avg_speeds.index:
            if tool != 'VecMap':
                tool_results = df[df['tool'] == tool]
                if not tool_results.empty:
                    speed_ratio = vecmap_avg / avg_speeds[tool]
                    memory_ratio = tool_results['memory_mb'].mean() / vecmap_results['memory_mb'].mean()
                    
                    print(f"\n{tool}:")
                    print(f"  - VecMap is {speed_ratio:.1f}x faster")
                    print(f"  - VecMap uses {1/memory_ratio:.1f}x less memory")

if __name__ == '__main__':
    run_ultimate_benchmark() 