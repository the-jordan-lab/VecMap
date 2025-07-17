#!/usr/bin/env python3
"""
Validate that VecMap benchmarks can be reproduced

This script checks that all dependencies are present and runs a quick
benchmark to verify performance is in the expected range.
"""

import sys
import os
import time

def check_dependencies():
    """Check that all required modules are available"""
    print("Checking dependencies...")
    errors = []
    
    # Check Python packages
    try:
        import numpy
        print("✓ numpy installed")
    except ImportError:
        errors.append("numpy not installed - run: pip install numpy")
    
    try:
        import pandas
        print("✓ pandas installed")
    except ImportError:
        errors.append("pandas not installed - run: pip install pandas")
    
    try:
        import vecmap
        print("✓ vecmap installed")
    except ImportError:
        errors.append("vecmap not installed - run: pip install -e . from project root")
    
    try:
        import test_geo_quick
        print("✓ test_geo_quick module found")
    except ImportError:
        errors.append("test_geo_quick.py missing from benchmarks/scripts/")
    
    # Check directories
    dirs_needed = ['benchmark_data', 'benchmark_results', 'benchmark_indices', 'benchmark_outputs']
    for d in dirs_needed:
        if not os.path.exists(d):
            os.makedirs(d, exist_ok=True)
            print(f"✓ Created {d}/ directory")
        else:
            print(f"✓ {d}/ directory exists")
    
    return errors


def run_quick_benchmark():
    """Run a quick benchmark to validate performance"""
    print("\nRunning quick benchmark...")
    
    from test_geo_quick import generate_transcriptome, simulate_rnaseq_reads
    from vecmap import vecmap
    
    # Small test
    ref_sequence, _, position_map = generate_transcriptome(100)
    reads = simulate_rnaseq_reads(ref_sequence, position_map, 5000)
    
    start = time.time()
    results = vecmap(ref_sequence, reads, 100)
    elapsed = time.time() - start
    
    speed = len(reads) / elapsed
    print(f"\nResults:")
    print(f"  Time: {elapsed:.3f} seconds")
    print(f"  Speed: {speed:,.0f} reads/second")
    print(f"  Expected range: 35,000-50,000 reads/second")
    
    if 35000 <= speed <= 50000:
        print("\n✓ Performance is in expected range!")
        return True
    else:
        print("\n⚠️  Performance is outside expected range")
        print("  This may be due to different hardware or system load")
        return False


def main():
    print("VecMap Benchmark Validation")
    print("="*50)
    
    # Check dependencies
    errors = check_dependencies()
    
    if errors:
        print("\n❌ Missing dependencies:")
        for error in errors:
            print(f"  - {error}")
        print("\nPlease install missing dependencies and try again.")
        sys.exit(1)
    
    print("\n✓ All dependencies satisfied!")
    
    # Run benchmark
    if run_quick_benchmark():
        print("\n" + "="*50)
        print("✓ Benchmarks are ready to run!")
        print("\nTo run full benchmarks:")
        print("  python benchmark_sota.py")
        print("\nTo run CRISPR benchmarks:")
        print("  python crispr_comprehensive_benchmark.py")
    else:
        print("\nBenchmarks can still be run, but results may vary.")


if __name__ == "__main__":
    main() 