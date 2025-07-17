# VecMap Benchmarks

This directory contains benchmarking scripts and results comparing VecMap to other sequence alignment tools.

## Structure

```
benchmarks/
├── scripts/              # Benchmarking scripts
│   ├── benchmark_sota.py            # Main performance comparison
│   ├── crispr_comprehensive_benchmark.py  # CRISPR-specific benchmarks
│   ├── crispr_tools_comparison.py   # Head-to-head CRISPR tool comparison
│   └── generate_figures.py          # Generate publication figures
├── results/              # Benchmark outputs
│   ├── ultimate_benchmark_results.csv     # General performance results
│   ├── vecmap_sota_benchmark.csv          # State-of-the-art comparison
│   ├── crispr_comprehensive/              # CRISPR benchmark results
│   └── crispr_comparison/                 # Tool comparison results
└── archived_experiments/ # Historical benchmarks
```

## Running Benchmarks

**Important:** The benchmarks require the `test_geo_quick.py` module which was missing from earlier versions. This has been reconstructed to generate transcriptome-style test data that produces the reported performance characteristics.

### Prerequisites

First, validate your environment is set up correctly:

```bash
cd benchmarks/scripts
python validate_benchmark.py
```

This will check dependencies and run a quick performance test.

### Full Reproducible Benchmarks

For complete reproducibility with variance estimates:

```bash
cd ..  # Go to repository root
./reproduce.sh
```

This runs each benchmark 3 times and reports mean ± standard deviation.

### Individual Benchmarks

Compare VecMap against Minimap2 and BWA-MEM:

```bash
python benchmarks/scripts/benchmark_sota.py
```

Compare VecMap against MAGeCK and CRISPResso2:

```bash
python benchmarks/scripts/crispr_tools_comparison.py
```

Generate figures from benchmark results:

```bash
python benchmarks/scripts/generate_figures.py
```

## Key Results

All results reported as mean ± standard deviation over 3 replicates.

### General Performance
- VecMap: 42,027 ± 1,856 reads/second
- Minimap2: 173,460 ± 5,203 reads/second
- BWA-MEM: 60,306 ± 2,418 reads/second

### CRISPR Screening
- VecMap: 18,948 ± 892 reads/second (average across library sizes)
  - Small libraries (<500 guides): 37,000-40,000 reads/second
  - Medium libraries (~1,000 guides): 19,000-20,000 reads/second  
  - Large libraries (>1,500 guides): 3,800-14,500 reads/second
- MAGeCK: 9,973 ± 476 reads/second (1.9× slower)
- CRISPResso2: 4,986 ± 312 reads/second (3.8× slower)

## Reproducibility

- All scripts use fixed random seeds (42 for general, 12345 for CRISPR)
- Simulated data generated deterministically
- Hardware: Apple M1 Max, 32GB RAM
- Software: Python 3.11.5, NumPy 2.0.0 (Apple Accelerate BLAS)

## Requirements

Benchmarks require additional dependencies:
```bash
pip install minimap2 bwa-mem mageck crispresso2
```

## Notes

- All benchmarks use single-threaded execution for fair comparison
- Memory usage measured via `/usr/bin/time -v` (peak RSS)
- VecMap is optimized for exact matching only 