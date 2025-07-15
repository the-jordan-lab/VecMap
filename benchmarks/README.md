# VecMap Benchmarks

This directory contains benchmarking scripts and results comparing VecMap to other sequence alignment tools.

## Structure

```
benchmarks/
├── scripts/              # Benchmarking scripts
│   ├── benchmark_sota.py            # Main performance comparison
│   ├── crispr_comprehensive_benchmark.py  # CRISPR-specific benchmarks
│   └── crispr_tools_comparison.py   # Head-to-head CRISPR tool comparison
├── results/              # Benchmark outputs
│   ├── ultimate_benchmark_results.csv     # General performance results
│   ├── crispr_comprehensive/              # CRISPR benchmark results
│   └── crispr_comparison/                 # Tool comparison results
└── archived_experiments/ # Historical benchmarks
```

## Running Benchmarks

### General Performance Benchmark

Compare VecMap against Minimap2 and BWA-MEM:

```bash
python benchmarks/scripts/benchmark_sota.py
```

### CRISPR-Specific Benchmarks

Compare VecMap against MAGeCK and CRISPResso2:

```bash
python benchmarks/scripts/crispr_tools_comparison.py
```

## Key Results

### General Performance
- VecMap: 42,027 reads/second
- Minimap2: 173,460 reads/second
- BWA-MEM: 60,306 reads/second

### CRISPR Screening
- VecMap: 18,948 reads/second (average)
- MAGeCK: 9,973 reads/second (1.9× slower)
- CRISPResso2: 4,986 reads/second (3.8× slower)

## Requirements

Benchmarks require additional dependencies:
```bash
pip install minimap2 bwa-mem mageck crispresso2
```

## Notes

- All benchmarks use real genomic data
- Results may vary based on hardware
- VecMap is optimized for exact matching only 