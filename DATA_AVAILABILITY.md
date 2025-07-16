# Data Availability

All data and code necessary to reproduce the results in the VecMap manuscript are available as follows:

## Code
- **VecMap software**: This repository (https://github.com/the-jordan-lab/VecMap)
- **Release version**: v1.0.0 
- **Git commit**: Use tag `v1.0.0` for exact reproduction

## Benchmark Data
- **Simulated reads**: Generated using scripts in `benchmarks/scripts/`
- **Reference sequences**: Downloaded from Ensembl (human transcriptome GRCh38)
- **CRISPR libraries**: Synthetic guides generated in benchmarking scripts

## Reproduction
To fully reproduce all benchmarks and figures:

```bash
git clone https://github.com/the-jordan-lab/VecMap.git
cd VecMap
git checkout v1.0.0
pip install -e .
./reproduce.sh
```

## Hardware/Software Requirements
Benchmarks were performed on:
- **CPU**: Apple M1 Max
- **RAM**: 32GB
- **OS**: macOS 14.5
- **Python**: 3.11.5
- **NumPy**: 2.0.0 (using Apple Accelerate BLAS)

## Archived Data
Upon publication, raw benchmark outputs will be deposited to Zenodo with DOI.

## Random Seeds
All benchmarking scripts use fixed random seeds for reproducibility:
- General benchmarks: `random_seed=42`
- CRISPR benchmarks: `random_seed=12345` 