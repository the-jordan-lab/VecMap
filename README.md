# VecMap

VecMap is a vectorized k-mer based short read mapper for accelerating short read alignment.

## Features
- 3.4x speedup over baseline through NumPy vectorization.
- Supports substitution errors.
- Simple Python implementation.

## Requirements
- Python 3.12+
- NumPy

## Installation
Install dependencies:
```
pip install -r requirements.txt
```

## Usage
The main function is `vecmap(ref, reads, read_len)` in `vecmap.py`.

Example (included in the script):
```python
# Run the benchmark
python vecmap.py
```

This generates a 1 Mbp reference, 100 reads, and prints mapping results and runtime.

Customize by providing your own reference and reads.

## Manuscript
See [manuscript.md](manuscript.md) for the full preprint.

## License
MIT
