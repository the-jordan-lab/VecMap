# VecMap

[![PyPI version](https://badge.fury.io/py/vecmap.svg)](https://badge.fury.io/py/vecmap)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

Ultrafast exact sequence matching using NumPy vectorization. Designed for CRISPR screens and barcode mapping where exact matching is biologically required.

**Paper**: [bioRxiv preprint](https://doi.org/10.1101/2025.XX.XX.XXXXXX)

## Installation

```bash
pip install vecmap
```

## Quick Start

```bash
# Command line
vecmap -r reference.fa -q reads.fq -o alignments.txt

# Python API
from vecmap import vecmap
alignments = vecmap(reference_seq, [(read_seq, read_id), ...])
```

## Performance

Benchmarks on human transcriptome (42K reads/second in pure Python):

| Tool | Reads/sec | Language | Notes |
|------|-----------|----------|-------|
| VecMap | 42,027 | Python | Exact matching only |
| Minimap2 | 173,460 | C | General purpose aligner |
| BWA-MEM | 60,306 | C | General purpose aligner |

For CRISPR screening specifically:
- VecMap: 18,948 reads/sec (average)
- 1.9× faster than MAGeCK
- 3.8× faster than CRISPResso2

## Applications

### CRISPR Guide Detection

```python
from vecmap.applications import CRISPRGuideDetector

guides = {
    "KRAS_sg1": "ACGTACGTACGTACGTACGT",
    "TP53_sg1": "GGCCGGCCGGCCGGCCGGCC"
}

detector = CRISPRGuideDetector(guides)
results = detector.detect_guides(reads)
counts = detector.summarize_detection(results)
```

### Barcode Demultiplexing

```python
from vecmap.applications import BarcodeProcessor

processor = BarcodeProcessor(
    barcode_whitelist=whitelist,
    barcode_length=16,
    umi_length=10
)

corrected = processor.correct_barcodes(processor.extract_barcodes(reads))
```

## Documentation

- [Examples](examples/) - Usage examples and tutorials
- [Benchmarks](benchmarks/) - Performance benchmarks and comparisons
- [API Reference](https://vecmap.readthedocs.io) - Full API documentation

## Citation

```bibtex
@article{jordan2025vecmap,
  title={VecMap: Ultrafast Exact Sequence Matching for CRISPR Screens},
  author={Jordan, James M},
  journal={bioRxiv},
  year={2025},
  doi={10.1101/2025.XX.XX.XXXXXX}
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.
