# VecMap Performance Guide

### 1. CRISPR Guide Detection Performance

**Key Factor: Guide Library Size**

| Library Size | Performance | Use Cases |
|-------------|-------------|-----------|
| <100 guides | ~100,000 reads/sec | Focused screens, validation |
| ~500 guides | ~40,000 reads/sec | Pathway screens |
| ~1,000 guides | ~20,000 reads/sec | Mid-scale screens |
| ~5,000 guides | ~5,000 reads/sec | Genome-wide screens |

**Why?** Each guide requires a k-mer index entry. More guides = more index lookups.

### 2. General Sequence Alignment Performance

**Key Factor: Reference Sequence Diversity**

| Reference Type | Performance | Explanation |
|----------------|-------------|-------------|
| Diverse sequences (transcriptome) | 40,000-50,000 reads/sec | Few k-mer collisions |
| Repetitive sequences | 500-2,000 reads/sec | Many k-mer collisions |
| Synthetic repeats (generate_reference) | ~500 reads/sec | Worst case scenario |

### 3. Optimal Use Cases

VecMap excels when:
- ✅ **Exact matching is required** (no mismatches/indels)
- ✅ **Reference has diverse k-mers** (transcriptomes, guide libraries)
- ✅ **Guide libraries are focused** (<1,000 guides)
- ✅ **Python integration needed** (single-cell pipelines)

VecMap struggles when:
- ❌ **Mismatches/indels needed** (use Minimap2)
- ❌ **Highly repetitive references** (use BWA-MEM)
- ❌ **Very large guide libraries** (>5,000 guides)

### 4. Performance Tips

1. **For CRISPR screens**: Keep guide libraries focused when possible
2. **For RNA-seq**: VecMap shines on transcriptomes (diverse sequences)
3. **For testing**: Don't use `generate_reference()` for benchmarks (creates worst-case repeats)

### 5. Hardware Considerations

- **CPU**: Single-threaded, benefits from high clock speed
- **Memory**: Efficient (~22MB for human transcriptome)
- **NumPy**: Uses Apple Accelerate on M-series Macs

## Benchmark Data

All benchmarks on Apple M3 Max, Python 3.12, NumPy 2.3.1:

### CRISPR Guide Detection (reads/second)
```
Guide Library Size vs Performance:
    10 guides: 121,305
   100 guides:  82,658  
   500 guides:  39,452
 1,000 guides:  23,526
 5,000 guides:   5,635
```

### Reference Type Impact (25,000 reads)
```
Transcriptome (diverse): 47,503 reads/sec
Continuous (repetitive):    468 reads/sec
```
