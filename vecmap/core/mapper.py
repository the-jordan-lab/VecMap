import random
import time
import numpy as np
from collections import defaultdict

def generate_reference(length):
    random.seed(0)
    unit = ''.join(random.choice('ACGT') for _ in range(100))
    times = length // 100
    return unit * times

def generate_reads(ref, num_reads, read_len, error_rate=0.01):
    reads = []
    for _ in range(num_reads):
        pos = random.randint(0, len(ref) - read_len)
        read = list(ref[pos:pos + read_len])
        for i in range(read_len):
            if random.random() < error_rate:
                read[i] = random.choice([b for b in 'ACGT' if b != read[i]])
        reads.append((''.join(read), pos))
    return reads

def build_seed_index(ref, seed_len):
    index = defaultdict(list)
    for i in range(len(ref) - seed_len + 1):
        index[ref[i:i+seed_len]].append(i)
    return index

def vecmap(ref, reads, read_len, seed_len=20, seed_offsets=[0,20,40,60,80]):
    """Vectorized short read mapping function.
    
    Args:
        ref (str): Reference sequence.
        reads (list): List of (read_seq, true_pos) tuples.
        read_len (int): Length of reads.
        seed_len (int): Seed length for indexing.
        seed_offsets (list): Offsets for multi-seed extraction.
    
    Returns:
        list: Mappings as (best_pos, min_mismatches, true_pos) tuples.
    """
    index = build_seed_index(ref, seed_len)
    ref_arr = np.array(list(ref))
    mappings = []
    for read, true_pos in reads:
        candidate_starts = set()
        for offset in seed_offsets:
            seed = read[offset:offset + seed_len]
            hits = index.get(seed, [])
            for hit in hits:
                start = hit - offset
                if start >= 0 and start + read_len <= len(ref):
                    candidate_starts.add(start)
        candidate_list = sorted(candidate_starts)
        num_candidates = len(candidate_list)
        if num_candidates == 0:
            best_pos = -1
            min_mismatches = -1
        else:
            read_arr = np.array(list(read))
            starts = np.array(candidate_list)
            substrs = ref_arr[starts[:, np.newaxis] + np.arange(read_len)]
            mismatches_arr = (substrs != read_arr).sum(axis=1)
            min_mismatches = mismatches_arr.min()
            best_idx = mismatches_arr.argmin()
            best_pos = starts[best_idx]
        mappings.append((best_pos, min_mismatches, true_pos))
    return mappings

# Example usage (benchmark)
if __name__ == '__main__':
    ref_len = 1000000
    num_reads = 100
    read_len = 100
    error_rate = 0.01

    random.seed(42)
    ref = generate_reference(ref_len)
    reads = generate_reads(ref, num_reads, read_len, error_rate)

    start = time.time()
    mappings = vecmap(ref, reads, read_len)
    end = time.time()

    print(f"Total time: {end - start} seconds")
    print("Sample mappings (best_pos, mismatches, true_pos):")
    print(mappings[:5])  # First 5 for brevity
