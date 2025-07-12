#!/usr/bin/env python3
"""
Quick test of VecMap on transcriptomic-like data
This simulates RNA-seq reads without downloading large datasets
"""

import time
import random
import numpy as np
from vecmap import vecmap

def generate_transcriptome(num_transcripts=100, min_length=200, max_length=5000):
    """Generate a simulated transcriptome with realistic features"""
    random.seed(42)
    
    transcripts = []
    transcript_info = []
    
    # Common sequence motifs found in human transcripts
    motifs = [
        "ATGGCG",  # Start codon region
        "TATAAA",  # TATA box
        "AATAAA",  # PolyA signal
        "GTAAGT",  # Splice donor
        "CAGGT",   # Splice acceptor
    ]
    
    for i in range(num_transcripts):
        length = random.randint(min_length, max_length)
        
        # Build transcript with some structure
        transcript = []
        
        # 5' UTR
        utr5_len = random.randint(50, 200)
        transcript.extend([random.choice('ACGT') for _ in range(utr5_len)])
        
        # Insert start codon
        transcript.extend(list("ATG"))
        
        # Coding sequence with occasional motifs
        cds_len = length - utr5_len - 100  # Leave room for 3' UTR
        for j in range(cds_len):
            if random.random() < 0.01:  # 1% chance of motif
                motif = random.choice(motifs)
                transcript.extend(list(motif))
                j += len(motif)
            else:
                transcript.append(random.choice('ACGT'))
        
        # 3' UTR with polyA signal
        transcript.extend([random.choice('ACGT') for _ in range(50)])
        transcript.extend(list("AATAAA"))  # PolyA signal
        transcript.extend(['A'] * 30)  # PolyA tail
        
        transcript_seq = ''.join(transcript[:length])
        transcripts.append(transcript_seq)
        transcript_info.append({
            'id': f'ENST{i:08d}',
            'length': len(transcript_seq),
            'gc_content': (transcript_seq.count('G') + transcript_seq.count('C')) / len(transcript_seq)
        })
    
    # Concatenate with padding
    padding = 'N' * 100
    ref_sequence = padding.join(transcripts)
    
    # Create position mapping
    position_map = []
    current_pos = 0
    for i, seq in enumerate(transcripts):
        position_map.append({
            'transcript_id': transcript_info[i]['id'],
            'start': current_pos,
            'end': current_pos + len(seq),
            'length': len(seq)
        })
        current_pos += len(seq) + len(padding)
    
    return ref_sequence, transcript_info, position_map

def simulate_rnaseq_reads(ref_sequence, position_map, num_reads, read_length=100, error_rate=0.01):
    """Simulate RNA-seq reads with realistic characteristics"""
    reads = []
    
    # Filter transcripts that can accommodate reads
    valid_transcripts = [t for t in position_map if t['length'] >= read_length]
    
    # Simulate expression levels (some transcripts more abundant)
    expression_levels = np.random.lognormal(0, 2, len(valid_transcripts))
    expression_levels = expression_levels / expression_levels.sum()
    
    for _ in range(num_reads):
        # Select transcript based on expression level
        transcript_idx = np.random.choice(len(valid_transcripts), p=expression_levels)
        transcript = valid_transcripts[transcript_idx]
        
        # Position bias: reads more likely from middle of transcript
        relative_pos = np.random.beta(2, 2)  # Beta distribution peaks in middle
        max_start_pos = transcript['end'] - read_length
        read_start = int(transcript['start'] + relative_pos * (max_start_pos - transcript['start']))
        
        # Extract read
        read_seq = ref_sequence[read_start:read_start + read_length]
        
        # Quality degradation towards 3' end
        error_rate_adjusted = error_rate * (1 + relative_pos * 0.5)
        
        # Introduce sequencing errors
        read_list = list(read_seq)
        for i in range(len(read_list)):
            if random.random() < error_rate_adjusted:
                if read_list[i] in 'ACGT':
                    read_list[i] = random.choice([b for b in 'ACGT' if b != read_list[i]])
        
        reads.append((''.join(read_list), read_start))
    
    return reads

def analyze_mapping_results(mappings, position_map):
    """Analyze mapping results in context of transcripts"""
    correct = 0
    transcript_hits = {}
    
    for mapping in mappings:
        mapped_pos, mismatches, true_pos = mapping
        
        if mapped_pos == true_pos:
            correct += 1
        
        # Find which transcript this maps to
        if mapped_pos != -1:
            for transcript in position_map:
                if transcript['start'] <= mapped_pos < transcript['end']:
                    tid = transcript['transcript_id']
                    if tid not in transcript_hits:
                        transcript_hits[tid] = 0
                    transcript_hits[tid] += 1
                    break
    
    return {
        'accuracy': correct / len(mappings) * 100,
        'transcript_hits': transcript_hits,
        'num_transcripts_hit': len(transcript_hits)
    }

def run_benchmark_suite():
    """Run comprehensive benchmarks"""
    print("=== VecMap Transcriptomic Data Benchmark ===\n")
    
    # Test different scales
    test_configs = [
        {'transcripts': 50, 'reads': 1000},
        {'transcripts': 100, 'reads': 5000},
        {'transcripts': 200, 'reads': 10000},
        {'transcripts': 500, 'reads': 25000},
    ]
    
    results = []
    
    for config in test_configs:
        print(f"\nTest: {config['transcripts']} transcripts, {config['reads']} reads")
        print("-" * 50)
        
        # Generate transcriptome
        ref_seq, transcript_info, position_map = generate_transcriptome(config['transcripts'])
        print(f"Generated transcriptome: {len(ref_seq):,} bp total")
        
        # Calculate transcriptome statistics
        lengths = [t['length'] for t in position_map]
        print(f"Transcript lengths: min={min(lengths)}, max={max(lengths)}, avg={np.mean(lengths):.0f}")
        
        # Generate reads
        reads = simulate_rnaseq_reads(ref_seq, position_map, config['reads'])
        
        # Run VecMap
        start_time = time.time()
        mappings = vecmap(ref_seq, reads, read_len=100)
        end_time = time.time()
        
        # Analyze results
        total_time = end_time - start_time
        reads_per_second = len(reads) / total_time
        
        mapped_count = sum(1 for m in mappings if m[0] != -1)
        mapping_rate = mapped_count / len(reads) * 100
        
        analysis = analyze_mapping_results(mappings, position_map)
        
        # Print results
        print(f"\nResults:")
        print(f"  Runtime: {total_time:.2f} seconds")
        print(f"  Speed: {reads_per_second:.0f} reads/second")
        print(f"  Mapping rate: {mapping_rate:.1f}%")
        print(f"  Accuracy: {analysis['accuracy']:.1f}%")
        print(f"  Transcripts hit: {analysis['num_transcripts_hit']} / {len(position_map)}")
        
        # Show top expressed transcripts
        if analysis['transcript_hits']:
            top_hits = sorted(analysis['transcript_hits'].items(), 
                            key=lambda x: x[1], reverse=True)[:5]
            print(f"\n  Top mapped transcripts:")
            for tid, count in top_hits:
                print(f"    {tid}: {count} reads")
        
        results.append({
            'transcripts': config['transcripts'],
            'reads': config['reads'],
            'ref_length': len(ref_seq),
            'runtime': total_time,
            'reads_per_sec': reads_per_second,
            'mapping_rate': mapping_rate,
            'accuracy': analysis['accuracy']
        })
    
    # Summary comparison
    print("\n" + "="*60)
    print("SUMMARY - VecMap Performance on Transcriptomic Data")
    print("="*60)
    print(f"{'Transcripts':<12} {'Reads':<8} {'Ref Size':<12} {'Time (s)':<10} {'Reads/s':<10} {'Accuracy':<10}")
    print("-"*60)
    
    for r in results:
        print(f"{r['transcripts']:<12} {r['reads']:<8} {r['ref_length']:<12,} "
              f"{r['runtime']:<10.2f} {r['reads_per_sec']:<10.0f} {r['accuracy']:<10.1f}%")
    
    # Performance scaling analysis
    if len(results) > 1:
        print(f"\nPerformance Scaling:")
        base_speed = results[0]['reads_per_sec']
        for r in results:
            scaling = r['reads_per_sec'] / base_speed
            print(f"  {r['reads']} reads: {scaling:.2f}x relative throughput")

def demonstrate_splice_aware_mapping():
    """Demonstrate handling of splice junctions (future enhancement)"""
    print("\n\n=== Splice Junction Handling Demo ===")
    print("Note: Current VecMap doesn't handle splice junctions.")
    print("This demonstrates the challenge for RNA-seq mapping:\n")
    
    # Create a simple gene with two exons
    exon1 = "ATGGCGAAATTTCCCCGGGAAATTT"
    intron = "GTAAGT" + "N"*100 + "CAGGT"  # Splice sites
    exon2 = "AAACCCGGGTTTAAA"
    
    gene = exon1 + intron + exon2
    transcript = exon1 + exon2  # Spliced transcript
    
    print(f"Gene structure: Exon1({len(exon1)}bp) - Intron({len(intron)}bp) - Exon2({len(exon2)}bp)")
    print(f"Spliced transcript: {len(transcript)}bp")
    
    # Create a read spanning the junction
    junction_read = transcript[15:35]  # Spans exon1-exon2 junction
    print(f"\nJunction-spanning read: {junction_read}")
    print("This read exists in the transcript but not in the genome!")
    print("\nFuture work: Implement splice-aware mapping for accurate RNA-seq alignment")

if __name__ == '__main__':
    # Run main benchmarks
    run_benchmark_suite()
    
    # Show splice junction challenge
    demonstrate_splice_aware_mapping()
    
    print("\n\nTo test with real GEO data, run: python test_geo_data.py")
    print("(This will download ~1GB of reference data)") 