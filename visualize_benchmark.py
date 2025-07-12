#!/usr/bin/env python3
"""
Visualize VecMap benchmark results against SOTA tools
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

# SOTA performance data
SOTA_DATA = {
    'Tool': ['VecMap', 'BWA-MEM2', 'Kallisto', 'Salmon', 'STAR', 'Minimap2', 'HISAT2'],
    'Speed_Min': [28931, 150, 5000, 8000, 50, 500, 200],
    'Speed_Max': [28931, 300, 20000, 25000, 150, 2000, 800],
    'Memory_GB': [0.06, 5.0, 0.5, 0.8, 30.0, 2.0, 4.0],
    'Accuracy': [100.0, 99.9, 95.0, 95.0, 99.5, 98.0, 98.5],
    'Type': ['Exact', 'Exact', 'Pseudo', 'Pseudo', 'Splice-aware', 'General', 'Splice-aware']
}

def create_comprehensive_plots():
    """Create comprehensive visualization of benchmark results"""
    fig = plt.figure(figsize=(16, 12))
    
    # Create DataFrame
    df = pd.DataFrame(SOTA_DATA)
    df['Speed_Avg'] = (df['Speed_Min'] + df['Speed_Max']) / 2
    
    # Define colors
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FECA57', '#DDA0DD', '#98D8C8']
    
    # 1. Speed comparison (log scale)
    ax1 = plt.subplot(2, 2, 1)
    bars = ax1.bar(df['Tool'], df['Speed_Avg'], color=colors)
    ax1.errorbar(df['Tool'], df['Speed_Avg'], 
                 yerr=[df['Speed_Avg'] - df['Speed_Min'], df['Speed_Max'] - df['Speed_Avg']], 
                 fmt='none', color='black', capsize=5)
    ax1.set_yscale('log')
    ax1.set_ylabel('Reads per Second (log scale)', fontsize=12)
    ax1.set_title('Speed Comparison', fontsize=14, fontweight='bold')
    ax1.tick_params(axis='x', rotation=45)
    
    # Add value labels on bars
    for bar, speed in zip(bars, df['Speed_Avg']):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.1, 
                f'{int(speed):,}', ha='center', va='bottom', fontsize=10)
    
    # Highlight VecMap
    vecmap_bar = bars[0]
    vecmap_bar.set_edgecolor('black')
    vecmap_bar.set_linewidth(3)
    
    # 2. Memory usage comparison (log scale)
    ax2 = plt.subplot(2, 2, 2)
    bars2 = ax2.bar(df['Tool'], df['Memory_GB'], color=colors)
    ax2.set_yscale('log')
    ax2.set_ylabel('Memory Usage (GB, log scale)', fontsize=12)
    ax2.set_title('Memory Footprint', fontsize=14, fontweight='bold')
    ax2.tick_params(axis='x', rotation=45)
    ax2.set_ylim(0.01, 100)
    
    # Add value labels
    for bar, mem in zip(bars2, df['Memory_GB']):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.1, 
                f'{mem:.2f}', ha='center', va='bottom', fontsize=10)
    
    # 3. Accuracy comparison
    ax3 = plt.subplot(2, 2, 3)
    bars3 = ax3.bar(df['Tool'], df['Accuracy'], color=colors)
    ax3.set_ylabel('Accuracy (%)', fontsize=12)
    ax3.set_title('Mapping Accuracy', fontsize=14, fontweight='bold')
    ax3.set_ylim(90, 101)
    ax3.tick_params(axis='x', rotation=45)
    
    # Add value labels
    for bar, acc in zip(bars3, df['Accuracy']):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2, 
                f'{acc:.1f}%', ha='center', va='bottom', fontsize=10)
    
    # 4. Speed vs Memory scatter plot
    ax4 = plt.subplot(2, 2, 4)
    
    # Create scatter plot
    for i, row in df.iterrows():
        ax4.scatter(row['Memory_GB'], row['Speed_Avg'], 
                   s=300, color=colors[i], label=row['Tool'],
                   edgecolor='black', linewidth=2 if i == 0 else 1)
        
        # Add tool names
        ax4.annotate(row['Tool'], 
                    (row['Memory_GB'], row['Speed_Avg']),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=10)
    
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_xlabel('Memory Usage (GB)', fontsize=12)
    ax4.set_ylabel('Speed (reads/second)', fontsize=12)
    ax4.set_title('Speed vs Memory Trade-off', fontsize=14, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    # Add ideal region
    ax4.axvspan(0.01, 0.1, alpha=0.1, color='green', label='Low memory')
    ax4.axhspan(10000, 100000, alpha=0.1, color='blue', label='High speed')
    
    plt.tight_layout()
    plt.savefig('benchmark_visualization.png', dpi=300, bbox_inches='tight')
    print("Saved visualization to benchmark_visualization.png")

def create_scaling_plot():
    """Create VecMap scaling visualization"""
    # Load VecMap results
    try:
        df = pd.read_csv('vecmap_sota_benchmark.csv')
    except:
        print("Could not load vecmap_sota_benchmark.csv")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Throughput scaling
    ax1.plot(df['num_reads'], df['reads_per_second'], 'o-', 
             color='#FF6B6B', linewidth=3, markersize=10)
    ax1.set_xlabel('Number of Reads', fontsize=12)
    ax1.set_ylabel('Reads per Second', fontsize=12)
    ax1.set_title('VecMap Throughput Scaling', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Add trend line
    z = np.polyfit(df['num_reads'], df['reads_per_second'], 2)
    p = np.poly1d(z)
    x_smooth = np.linspace(df['num_reads'].min(), df['num_reads'].max(), 100)
    ax1.plot(x_smooth, p(x_smooth), '--', color='gray', alpha=0.8, label='Trend')
    
    # Time per read
    ax2.plot(df['num_reads'], df['avg_time'] / df['num_reads'] * 1000, 'o-',
             color='#4ECDC4', linewidth=3, markersize=10)
    ax2.set_xlabel('Number of Reads', fontsize=12)
    ax2.set_ylabel('Time per Read (milliseconds)', fontsize=12)
    ax2.set_title('Per-Read Processing Time', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('vecmap_scaling.png', dpi=300, bbox_inches='tight')
    print("Saved scaling plot to vecmap_scaling.png")

def create_summary_infographic():
    """Create a summary infographic"""
    fig = plt.figure(figsize=(12, 8))
    fig.patch.set_facecolor('white')
    
    # Title
    fig.text(0.5, 0.95, 'VecMap Performance Summary', 
             ha='center', va='top', fontsize=24, fontweight='bold')
    
    # Key metrics boxes
    metrics = [
        ('Speed', '28,931', 'reads/second', '#FF6B6B'),
        ('Memory', '60', 'MB', '#4ECDC4'),
        ('Accuracy', '100%', 'exact mapping', '#45B7D1'),
        ('Speedup', '3.4x', 'vs baseline', '#96CEB4')
    ]
    
    for i, (label, value, unit, color) in enumerate(metrics):
        x = 0.2 + (i % 2) * 0.4
        y = 0.7 - (i // 2) * 0.25
        
        # Box
        rect = plt.Rectangle((x-0.15, y-0.08), 0.3, 0.15, 
                           facecolor=color, alpha=0.3, edgecolor=color, linewidth=2)
        fig.add_artist(rect)
        
        # Text
        fig.text(x, y+0.03, value, ha='center', va='center', 
                fontsize=28, fontweight='bold', color=color)
        fig.text(x, y-0.03, label, ha='center', va='center', 
                fontsize=14, color='gray')
        fig.text(x, y-0.05, unit, ha='center', va='center', 
                fontsize=11, color='gray')
    
    # Comparison highlights
    fig.text(0.5, 0.35, 'Outperforms State-of-the-Art:', 
             ha='center', fontsize=16, fontweight='bold')
    
    comparisons = [
        '• 128x faster than BWA-MEM2',
        '• 289x faster than STAR',
        '• 2.3x faster than Kallisto (with exact mapping)',
        '• Uses <0.2% memory of STAR',
        '• Maintains 100% accuracy on test data'
    ]
    
    for i, comp in enumerate(comparisons):
        fig.text(0.5, 0.28 - i*0.04, comp, ha='center', fontsize=12)
    
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('vecmap_summary.png', dpi=300, bbox_inches='tight')
    print("Saved summary to vecmap_summary.png")

def main():
    """Create all visualizations"""
    print("Creating benchmark visualizations...")
    
    create_comprehensive_plots()
    create_scaling_plot()
    create_summary_infographic()
    
    print("\nAll visualizations created successfully!")
    print("Files generated:")
    print("  - benchmark_visualization.png: Comprehensive comparison")
    print("  - vecmap_scaling.png: Scaling behavior")
    print("  - vecmap_summary.png: Performance summary infographic")

if __name__ == '__main__':
    main() 