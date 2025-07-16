#!/bin/bash
# VecMap Reproduction Script
# This script reproduces all benchmarks and figures from the VecMap paper

set -euo pipefail

echo "VecMap Benchmark Reproduction Script"
echo "===================================="
echo ""

# Record system information
echo "System Information:"
echo "-------------------"
echo "Date: $(date)"
echo "Hostname: $(hostname)"
echo "OS: $(uname -a)"
echo "Python: $(python --version)"
echo "NumPy version and BLAS info:"
python -c "import numpy; numpy.show_config(); print(f'NumPy version: {numpy.__version__}')"
echo ""

# Create output directories
mkdir -p benchmark_results/{general,crispr}
mkdir -p figures

# Function to run benchmark with replicates
run_benchmark_with_replicates() {
    local script=$1
    local output_prefix=$2
    local n_replicates=3
    
    echo "Running $script (${n_replicates} replicates)..."
    
    for i in $(seq 1 $n_replicates); do
        echo "  Replicate $i of $n_replicates"
        /usr/bin/time -v python $script > "${output_prefix}_rep${i}.log" 2>&1
    done
    
    # Aggregate results
    python -c "
import pandas as pd
import numpy as np
import glob

# Read all replicate CSV files
files = glob.glob('${output_prefix}_rep*.csv')
if files:
    dfs = [pd.read_csv(f) for f in files]
    combined = pd.concat(dfs)
    
    # Calculate mean and std for numeric columns
    numeric_cols = combined.select_dtypes(include=[np.number]).columns
    summary = combined.groupby(combined.columns.difference(numeric_cols).tolist())[numeric_cols].agg(['mean', 'std'])
    
    summary.to_csv('${output_prefix}_summary.csv')
    print(f'Summary statistics saved to ${output_prefix}_summary.csv')
"
}

# 1. General performance benchmarks
echo "1. Running general performance benchmarks..."
echo "==========================================="
run_benchmark_with_replicates \
    "benchmarks/scripts/benchmark_sota.py" \
    "benchmark_results/general/vecmap_general"

# 2. CRISPR-specific benchmarks
echo ""
echo "2. Running CRISPR-specific benchmarks..."
echo "========================================"
run_benchmark_with_replicates \
    "benchmarks/scripts/crispr_comprehensive_benchmark.py" \
    "benchmark_results/crispr/vecmap_crispr"

# 3. Tool comparison benchmarks
echo ""
echo "3. Running CRISPR tool comparisons..."
echo "====================================="
run_benchmark_with_replicates \
    "benchmarks/scripts/crispr_tools_comparison.py" \
    "benchmark_results/crispr/tools_comparison"

# 4. Generate figures
echo ""
echo "4. Generating figures from benchmark data..."
echo "==========================================="
if [ -f "benchmarks/scripts/generate_figures.py" ]; then
    python benchmarks/scripts/generate_figures.py
    echo "Figures generated in docs/figures/"
else
    echo "Warning: Figure generation script not found"
fi

# 5. Generate final report
echo ""
echo "5. Generating reproducibility report..."
echo "======================================"
cat > reproducibility_report.txt << EOF
VecMap Reproducibility Report
Generated: $(date)

System Configuration:
- OS: $(uname -s) $(uname -r)
- CPU: $(grep "model name" /proc/cpuinfo | head -1 | cut -d: -f2 || sysctl -n machdep.cpu.brand_string 2>/dev/null || echo "Unknown")
- Memory: $(free -h | grep Mem | awk '{print $2}' || echo "Unknown")
- Python: $(python --version 2>&1)
- NumPy: $(python -c "import numpy; print(numpy.__version__)")
- Git commit: $(git rev-parse HEAD)

Benchmark Results Summary:
EOF

# Add summary statistics to report
if [ -f "benchmark_results/general/vecmap_general_summary.csv" ]; then
    echo "" >> reproducibility_report.txt
    echo "General Performance (reads/second, mean ± std):" >> reproducibility_report.txt
    cat benchmark_results/general/vecmap_general_summary.csv >> reproducibility_report.txt
fi

if [ -f "benchmark_results/crispr/vecmap_crispr_summary.csv" ]; then
    echo "" >> reproducibility_report.txt
    echo "CRISPR Performance (reads/second, mean ± std):" >> reproducibility_report.txt
    cat benchmark_results/crispr/vecmap_crispr_summary.csv >> reproducibility_report.txt
fi

echo ""
echo "Reproduction complete!"
echo "====================="
echo "Results saved to:"
echo "  - benchmark_results/     (raw benchmark data)"
echo "  - figures/              (generated figures)"
echo "  - reproducibility_report.txt (summary report)"
echo ""
echo "To archive results for Zenodo:"
echo "  tar -czf vecmap_reproduction_$(date +%Y%m%d).tar.gz benchmark_results figures reproducibility_report.txt" 