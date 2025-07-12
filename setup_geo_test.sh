#!/bin/bash
# Setup script for VecMap GEO testing

echo "Setting up VecMap GEO testing environment..."

# Create virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
source venv/bin/activate

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip

# Install requirements
echo "Installing Python dependencies..."
pip install -r requirements.txt

# Check if conda is available for SRA tools
if command -v conda &> /dev/null; then
    echo "Installing SRA tools via conda..."
    conda install -y -c bioconda sra-tools
else
    echo "Warning: Conda not found. You may need to install SRA tools manually for full functionality."
    echo "Visit: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit"
fi

echo ""
echo "Setup complete! To run the tests:"
echo ""
echo "1. Quick test (no downloads required):"
echo "   python test_geo_quick.py"
echo ""
echo "2. Full GEO test (downloads ~1GB reference data):"
echo "   python test_geo_data.py"
echo ""
echo "Make sure to activate the virtual environment first:"
echo "   source venv/bin/activate" 