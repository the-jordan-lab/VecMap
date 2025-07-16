"""
VecMap: Vectorized K-mer Alignment for Python
Setup script for distribution
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="vecmap",
    version="1.0.0",
    author="James M. Jordan",
    author_email="jjordan@bio.fsu.edu",
    description="NumPy-vectorized sequence alignment for transcriptomes and exact matching",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/the-jordan-lab/VecMap",
    project_urls={
        "Bug Reports": "https://github.com/the-jordan-lab/VecMap/issues",
        "Source": "https://github.com/the-jordan-lab/VecMap",
        "Paper": "https://www.biorxiv.org/content/10.1101/2025.XX.XX.XXXXXX",
    },
    packages=find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.20.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
        ],
        "viz": [
            "matplotlib>=3.3.0",
            "seaborn>=0.11.0",
            "pandas>=1.3.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "vecmap=vecmap.cli.main:main",
        ],
    },
    keywords=[
        "bioinformatics",
        "sequence alignment", 
        "RNA-seq",
        "transcriptomics",
        "CRISPR",
        "single-cell",
        "barcode demultiplexing",
        "numpy",
        "vectorization"
    ],
) 