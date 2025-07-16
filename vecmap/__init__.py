"""
VecMap: Vectorized K-mer Alignment for Python
==============================================

A NumPy-vectorized approach to ultrafast exact sequence alignment.

Author: James M. Jordan
License: MIT
"""

__version__ = "1.0.0"
__author__ = "James M. Jordan"
__email__ = "jjordan@bio.fsu.edu"

from .core.mapper import vecmap, generate_reference, generate_reads
from .core.mapper_numba import vecmap_numba

__all__ = ["vecmap", "vecmap_numba", "generate_reference", "generate_reads"]
