"""Specialized applications of VecMap for single-cell and CRISPR analysis."""

from .crispr import CRISPRGuideDetector, BarcodeGuideMatcher, detect_crispr_guides
from .barcode import BarcodeProcessor, HashtagDemultiplexer, FeatureBarcodeDetector

__all__ = [
    "CRISPRGuideDetector",
    "BarcodeGuideMatcher", 
    "detect_crispr_guides",
    "BarcodeProcessor",
    "HashtagDemultiplexer",
    "FeatureBarcodeDetector"
] 