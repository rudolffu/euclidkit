"""
Core functionality for euclidkit package.

This module provides the fundamental data access and processing capabilities.
"""

from .data_access import EuclidArchive
from .spectra import SpectrumLoader, SpectrumProcessor, SpectrumCompiler
from .cutouts import CutoutGenerator, make_cutouts

__all__ = [
    "EuclidArchive",
    "SpectrumLoader", "SpectrumProcessor", "SpectrumCompiler",
    "CutoutGenerator", "make_cutouts",
]