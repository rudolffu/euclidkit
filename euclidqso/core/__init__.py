"""
Core functionality for euclidqso package.

This module provides the fundamental data access and processing capabilities.
"""

from .data_access import EuclidArchive
from .spectra import SpectrumLoader, SpectrumProcessor, SpectrumCompiler

__all__ = [
    "EuclidArchive",
    "SpectrumLoader", "SpectrumProcessor", "SpectrumCompiler",
]