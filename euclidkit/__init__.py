"""
euclidkit: Euclid QSO Analysis Package

A comprehensive toolkit for analyzing QSO candidates from Euclid data,
including spectroscopic and photometric analysis, template generation,
and integration with external surveys like DESI.
"""

from euclidkit.version import __version__
import importlib

__all__ = ["__version__", "data_access", "spectra", "io"]

_lazy_modules = {
    "data_access": "euclidkit.core.data_access",
    "spectra": "euclidkit.core.spectra",
    "io": "euclidkit.utils.io",
}

# Package metadata
__author__ = "Yuming Fu"
__email__ = "fuympku@outlook.com" 
__license__ = "GPLv3"
__url__ = "https://github.com/rudolffu/euclidkit"

# Lazy module loading to avoid heavy network-bound imports at CLI startup.
def __getattr__(name):
    if name in _lazy_modules:
        module = importlib.import_module(_lazy_modules[name])
        globals()[name] = module
        return module
    raise AttributeError(f"module 'euclidkit' has no attribute '{name}'")

# Set up logging
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())
