"""
euclidqso: Euclid QSO Analysis Package

A comprehensive toolkit for analyzing QSO candidates from Euclid data,
including spectroscopic and photometric analysis, template generation,
and integration with external surveys like DESI.
"""

from euclidqso.version import __version__
import importlib

__all__ = ["__version__", "data_access", "spectra", "io"]

_lazy_modules = {
    "data_access": "euclidqso.core.data_access",
    "spectra": "euclidqso.core.spectra",
    "io": "euclidqso.utils.io",
}

# Package metadata
__author__ = "Yuming Fu"
__email__ = "fuympku@outlook.com" 
__license__ = "GPLv3"
__url__ = "https://github.com/rudolffu/euclidqso"

# Lazy module loading to avoid heavy network-bound imports at CLI startup.
def __getattr__(name):
    if name in _lazy_modules:
        module = importlib.import_module(_lazy_modules[name])
        globals()[name] = module
        return module
    raise AttributeError(f"module 'euclidqso' has no attribute '{name}'")

# Set up logging
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())
