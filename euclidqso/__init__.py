"""
euclidqso: Euclid QSO Analysis Package

A comprehensive toolkit for analyzing QSO candidates from Euclid data,
including spectroscopic and photometric analysis, template generation,
and integration with external surveys like DESI.
"""

from euclidqso.version import __version__

# Core functionality
from euclidqso.core import data_access, spectra

# Analysis modules (will be available as modules are implemented)
# from euclidqso.analysis import (
#     photometry, spectroscopy, redshift, templates, composite, qso_fitting
# )

# Visualization (will be available as modules are implemented)
# from euclidqso.visualization import plots, specplots, photplots

# External integrations (will be available as modules are implemented)
# from euclidqso.external import desi, wise, gaia, des

# Utilities
from euclidqso.utils import io

# Configuration
# from euclidqso import config

__all__ = [
    "__version__",
    # Core
    "data_access", "spectra",
    # Utils
    "io",
]

# Package metadata
__author__ = "Yuming Fu"
__email__ = "fuympku@outlook.com" 
__license__ = "GPLv3"
__url__ = "https://github.com/rudolffu/euclidqso"

# Set up logging
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())