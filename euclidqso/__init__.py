"""
euclidqso: Euclid QSO Analysis Package

A comprehensive toolkit for analyzing QSO candidates from Euclid data,
including spectroscopic and photometric analysis, template generation,
and integration with external surveys like DESI.
"""

from euclidqso.version import __version__

# Core functionality (will be available as modules are implemented)
# from euclidqso.core import data_access, cutouts, spectra, catalogs

# Analysis modules (will be available as modules are implemented)
# from euclidqso.analysis import (
#     photometry, spectroscopy, redshift, templates, composite, qso_fitting
# )

# Visualization (will be available as modules are implemented)
# from euclidqso.visualization import plots, specplots, photplots

# External integrations (will be available as modules are implemented)
# from euclidqso.external import desi, wise, gaia, des

# Utilities (will be available as modules are implemented)  
# from euclidqso.utils import io, math, astro, validation

# Configuration
# from euclidqso import config

__all__ = [
    "__version__",
    # Core modules will be added as they are implemented
    # Analysis modules will be added as they are implemented
    # Visualization modules will be added as they are implemented
    # External modules will be added as they are implemented
    # Utility modules will be added as they are implemented
]

# Package metadata
__author__ = "Yuming Fu"
__email__ = "fuympku@outlook.com" 
__license__ = "GPLv3"
__url__ = "https://github.com/rudolffu/euclidqso"

# Set up logging
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())