#!/usr/bin/env python3
"""
Setup script for euclidkit.
This is a fallback for older pip versions that don't support pyproject.toml.
"""

from setuptools import setup, find_packages
import os
import sys

# Ensure Python 3.9+
if sys.version_info < (3, 9):
    raise RuntimeError("euclidkit requires Python 3.9 or later")

# Read version from euclidkit/version.py
version_file = os.path.join(os.path.dirname(__file__), "euclidkit", "version.py")
version_info = {}
if os.path.exists(version_file):
    with open(version_file) as f:
        exec(f.read(), version_info)
    version = version_info["__version__"]
else:
    version = "0.1.0"

# Read README
readme_file = os.path.join(os.path.dirname(__file__), "README.md")
if os.path.exists(readme_file):
    with open(readme_file, "r", encoding="utf-8") as f:
        long_description = f.read()
else:
    long_description = "Euclid QSO Analysis Package"

# Core dependencies
install_requires = [
    "numpy>=1.19.0",
    "scipy>=1.6.0",
    "matplotlib>=3.3.0",
    "astropy>=4.0.0",
    "mocpy>=0.13.0",
    "astroquery>=0.4.0",
    "pandas>=1.2.0",
    "photutils>=1.0.0",
    "sep>=1.2.0",
    "scikit-image>=0.18.0",
    "tqdm>=4.50.0",
    "pyyaml>=5.4.0",
    "click>=8.0.0",
    "requests>=2.25.0",
]

# Optional dependencies
extras_require = {
    "desi": [
        "sparcl-client>=1.0.0",
    ],
    "visualization": [
        "jdaviz>=2.0.0",
        "bokeh>=2.0.0",
        "plotly>=5.0.0",
    ],
    "dev": [
        "pytest>=6.0.0",
        "pytest-cov>=2.10.0",
        "black>=21.0.0",
        "isort>=5.0.0",
        "flake8>=3.8.0",
        "mypy>=0.800",
        "pre-commit>=2.10.0",
        "sphinx>=4.0.0",
        "sphinx-rtd-theme>=0.5.0",
        "jupyter>=1.0.0",
        "ipython>=7.20.0",
    ],
}

# Add 'complete' option
extras_require["complete"] = [
    dep for extra in ["desi", "visualization", "dev"] 
    for dep in extras_require[extra]
]

# Entry points for command line scripts
entry_points = {
    "console_scripts": [
        "euclidkit=euclidkit.cli.main:main",
        "euclidkit-cutouts=euclidkit.cli.cutout_cli:cutouts",
        "euclidkit-spectra=euclidkit.cli.spec_cli:main", 
        "euclidkit-composite=euclidkit.cli.composite_cli:main",
        "euclidkit-pipeline=euclidkit.cli.main:pipeline",
    ],
}

# Package data
package_data = {
    "euclidkit": [
        "data/templates/*",
        "data/filters/*",
        "data/config/*",
        "data/mocs/*",
    ],
}

if __name__ == "__main__":
    setup(
        name="euclidkit",
        version=version,
        description="Euclid QSO Analysis Package",
        long_description=long_description,
        long_description_content_type="text/markdown",
        author="Yuming Fu",
        author_email="fuympku@outlook.com",
        maintainer="Yuming Fu", 
        maintainer_email="fuympku@outlook.com",
        url="https://github.com/rudolffu/euclidkit",
        project_urls={
            "Documentation": "https://euclidkit.readthedocs.io",
            "Source": "https://github.com/rudolffu/euclidkit",
            "Tracker": "https://github.com/rudolffu/euclidkit/issues",
        },
        packages=find_packages(),
        include_package_data=True,
        package_data=package_data,
        install_requires=install_requires,
        extras_require=extras_require,
        entry_points=entry_points,
        python_requires=">=3.9",
        classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Astronomy",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
            "Operating System :: OS Independent",
        ],
        keywords="astronomy euclid qso spectroscopy photometry",
        zip_safe=False,
    )
