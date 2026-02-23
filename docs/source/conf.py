"""Sphinx configuration for euclidkit documentation."""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Ensure package import works during docs build.
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

project = "euclidkit"
author = "Yuming Fu"
copyright = "2026, Yuming Fu"

try:
    from euclidkit.version import __version__ as release
except Exception:
    release = "0.0.0"

version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = os.environ.get("SPHINX_THEME", "sphinx_rtd_theme")
html_static_path = ["_static"]

master_doc = "index"
