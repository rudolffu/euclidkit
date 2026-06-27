#!/usr/bin/env python3
"""Compatibility shim for tools that still invoke setup.py directly."""

import sys

from setuptools import setup


if sys.version_info < (3, 9):
    raise RuntimeError("euclidkit requires Python 3.9 or later")


if __name__ == "__main__":
    setup()
