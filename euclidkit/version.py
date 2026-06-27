"""Version information for euclidkit."""

from __future__ import annotations

import re
from importlib.metadata import PackageNotFoundError, version


def _get_version() -> str:
    """Return the installed package version, falling back for source checkouts."""
    try:
        from setuptools_scm import get_version

        return get_version(root="..", relative_to=__file__)
    except Exception:
        pass

    try:
        return version("euclidkit")
    except PackageNotFoundError:
        return "0+unknown"


__version__ = _get_version()

# Extract numeric release tuple from PEP 440 versions (e.g. 0.2.0rc1 -> (0, 2, 0)).
_release_match = re.match(r"^(\d+)\.(\d+)\.(\d+)", __version__)
if _release_match:
    __version_info__ = tuple(int(x) for x in _release_match.groups())
else:
    __version_info__ = (0, 0, 0)

# Version metadata
__author__ = "Yuming Fu"
__email__ = "fuympku@outlook.com"
__license__ = "BSD-3-Clause"
__copyright__ = "Copyright 2026, Yuming Fu"
