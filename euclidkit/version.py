"""Version information for euclidkit package."""

import re

__version__ = "0.2.2"

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
