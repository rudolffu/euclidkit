"""Test version information."""

import euclidkit


def test_version_exists():
    """Test that version information is available."""
    assert hasattr(euclidkit, '__version__')
    assert isinstance(euclidkit.__version__, str)
    

def test_author_info():
    """Test that author information is available."""
    assert hasattr(euclidkit, '__author__')
    assert euclidkit.__author__ == "Yuming Fu"
    assert euclidkit.__email__ == "fuympku@outlook.com"