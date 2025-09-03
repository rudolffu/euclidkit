"""Test version information."""

import euclidqso


def test_version_exists():
    """Test that version information is available."""
    assert hasattr(euclidqso, '__version__')
    assert isinstance(euclidqso.__version__, str)
    

def test_author_info():
    """Test that author information is available."""
    assert hasattr(euclidqso, '__author__')
    assert euclidqso.__author__ == "Yuming Fu"
    assert euclidqso.__email__ == "fuympku@outlook.com"