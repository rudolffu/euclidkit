"""Explicit fixtures for tests that must not read personal credentials."""
from types import SimpleNamespace

import pytest


@pytest.fixture
def isolated_credentials(monkeypatch, tmp_path):
    from euclidkit.core import data_access

    for name in ['EUCLID_USER', 'EUCLID_PASSWORD', 'EUCLIDKIT_CREDENTIALS_FILE']:
        monkeypatch.delenv(name, raising=False)
    config = SimpleNamespace(get=lambda key, default=None: (
        str(tmp_path/'absent-credentials.txt') if key == 'data.credentials_file' else default))
    monkeypatch.setattr(data_access, 'config', config)
    return config
