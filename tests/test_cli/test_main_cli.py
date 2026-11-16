"""Tests for the top-level euclidkit CLI help."""

from click.testing import CliRunner

from euclidkit.cli.main import main


def test_top_level_help_uses_complete_command_summaries():
    """Top-level help should not ellipsize command descriptions."""
    result = CliRunner().invoke(main, ['--help'])

    assert result.exit_code == 0
    commands_section = result.output.split('Commands:', maxsplit=1)[1]
    assert '...' not in commands_section
    # Command names can grow; Click adjusts padding and wraps descriptions.
    normalised = ' '.join(commands_section.split())
    assert 'compile-segmap Cut out MER segmentation maps.' in normalised
    assert 'compile-spectra Export spectra to Parquet or FITS.' in normalised
    assert 'crossmatch Crossmatch sources with Euclid MER.' in normalised
    assert 'query-segmap Query MER segmentation-map metadata.' in normalised
    assert 'query-spectra Query spectra-source rows by object ID or position.' in normalised
    assert 'select-footprint Filter a catalog by Euclid footprint.' in normalised

    assert 'download-cutana-alternatives Download alternate Cutana mosaics.' in normalised


def test_diagnostics_finds_configured_credentials(monkeypatch, tmp_path):
    from types import SimpleNamespace
    from euclidkit.cli import main as main_module

    cred = tmp_path / 'cred.txt'
    cred.write_text('user\nsecret\n')
    monkeypatch.delenv('EUCLIDKIT_CREDENTIALS_FILE', raising=False)
    monkeypatch.setattr(main_module, 'config', SimpleNamespace(
        get=lambda key, default=None: str(cred) if key == 'data.credentials_file' else default))

    result = CliRunner().invoke(main, ['diagnostics', '--check-data'])

    assert result.exit_code == 0
    assert f'credentials: Credentials file found at {cred}' in result.output
    assert 'secret' not in result.output


def test_declared_console_scripts_resolve():
    import importlib
    import tomllib
    from pathlib import Path

    pyproject = Path(__file__).resolve().parents[2] / 'pyproject.toml'
    scripts = tomllib.loads(pyproject.read_text())['project']['scripts']
    for target in scripts.values():
        module, attr = target.split(':')
        assert callable(getattr(importlib.import_module(module), attr))
