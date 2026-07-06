"""Tests for the top-level euclidkit CLI help."""

from click.testing import CliRunner

from euclidkit.cli.main import main


def test_top_level_help_uses_complete_command_summaries():
    """Top-level help should not ellipsize command descriptions."""
    result = CliRunner().invoke(main, ['--help'])

    assert result.exit_code == 0
    commands_section = result.output.split('Commands:', maxsplit=1)[1]
    assert '...' not in commands_section
    assert 'compile-segmap      Cut out MER segmentation maps.' in result.output
    assert 'compile-spectra     Export spectra to Parquet or FITS.' in result.output
    assert 'crossmatch          Crossmatch sources with Euclid MER.' in result.output
    assert 'query-segmap        Query MER segmentation-map metadata.' in result.output
    assert 'query-spectra       Query spectra-source rows by object ID or position.' in result.output
    assert 'select-footprint    Filter a catalog by Euclid footprint.' in result.output
