"""Tests for compile-spectra datalink mode."""

import os
import tempfile
from unittest.mock import Mock, patch

from astropy.table import Table
from click.testing import CliRunner

from euclidkit.cli.crossmatch_cli import compile_spectra


def _make_temp_spectra_file(table: Table) -> str:
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.fits', delete=False)
    table.write(tmp.name, format='fits', overwrite=True)
    return tmp.name


def test_compile_spectra_use_datalink_mode():
    """compile-spectra should use datalink path and archive login when requested."""
    runner = CliRunner()
    spectra_table = Table({
        'source_id': [100001, 100002],
        'object_id': [100001, 100002],
    })
    spectra_file = _make_temp_spectra_file(spectra_table)

    try:
        with tempfile.TemporaryDirectory() as output_dir:
            mock_compiler = Mock()
            mock_compiler.compile_spectra_datalink.return_value = [f"{output_dir}/dl_compiled_chunk_001.fits"]
            mock_compiler.create_metadata_table.return_value = f"{output_dir}/dl_compiled_metadata.fits"

            mock_archive = Mock()
            mock_archive.euclid = Mock()

            with patch('euclidkit.core.spectra.SpectrumCompiler', return_value=mock_compiler):
                with patch('euclidkit.core.data_access.EuclidArchive', return_value=mock_archive):
                    with patch('euclidkit.utils.io.load_table', return_value=spectra_table):
                        result = runner.invoke(compile_spectra, [
                            '--spectra-table', spectra_file,
                            '--output-dir', output_dir,
                            '--prefix', 'dl_compiled',
                            '--use-datalink',
                            '--environment', 'IDR',
                            '--retrieval-type', 'SPECTRA_BGS',
                            '--schema', 'dr1',
                        ])

            assert result.exit_code == 0
            mock_archive.login.assert_called_once()
            mock_archive.logout.assert_called_once()
            mock_compiler.compile_spectra_datalink.assert_called_once()
            kwargs = mock_compiler.compile_spectra_datalink.call_args.kwargs
            assert kwargs['retrieval_type'] == 'SPECTRA_BGS'
            assert kwargs['schema'] == 'dr1'
            assert kwargs['euclid_client'] is mock_archive.euclid
    finally:
        if os.path.exists(spectra_file):
            os.unlink(spectra_file)
