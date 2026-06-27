"""Test CLI crossmatch functionality."""

import pytest
import tempfile
import os
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from click.testing import CliRunner

from astropy.table import Table
import numpy as np

from euclidkit.cli.crossmatch_cli import (
    crossmatch,
    query_spectra,
    query_zspe,
    query_segmap,
    compile_spectra,
    upload_table,
)


class TestCrossmatchCLI:
    """Test cases for crossmatch CLI command."""

    def setup_method(self):
        """Set up test fixtures."""
        self.runner = CliRunner()
        
        # Create sample input table
        self.sample_table = Table({
            'ra': [150.0, 151.0, 152.0],
            'dec': [2.0, 2.1, 2.2],
            'source_id': [1, 2, 3]
        })

    def create_temp_input_file(self):
        """Create a temporary input file."""
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        self.sample_table.write(temp_file.name, format='csv', overwrite=True)
        return temp_file.name

    def test_crossmatch_basic_usage(self):
        """Test basic crossmatch command usage."""
        input_file = self.create_temp_input_file()
        
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    # Mock the archive instance
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    
                    # Mock crossmatch results
                    mock_results = Table({
                        'object_id': [100001, 100002],
                        'separation_arcsec': [0.1, 0.2],
                        'ra': [150.001, 151.001],
                        'dec': [2.001, 2.101]
                    })
                    mock_archive.crossmatch_sources.return_value = mock_results
                    
                    result = self.runner.invoke(crossmatch, [
                        '--input', input_file,
                        '--output', output_file.name,
                        '--radius', '1.0',
                        '--environment', 'REG'
                    ])
                    
                    assert result.exit_code == 0
                    assert "Crossmatch completed: 2 matches found" in result.output
                    assert "Results saved to:" in result.output
                    
                    # Verify archive methods were called
                    mock_archive.login.assert_called_once()
                    mock_archive.crossmatch_sources.assert_called_once()
                    mock_archive.logout.assert_called_once()
                    
        finally:
            os.unlink(input_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_crossmatch_with_credentials(self):
        """Test crossmatch command with custom credentials."""
        input_file = self.create_temp_input_file()
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as cred_file:
            cred_file.write("username\npassword")
            cred_file.flush()
            
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_archive.crossmatch_sources.return_value = Table({'object_id': [1]})
                    
                    result = self.runner.invoke(crossmatch, [
                        '--input', input_file,
                        '--output', output_file.name,
                        '--credentials', cred_file.name
                    ])
                    
                    assert result.exit_code == 0
                    mock_archive.login.assert_called_once_with(credentials_file=cred_file.name)
                    
        finally:
            os.unlink(input_file)
            os.unlink(cred_file.name)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_crossmatch_with_custom_columns(self):
        """Test crossmatch command with custom RA/Dec columns."""
        # Create table with different column names
        custom_table = Table({
            'right_ascension': [150.0, 151.0],
            'declination': [2.0, 2.1],
            'id': [1, 2]
        })
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as input_file:
            custom_table.write(input_file.name, format='csv', overwrite=True)
            
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_archive.crossmatch_sources.return_value = Table({'object_id': [1]})
                    
                    result = self.runner.invoke(crossmatch, [
                        '--input', input_file.name,
                        '--output', output_file.name,
                        '--ra-col', 'right_ascension',
                        '--dec-col', 'declination'
                    ])
                    
                    assert result.exit_code == 0
                    
                    # Check that custom column names were passed
                    call_args = mock_archive.crossmatch_sources.call_args
                    assert call_args[1]['ra_col'] == 'right_ascension'
                    assert call_args[1]['dec_col'] == 'declination'
                    
        finally:
            os.unlink(input_file.name)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_crossmatch_with_max_sources(self):
        """Test crossmatch command with max sources limit."""
        input_file = self.create_temp_input_file()
        
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_archive.crossmatch_sources.return_value = Table({'object_id': [1]})
                    
                    result = self.runner.invoke(crossmatch, [
                        '--input', input_file,
                        '--output', output_file.name,
                        '--max-sources', '10'
                    ])
                    
                    assert result.exit_code == 0
                    
                    # Check that max_sources was passed
                    call_args = mock_archive.crossmatch_sources.call_args
                    assert call_args[1]['max_sources'] == 10
                    
        finally:
            os.unlink(input_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_crossmatch_verbose_output(self):
        """Test crossmatch command with verbose output."""
        input_file = self.create_temp_input_file()
        
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_results = Table({
                        'object_id': [100001],
                        'separation_arcsec': np.array([0.5])
                    })
                    mock_archive.crossmatch_sources.return_value = mock_results
                    
                    result = self.runner.invoke(crossmatch, [
                        '--input', input_file,
                        '--output', output_file.name,
                        '--verbose'
                    ])
                    
                    assert result.exit_code == 0
                    assert "Connected to PDR environment" in result.output
                    assert "Input table:" in result.output
                    assert "Search radius:" in result.output
                    assert "Separation statistics" in result.output
                    
        finally:
            os.unlink(input_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_crossmatch_full_async_option(self):
        """Full async flag should disable batching and propagate to archive."""
        input_file = self.create_temp_input_file()
        
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_archive.crossmatch_sources.return_value = {'job_id': 'ABC123'}
                    
                    result = self.runner.invoke(crossmatch, [
                        '--input', input_file,
                        '--output', output_file.name,
                        '--full-async',
                        '--verbose'
                    ])
                    
                    assert result.exit_code == 0
                    assert "Full-table async mode enabled" in result.output
                    assert "Crossmatch job submitted asynchronously." in result.output
                    assert "Job ID: ABC123" in result.output
                    assert "Job info saved to:" in result.output
                    call_args = mock_archive.crossmatch_sources.call_args
                    assert call_args[1]['full_async'] is True
                    
        finally:
            os.unlink(input_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_crossmatch_drop_empty_columns_option_local_input(self):
        """Drop-empty-columns flag should propagate for local input mode."""
        input_file = self.create_temp_input_file()
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_archive.crossmatch_sources.return_value = Table({'object_id': [1]})

                    result = self.runner.invoke(crossmatch, [
                        '--input', input_file,
                        '--output', output_file.name,
                        '--drop-empty-columns',
                    ])

                    assert result.exit_code == 0
                    call_args = mock_archive.crossmatch_sources.call_args
                    assert call_args.kwargs['drop_empty_columns'] is True
        finally:
            os.unlink(input_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_crossmatch_drop_empty_columns_option_user_table(self):
        """Drop-empty-columns flag should propagate for archive user-table mode."""
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
            try:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_archive.crossmatch_user_table.return_value = Table({'object_id': [1]})

                    result = self.runner.invoke(crossmatch, [
                        '--user-table-name', 'my_table',
                        '--output', output_file.name,
                        '--drop-empty-columns',
                    ])

                    assert result.exit_code == 0
                    call_args = mock_archive.crossmatch_user_table.call_args
                    assert call_args.kwargs['drop_empty_columns'] is True
            finally:
                if os.path.exists(output_file.name):
                    os.unlink(output_file.name)

    def test_upload_table_cli_job_submission(self):
        """Uploading a table should call archive helper and report job info."""
        input_file = self.create_temp_input_file()
        try:
            with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                mock_archive = Mock()
                mock_archive_class.return_value = mock_archive
                mock_archive.upload_user_table.return_value = {'job_id': 'job-77', 'format': 'csv'}

                result = self.runner.invoke(upload_table, [
                    '--input', input_file,
                    '--table-name', 'user_table',
                    '--description', 'demo',
                    '--format', 'csv',
                    '--overwrite',
                    '--environment', 'PDR'
                ])

                assert result.exit_code == 0
                assert "Upload job submitted (ID: job-77)." in result.output
                assert "Table name: user_table" in result.output
                mock_archive.upload_user_table.assert_called_once()
        finally:
            os.unlink(input_file)

    def test_upload_table_cli_success(self):
        """Synchronous uploads should report success."""
        input_file = self.create_temp_input_file()
        try:
            with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                mock_archive = Mock()
                mock_archive_class.return_value = mock_archive
                mock_archive.upload_user_table.return_value = {'job_id': None, 'format': 'fits'}

                result = self.runner.invoke(upload_table, [
                    '--input', input_file,
                    '--table-name', 'table_sync',
                    '--format', 'fits'
                ])

                assert result.exit_code == 0
                assert "Table uploaded successfully." in result.output
                assert "Format: fits" in result.output
        finally:
            os.unlink(input_file)
    def test_crossmatch_idr_wide_prefix(self):
        """IDR environment should prefix output filename and default to WIDE."""
        input_file = self.create_temp_input_file()
        
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_archive.crossmatch_sources.return_value = Table({'object_id': [1]})
                    
                    result = self.runner.invoke(crossmatch, [
                        '--input', input_file,
                        '--output', output_file.name,
                        '--environment', 'IDR'
                    ])
                    
                    assert result.exit_code == 0
                    call_args = mock_archive.crossmatch_sources.call_args
                    assert call_args[1]['idr_field'] == 'WIDE'
                    assert call_args[1]['idr_deep_partition'] == 'survey'
                    output_path = Path(call_args[1]['output_file'])
                    assert output_path.name.startswith('wide_')
                    assert f"Results saved to: {output_path}" in result.output
                    
        finally:
            os.unlink(input_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_crossmatch_idr_deep_selection(self):
        """IDR deep queries should use deep catalogue and prefix."""
        input_file = self.create_temp_input_file()
        
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_archive.crossmatch_sources.return_value = Table({'object_id': [1]})
                    
                    result = self.runner.invoke(crossmatch, [
                        '--input', input_file,
                        '--output', output_file.name,
                        '--environment', 'IDR',
                        '--idr-field', 'DEEP'
                    ])
                    
                    assert result.exit_code == 0
                    call_args = mock_archive.crossmatch_sources.call_args
                    assert call_args[1]['idr_field'] == 'DEEP'
                    assert call_args[1]['idr_deep_partition'] == 'survey'
                    output_path = Path(call_args[1]['output_file'])
                    assert output_path.name.startswith('deep_')
                    assert "Results saved to:" in result.output
                    
        finally:
            os.unlink(input_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_crossmatch_idr_deep_partition_option(self):
        """IDR deep partition option should be accepted and forwarded."""
        input_file = self.create_temp_input_file()

        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_archive.crossmatch_sources.return_value = Table({'object_id': [1]})

                    result = self.runner.invoke(crossmatch, [
                        '--input', input_file,
                        '--output', output_file.name,
                        '--environment', 'IDR',
                        '--idr-field', 'DEEP',
                        '--idr-deep-partition', 'mode',
                        '--verbose',
                    ])

                    assert result.exit_code == 0
                    call_args = mock_archive.crossmatch_sources.call_args
                    assert call_args[1]['idr_field'] == 'DEEP'
                    assert call_args[1]['idr_deep_partition'] == 'mode'
                    assert "IDR DEEP partition: mode" in result.output

        finally:
            os.unlink(input_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_crossmatch_error_handling(self):
        """Test error handling in crossmatch command."""
        input_file = self.create_temp_input_file()
        
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_archive.crossmatch_sources.side_effect = Exception("Connection failed")
                    
                    result = self.runner.invoke(crossmatch, [
                        '--input', input_file,
                        '--output', output_file.name
                    ])
                    
                    assert result.exit_code == 1
                    assert "Error in crossmatch: Connection failed" in result.output
                    mock_archive.logout.assert_called_once()
                    
        finally:
            os.unlink(input_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_crossmatch_match_mode_object_id(self):
        """Test crossmatch with object-id match mode and no separation stats."""
        input_file = self.create_temp_input_file()
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_results = Table({
                        'object_id': [100001, 100002],
                    })
                    mock_archive.crossmatch_sources.return_value = mock_results

                    result = self.runner.invoke(crossmatch, [
                        '--input', input_file,
                        '--output', output_file.name,
                        '--match-mode', 'object-id'
                    ])

                    assert result.exit_code == 0
                    # Ensure flag propagated
                    call_args = mock_archive.crossmatch_sources.call_args
                    assert call_args[1]['use_object_id'] is True
                    # No separation stats printed
                    assert 'Separation statistics' not in result.output
        finally:
            os.unlink(input_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_crossmatch_match_mode_spatial(self):
        """Test crossmatch with spatial match mode forces spatial path."""
        input_file = self.create_temp_input_file()
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_results = Table({
                        'object_id': [100001],
                        'separation_arcsec': [0.3],
                    })
                    mock_archive.crossmatch_sources.return_value = mock_results

                    result = self.runner.invoke(crossmatch, [
                        '--input', input_file,
                        '--output', output_file.name,
                        '--match-mode', 'spatial'
                    ])

                    assert result.exit_code == 0
                    call_args = mock_archive.crossmatch_sources.call_args
                    assert call_args[1]['use_object_id'] is False
                    # Separation stats present in output in this path
                    assert 'Separation statistics' in result.output
        finally:
            os.unlink(input_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)


class TestQuerySpectraCLI:
    """Test cases for query-spectra CLI command."""

    def setup_method(self):
        """Set up test fixtures."""
        self.runner = CliRunner()
        
        # Sample crossmatch table
        self.crossmatch_table = Table({
            'object_id': [100001, 100002, 100003],
            'ra': [150.0, 151.0, 152.0],
            'dec': [2.0, 2.1, 2.2]
        })

    def create_temp_crossmatch_file(self):
        """Create a temporary crossmatch file."""
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fits', delete=False)
        self.crossmatch_table.write(temp_file.name, format='fits', overwrite=True)
        return temp_file.name

    def test_query_spectra_with_crossmatch_file(self):
        """Test query-spectra with crossmatch file input."""
        crossmatch_file = self.create_temp_crossmatch_file()
        
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    
                    mock_spectra = Table({
                        'object_id': [100001, 100002],
                        'spectrum_id': ['spec_1', 'spec_2'],
                        'instrument_name': ['NISP', 'VIS']
                    })
                    mock_archive.query_spectra_sources.return_value = mock_spectra
                    
                    with patch('euclidkit.utils.io.load_table', return_value=self.crossmatch_table):
                        result = self.runner.invoke(query_spectra, [
                            '--crossmatch', crossmatch_file,
                            '--output', output_file.name
                        ])
                    
                    assert result.exit_code == 0
                    assert "Spectral query completed: 2 spectra found" in result.output
                    assert "Unique objects with spectra: 2" in result.output
                    
        finally:
            os.unlink(crossmatch_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_query_spectra_with_crossmatch_only(self):
        """Test query-spectra requires crossmatch and loads it."""
        crossmatch_file = self.create_temp_crossmatch_file()
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive

                    mock_spectra = Table({
                        'object_id': [100001, 100002, 100003],
                        'spectrum_id': ['spec_1', 'spec_2', 'spec_3']
                    })
                    mock_archive.query_spectra_sources.return_value = mock_spectra

                    with patch('euclidkit.utils.io.load_table', return_value=self.crossmatch_table):
                        result = self.runner.invoke(query_spectra, [
                            '--crossmatch', crossmatch_file,
                            '--output', output_file.name
                        ])

                    assert result.exit_code == 0
                    assert "Spectral query completed: 3 spectra found" in result.output
        finally:
            os.unlink(crossmatch_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_query_spectra_no_input_error(self):
        """Test error when no input provided."""
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
            result = self.runner.invoke(query_spectra, [
                '--output', output_file.name
            ])
            
            assert result.exit_code == 2
            assert "Missing option '--crossmatch' / '-x'" in result.output
            
        if os.path.exists(output_file.name):
            os.unlink(output_file.name)

    def test_query_spectra_instrument_breakdown(self):
        """Test instrument breakdown display."""
        crossmatch_file = self.create_temp_crossmatch_file()
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
            try:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive

                    mock_spectra = Table({
                        'object_id': [100001, 100002, 100003, 100004],
                        'instrument_name': ['NISP', 'NISP', 'VIS', 'NISP'],
                        'spectrum_id': ['spec_1', 'spec_2', 'spec_3', 'spec_4']
                    })
                    mock_archive.query_spectra_sources.return_value = mock_spectra

                    with patch('euclidkit.utils.io.load_table', return_value=self.crossmatch_table):
                        result = self.runner.invoke(query_spectra, [
                            '--crossmatch', crossmatch_file,
                            '--output', output_file.name,
                            '--verbose'
                        ])

                    assert result.exit_code == 0
                    assert "Spectra by instrument:" in result.output
                    assert "NISP: 3" in result.output
                    assert "VIS: 1" in result.output
            finally:
                os.unlink(crossmatch_file)
                if os.path.exists(output_file.name):
                    os.unlink(output_file.name)


class TestQueryZspeCLI:
    """Test cases for query-zspe CLI command."""

    def setup_method(self):
        """Set up test fixtures."""
        self.runner = CliRunner()
        self.crossmatch_table = Table({
            'object_id': [100001, 100002, 100003],
        })

    def create_temp_crossmatch_file(self):
        """Create a temporary crossmatch file."""
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fits', delete=False)
        self.crossmatch_table.write(temp_file.name, format='fits', overwrite=True)
        return temp_file.name

    def test_query_zspe_defaults_to_qso_wide_idr(self):
        """query-zspe should default to IDR QSO WIDE candidate lookup."""
        crossmatch_file = self.create_temp_crossmatch_file()
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_archive.query_zspe_candidates.return_value = Table({
                        'object_id': [100001, 100002],
                        'spe_rank': [1, 1],
                        'spe_z': [2.1, 0.8],
                        'spe_z_err': [0.01, 0.02],
                        'source_table': ['wide_survey', 'wide'],
                    })

                    with patch('euclidkit.utils.io.load_table', return_value=self.crossmatch_table):
                        result = self.runner.invoke(query_zspe, [
                            '--crossmatch', crossmatch_file,
                            '--output', output_file.name,
                        ])

                    assert result.exit_code == 0
                    mock_archive_class.assert_called_once_with(environment='IDR')
                    mock_archive.query_zspe_candidates.assert_called_once()
                    kwargs = mock_archive.query_zspe_candidates.call_args.kwargs
                    assert kwargs['object_type'] == 'qso'
                    assert kwargs['idr_field'] == 'WIDE'
                    assert kwargs['full_async'] is False
                    assert "SPE redshift query completed: 2 candidates found" in result.output
                    assert "Unique objects with SPE redshifts: 2" in result.output
        finally:
            os.unlink(crossmatch_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_query_zspe_full_async_handles_metadata(self):
        """query-zspe --full-async should pass async options and print metadata."""
        crossmatch_file = self.create_temp_crossmatch_file()
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_archive.query_zspe_candidates.return_value = {
                        'results_downloaded': True,
                        'result_row_count': 3,
                        'chunk_count': 2,
                        'chunk_size': 2,
                        'job_id': 'zs-job',
                    }

                    with patch('euclidkit.utils.io.load_table', return_value=self.crossmatch_table):
                        result = self.runner.invoke(query_zspe, [
                            '--crossmatch', crossmatch_file,
                            '--output', output_file.name,
                            '--full-async',
                            '--async-chunk-size', '2',
                        ])

                    assert result.exit_code == 0
                    kwargs = mock_archive.query_zspe_candidates.call_args.kwargs
                    assert kwargs['full_async'] is True
                    assert kwargs['async_chunk_size'] == 2
                    assert "SPE redshift async query completed and results were downloaded." in result.output
                    assert "Rows downloaded: 3" in result.output
                    assert "Chunks processed: 2 (chunk size: 2)" in result.output
        finally:
            os.unlink(crossmatch_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)


class TestQuerySegmapCLI:
    """Test cases for query-segmap CLI command."""

    def setup_method(self):
        """Set up test fixtures."""
        self.runner = CliRunner()
        self.input_table = Table({
            'object_id': [100001, 100002],
            'SEGMENTATION_MAP_ID': [1234567890000, 9876543210000],
            'ra': [150.0, 151.0],
            'dec': [2.0, 2.1],
        })

    def create_temp_input_file(self):
        """Create a temporary query-segmap input file."""
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fits', delete=False)
        self.input_table.write(temp_file.name, format='fits', overwrite=True)
        return temp_file.name

    def test_query_segmap_defaults_to_pdr_and_prints_counts(self):
        """query-segmap should load input, query archive, and print row-count details."""
        input_file = self.create_temp_input_file()
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive._last_segmap_valid_tile_count = 2
                    mock_archive_class.return_value = mock_archive
                    mock_archive.query_segmentation_maps.return_value = Table({
                        'object_id': [100001],
                        'segmentation_map_id': [1234567890000],
                    })

                    with patch('euclidkit.utils.io.load_table', return_value=self.input_table):
                        result = self.runner.invoke(query_segmap, [
                            '--input', input_file,
                            '--output', output_file.name,
                        ])

                    assert result.exit_code == 0
                    mock_archive_class.assert_called_once_with(environment='PDR')
                    mock_archive.login.assert_called_once()
                    mock_archive.query_segmentation_maps.assert_called_once()
                    kwargs = mock_archive.query_segmentation_maps.call_args.kwargs
                    assert kwargs['source_table'] is self.input_table
                    assert kwargs['output_file'] == output_file.name
                    assert "Segmentation map query completed: 1 matches found" in result.output
                    assert "Input rows: 2" in result.output
                    assert "Valid tile rows: 2" in result.output
        finally:
            os.unlink(input_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_query_segmap_accepts_idr_without_field_options(self):
        """query-segmap -e IDR should forward only the environment selection."""
        input_file = self.create_temp_input_file()
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive._last_segmap_valid_tile_count = 2
                    mock_archive_class.return_value = mock_archive
                    mock_archive.query_segmentation_maps.return_value = Table()

                    with patch('euclidkit.utils.io.load_table', return_value=self.input_table):
                        result = self.runner.invoke(query_segmap, [
                            '--input', input_file,
                            '--output', output_file.name,
                            '--environment', 'IDR',
                        ])

                    assert result.exit_code == 0
                    mock_archive_class.assert_called_once_with(environment='IDR')
                    call_kwargs = mock_archive.query_segmentation_maps.call_args.kwargs
                    assert 'idr_field' not in call_kwargs
                    assert 'idr_deep_partition' not in call_kwargs
        finally:
            os.unlink(input_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)

    def test_query_segmap_surfaces_missing_segmentation_id(self):
        """Missing SEGMENTATION_MAP_ID should be reported as a CLI error."""
        input_file = self.create_temp_input_file()
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidkit.core.data_access.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    mock_archive.query_segmentation_maps.side_effect = ValueError(
                        "Input table must contain SEGMENTATION_MAP_ID. Run euclidkit crossmatch first "
                        "to add segmentation_map_id before query-segmap."
                    )

                    with patch('euclidkit.utils.io.load_table', return_value=self.input_table):
                        result = self.runner.invoke(query_segmap, [
                            '--input', input_file,
                            '--output', output_file.name,
                        ])

                    assert result.exit_code != 0
                    assert "Run euclidkit crossmatch first" in result.output
        finally:
            os.unlink(input_file)
            if os.path.exists(output_file.name):
                os.unlink(output_file.name)


class TestCompileSpectraCLI:
    """Test cases for compile-spectra CLI command."""

    def setup_method(self):
        """Set up test fixtures."""
        self.runner = CliRunner()
        
        # Sample spectra table
        self.spectra_table = Table({
            'spectrum_id': ['spec_1', 'spec_2', 'spec_3'],
            'object_id': [100001, 100002, 100003],
            'file_path': ['/data/spec1.fits', '/data/spec2.fits', '/data/spec3.fits'],
            'instrument_name': ['NISP', 'VIS', 'NISP']
        })

    def create_temp_spectra_file(self):
        """Create a temporary spectra file."""
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fits', delete=False)
        self.spectra_table.write(temp_file.name, format='fits', overwrite=True)
        return temp_file.name

    def test_compile_spectra_basic_usage(self):
        """Test basic compile-spectra command usage."""
        spectra_file = self.create_temp_spectra_file()
        
        try:
            with tempfile.TemporaryDirectory() as output_dir:
                with patch('euclidkit.core.spectra.SpectrumCompiler') as mock_compiler_class:
                    mock_compiler = Mock()
                    mock_compiler_class.return_value = mock_compiler
                    
                    # Mock compilation results
                    output_files = [
                        f"{output_dir}/compiled_spectra_001.fits"
                    ]
                    mock_compiler.compile_spectra.return_value = output_files
                    mock_compiler.create_metadata_table.return_value = f"{output_dir}/compiled_spectra_metadata.fits"
                    
                    with patch('euclidkit.utils.io.load_table', return_value=self.spectra_table):
                        result = self.runner.invoke(compile_spectra, [
                            '--spectra-table', spectra_file,
                            '--output-dir', output_dir,
                            '--output-format', 'fits',
                            '--prefix', 'test_compiled'
                        ])
                    
                    assert result.exit_code == 0
                    assert "Compilation completed successfully!" in result.output
                    assert "Created 1 FITS files:" in result.output
                    assert "Total spectra processed: 3" in result.output
                    
                    # Verify compiler methods were called
                    mock_compiler.compile_spectra.assert_called_once()
                    mock_compiler.create_metadata_table.assert_called_once()
                    
        finally:
            os.unlink(spectra_file)

    def test_compile_spectra_with_options(self):
        """Test compile-spectra with custom options."""
        spectra_file = self.create_temp_spectra_file()
        
        try:
            with tempfile.TemporaryDirectory() as output_dir:
                with patch('euclidkit.core.spectra.SpectrumCompiler') as mock_compiler_class:
                    mock_compiler = Mock()
                    mock_compiler_class.return_value = mock_compiler
                    mock_compiler.compile_spectra.return_value = []
                    mock_compiler.create_metadata_table.return_value = "metadata.fits"
                    
                    with patch('euclidkit.utils.io.load_table', return_value=self.spectra_table):
                        result = self.runner.invoke(compile_spectra, [
                            '--spectra-table', spectra_file,
                            '--output-dir', output_dir,
                            '--output-format', 'fits',
                            '--prefix', 'custom_prefix',
                            '--max-extensions', '1000',
                            '--overwrite',
                            '--verbose'
                        ])
                    
                    assert result.exit_code == 0
                    
                    # Check that options were passed correctly
                    compile_call_args = mock_compiler.compile_spectra.call_args
                    assert compile_call_args[1]['output_prefix'] == 'custom_prefix'
                    assert compile_call_args[1]['overwrite'] is True
                    
                    # Check compiler initialization
                    mock_compiler_class.assert_called_once_with(max_extensions=1000)
                    
        finally:
            os.unlink(spectra_file)


    def test_compile_spectra_multiple_files(self):
        """Test compile-spectra creating multiple output files."""
        spectra_file = self.create_temp_spectra_file()
        
        try:
            with tempfile.TemporaryDirectory() as output_dir:
                with patch('euclidkit.core.spectra.SpectrumCompiler') as mock_compiler_class:
                    mock_compiler = Mock()
                    mock_compiler_class.return_value = mock_compiler
                    
                    # Mock multiple output files
                    output_files = [
                        f"{output_dir}/compiled_001.fits",
                        f"{output_dir}/compiled_002.fits",
                        f"{output_dir}/compiled_003.fits"
                    ]
                    mock_compiler.compile_spectra.return_value = output_files
                    mock_compiler.create_metadata_table.return_value = f"{output_dir}/metadata.fits"
                    
                    with patch('euclidkit.utils.io.load_table', return_value=self.spectra_table):
                        result = self.runner.invoke(compile_spectra, [
                            '--spectra-table', spectra_file,
                            '--output-dir', output_dir,
                            '--output-format', 'fits'
                        ])
                    
                    assert result.exit_code == 0
                    assert "Created 3 FITS files:" in result.output
                    assert "1. compiled_001.fits" in result.output
                    assert "2. compiled_002.fits" in result.output
                    assert "3. compiled_003.fits" in result.output
                    
        finally:
            os.unlink(spectra_file)

    def test_compile_spectra_error_handling(self):
        """Test error handling in compile-spectra command."""
        spectra_file = self.create_temp_spectra_file()
        
        try:
            with tempfile.TemporaryDirectory() as output_dir:
                with patch('euclidkit.core.spectra.SpectrumCompiler') as mock_compiler_class:
                    mock_compiler = Mock()
                    mock_compiler_class.return_value = mock_compiler
                    mock_compiler.compile_spectra.side_effect = Exception("Compilation failed")
                    
                    with patch('euclidkit.utils.io.load_table', return_value=self.spectra_table):
                        result = self.runner.invoke(compile_spectra, [
                            '--spectra-table', spectra_file,
                            '--output-dir', output_dir,
                            '--output-format', 'fits'
                        ])
                    
                    assert result.exit_code == 1
                    assert "Error compiling spectra: Compilation failed" in result.output
                    
        finally:
            os.unlink(spectra_file)

    def test_compile_spectra_verbose_output(self):
        """Test verbose output in compile-spectra."""
        spectra_file = self.create_temp_spectra_file()
        
        try:
            with tempfile.TemporaryDirectory() as output_dir:
                with patch('euclidkit.core.spectra.SpectrumCompiler') as mock_compiler_class:
                    mock_compiler = Mock()
                    mock_compiler_class.return_value = mock_compiler
                    mock_compiler.compile_spectra.return_value = []
                    mock_compiler.create_metadata_table.return_value = "metadata.fits"
                    
                    with patch('euclidkit.utils.io.load_table', return_value=self.spectra_table):
                        result = self.runner.invoke(compile_spectra, [
                            '--spectra-table', spectra_file,
                            '--output-dir', output_dir,
                            '--output-format', 'fits',
                            '--verbose'
                        ])
                    
                    assert result.exit_code == 0
                    assert "Loaded 3 spectral sources" in result.output
                    assert "Max extensions per file:" in result.output
                    assert "Output directory:" in result.output
                    
        finally:
            os.unlink(spectra_file)
