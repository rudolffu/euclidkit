"""Test CLI crossmatch functionality."""

import pytest
import tempfile
import os
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from click.testing import CliRunner

from astropy.table import Table
import numpy as np

from euclidqso.cli.crossmatch_cli import crossmatch, query_spectra, compile_spectra


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
                with patch('euclidqso.cli.crossmatch_cli.EuclidArchive') as mock_archive_class:
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
                with patch('euclidqso.cli.crossmatch_cli.EuclidArchive') as mock_archive_class:
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
                with patch('euclidqso.cli.crossmatch_cli.EuclidArchive') as mock_archive_class:
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
                with patch('euclidqso.cli.crossmatch_cli.EuclidArchive') as mock_archive_class:
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
                with patch('euclidqso.cli.crossmatch_cli.EuclidArchive') as mock_archive_class:
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

    def test_crossmatch_error_handling(self):
        """Test error handling in crossmatch command."""
        input_file = self.create_temp_input_file()
        
        try:
            with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
                with patch('euclidqso.cli.crossmatch_cli.EuclidArchive') as mock_archive_class:
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
                with patch('euclidqso.cli.crossmatch_cli.EuclidArchive') as mock_archive_class:
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
                with patch('euclidqso.cli.crossmatch_cli.EuclidArchive') as mock_archive_class:
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
                with patch('euclidqso.cli.crossmatch_cli.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive
                    
                    mock_spectra = Table({
                        'object_id': [100001, 100002],
                        'spectrum_id': ['spec_1', 'spec_2'],
                        'instrument_name': ['NISP', 'VIS']
                    })
                    mock_archive.query_spectra_sources.return_value = mock_spectra
                    
                    with patch('euclidqso.cli.crossmatch_cli.load_table', return_value=self.crossmatch_table):
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
                with patch('euclidqso.cli.crossmatch_cli.EuclidArchive') as mock_archive_class:
                    mock_archive = Mock()
                    mock_archive_class.return_value = mock_archive

                    mock_spectra = Table({
                        'object_id': [100001, 100002, 100003],
                        'spectrum_id': ['spec_1', 'spec_2', 'spec_3']
                    })
                    mock_archive.query_spectra_sources.return_value = mock_spectra

                    with patch('euclidqso.cli.crossmatch_cli.load_table', return_value=self.crossmatch_table):
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
            
            assert result.exit_code == 1
            assert "Must provide --crossmatch file" in result.output
            
        if os.path.exists(output_file.name):
            os.unlink(output_file.name)

    def test_query_spectra_instrument_breakdown(self):
        """Test instrument breakdown display."""
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as output_file:
            with patch('euclidqso.cli.crossmatch_cli.EuclidArchive') as mock_archive_class:
                mock_archive = Mock()
                mock_archive_class.return_value = mock_archive
                
                mock_spectra = Table({
                    'object_id': [100001, 100002, 100003, 100004],
                    'instrument_name': ['NISP', 'NISP', 'VIS', 'NISP'],
                    'spectrum_id': ['spec_1', 'spec_2', 'spec_3', 'spec_4']
                })
                mock_archive.query_spectra_sources.return_value = mock_spectra
                
                with patch('euclidqso.cli.crossmatch_cli.load_table', return_value=self.crossmatch_table):
                    result = self.runner.invoke(query_spectra, [
                        '--crossmatch', 'dummy.fits',
                        '--output', output_file.name,
                        '--verbose'
                    ])
                
                assert result.exit_code == 0
                assert "Spectra by instrument:" in result.output
                assert "NISP: 3" in result.output
                assert "VIS: 1" in result.output
                
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
                with patch('euclidqso.cli.crossmatch_cli.SpectrumCompiler') as mock_compiler_class:
                    mock_compiler = Mock()
                    mock_compiler_class.return_value = mock_compiler
                    
                    # Mock compilation results
                    output_files = [
                        f"{output_dir}/compiled_spectra_001.fits"
                    ]
                    mock_compiler.compile_spectra.return_value = output_files
                    mock_compiler.create_metadata_table.return_value = f"{output_dir}/compiled_spectra_metadata.fits"
                    
                    with patch('euclidqso.cli.crossmatch_cli.load_table', return_value=self.spectra_table):
                        result = self.runner.invoke(compile_spectra, [
                            '--spectra-table', spectra_file,
                            '--output-dir', output_dir,
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
                with patch('euclidqso.cli.crossmatch_cli.SpectrumCompiler') as mock_compiler_class:
                    mock_compiler = Mock()
                    mock_compiler_class.return_value = mock_compiler
                    mock_compiler.compile_spectra.return_value = []
                    mock_compiler.create_metadata_table.return_value = "metadata.fits"
                    
                    with patch('euclidqso.cli.crossmatch_cli.load_table', return_value=self.spectra_table):
                        result = self.runner.invoke(compile_spectra, [
                            '--spectra-table', spectra_file,
                            '--output-dir', output_dir,
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
                with patch('euclidqso.cli.crossmatch_cli.SpectrumCompiler') as mock_compiler_class:
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
                    
                    with patch('euclidqso.cli.crossmatch_cli.load_table', return_value=self.spectra_table):
                        result = self.runner.invoke(compile_spectra, [
                            '--spectra-table', spectra_file,
                            '--output-dir', output_dir
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
                with patch('euclidqso.cli.crossmatch_cli.SpectrumCompiler') as mock_compiler_class:
                    mock_compiler = Mock()
                    mock_compiler_class.return_value = mock_compiler
                    mock_compiler.compile_spectra.side_effect = Exception("Compilation failed")
                    
                    with patch('euclidqso.cli.crossmatch_cli.load_table', return_value=self.spectra_table):
                        result = self.runner.invoke(compile_spectra, [
                            '--spectra-table', spectra_file,
                            '--output-dir', output_dir
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
                with patch('euclidqso.cli.crossmatch_cli.SpectrumCompiler') as mock_compiler_class:
                    mock_compiler = Mock()
                    mock_compiler_class.return_value = mock_compiler
                    mock_compiler.compile_spectra.return_value = []
                    mock_compiler.create_metadata_table.return_value = "metadata.fits"
                    
                    with patch('euclidqso.cli.crossmatch_cli.load_table', return_value=self.spectra_table):
                        result = self.runner.invoke(compile_spectra, [
                            '--spectra-table', spectra_file,
                            '--output-dir', output_dir,
                            '--verbose'
                        ])
                    
                    assert result.exit_code == 0
                    assert "Loaded 3 spectral sources" in result.output
                    assert "Max extensions per file:" in result.output
                    assert "Output directory:" in result.output
                    
        finally:
            os.unlink(spectra_file)
