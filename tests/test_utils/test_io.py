"""Test I/O utility functionality."""

import pytest
import numpy as np
import pandas as pd
import tempfile
import os
from pathlib import Path

from astropy.table import Table
from astropy.io import fits

from euclidkit.utils.io import (
    load_table, save_table, DataLoader, FileManager, FormatConverter
)


class TestLoadSaveTable:
    """Test basic table loading and saving functions."""

    def setup_method(self):
        """Set up test fixtures."""
        self.sample_table = Table({
            'ra': [150.0, 151.0, 152.0],
            'dec': [2.0, 2.1, 2.2],
            'mag': [20.5, 21.0, 19.8],
            'object_id': [1001, 1002, 1003]
        })

    def test_load_table_csv(self):
        """Test loading CSV files."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            self.sample_table.write(f.name, format='csv', overwrite=True)
            
        try:
            loaded_table = load_table(f.name)
            assert isinstance(loaded_table, Table)
            assert len(loaded_table) == 3
            assert 'ra' in loaded_table.colnames
        finally:
            os.unlink(f.name)

    def test_load_table_fits(self):
        """Test loading FITS files."""
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as f:
            self.sample_table.write(f.name, format='fits', overwrite=True)
            
        try:
            loaded_table = load_table(f.name)
            assert isinstance(loaded_table, Table)
            assert len(loaded_table) == 3
        finally:
            os.unlink(f.name)

    def test_load_table_parquet(self):
        """Test loading Parquet files."""
        with tempfile.NamedTemporaryFile(suffix='.parquet', delete=False) as f:
            self.sample_table.to_pandas().to_parquet(f.name)

        try:
            loaded_table = load_table(f.name)
            assert isinstance(loaded_table, Table)
            assert len(loaded_table) == 3
            assert 'ra' in loaded_table.colnames
        finally:
            os.unlink(f.name)

    def test_load_table_votable(self):
        """Test loading VOTable files."""
        with tempfile.NamedTemporaryFile(suffix='.xml', delete=False) as f:
            self.sample_table.write(f.name, format='votable', overwrite=True)
            
        try:
            loaded_table = load_table(f.name)
            assert isinstance(loaded_table, Table)
            assert len(loaded_table) == 3
        finally:
            os.unlink(f.name)

    def test_load_table_auto_format_detection(self):
        """Test automatic format detection."""
        # Test with different extensions
        formats_and_extensions = [
            ('csv', '.csv'),
            ('fits', '.fits'),
            ('votable', '.xml')
        ]
        
        for fmt, ext in formats_and_extensions:
            with tempfile.NamedTemporaryFile(suffix=ext, delete=False) as f:
                self.sample_table.write(f.name, format=fmt, overwrite=True)
                
            try:
                loaded_table = load_table(f.name)  # No format specified
                assert isinstance(loaded_table, Table)
                assert len(loaded_table) == 3
            finally:
                os.unlink(f.name)

        with tempfile.NamedTemporaryFile(suffix='.parquet', delete=False) as f:
            self.sample_table.to_pandas().to_parquet(f.name)

        try:
            loaded_table = load_table(f.name)
            assert isinstance(loaded_table, Table)
            assert len(loaded_table) == 3
        finally:
            os.unlink(f.name)

    def test_save_table_csv(self):
        """Test saving to CSV format."""
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as f:
            output_path = f.name
            
        try:
            save_table(self.sample_table, output_path, format='csv')
            
            # Verify file was created and can be loaded
            loaded_table = load_table(output_path)
            assert len(loaded_table) == 3
        finally:
            if os.path.exists(output_path):
                os.unlink(output_path)

    def test_save_table_auto_format(self):
        """Test automatic format detection for saving."""
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as f:
            output_path = f.name
            
        try:
            save_table(self.sample_table, output_path)  # No format specified
            
            # Verify file was created
            assert os.path.exists(output_path)
            loaded_table = load_table(output_path)
            assert len(loaded_table) == 3
        finally:
            if os.path.exists(output_path):
                os.unlink(output_path)

    def test_save_table_parquet(self):
        """Test saving to Parquet format."""
        with tempfile.NamedTemporaryFile(suffix='.parquet', delete=False) as f:
            output_path = f.name

        try:
            save_table(self.sample_table, output_path)
            loaded_table = load_table(output_path)
            assert len(loaded_table) == 3
            assert 'object_id' in loaded_table.colnames
        finally:
            if os.path.exists(output_path):
                os.unlink(output_path)

    def test_load_table_invalid_format(self):
        """Test error handling for invalid format."""
        with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as f:
            f.write("not a valid table format")
            f.flush()
            
        try:
            with pytest.raises(ValueError, match="Could not determine file format"):
                load_table(f.name)
        finally:
            os.unlink(f.name)

    def test_load_table_nonexistent_file(self):
        """Test error handling for non-existent files."""
        with pytest.raises(FileNotFoundError):
            load_table('/nonexistent/file.csv')


class TestDataLoader:
    """Test DataLoader class functionality."""

    def setup_method(self):
        """Set up test fixtures."""
        self.loader = DataLoader()
        self.sample_data = {
            'ra': [150.0, 151.0, 152.0],
            'dec': [2.0, 2.1, 2.2],
            'mag': [20.5, 21.0, 19.8]
        }

    def test_init(self):
        """Test initialization."""
        assert self.loader.cache_enabled is True
        assert self.loader.cache_dir is None

    def test_init_with_cache_dir(self):
        """Test initialization with cache directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            loader = DataLoader(cache_dir=temp_dir, cache_enabled=False)
            assert loader.cache_dir == Path(temp_dir)
            assert loader.cache_enabled is False

    def test_load_from_dict(self):
        """Test loading from dictionary."""
        table = self.loader.load_from_dict(self.sample_data)
        assert isinstance(table, Table)
        assert len(table) == 3
        assert 'ra' in table.colnames

    def test_load_from_pandas(self):
        """Test loading from pandas DataFrame."""
        df = pd.DataFrame(self.sample_data)
        table = self.loader.load_from_pandas(df)
        assert isinstance(table, Table)
        assert len(table) == 3
        assert 'ra' in table.colnames

    def test_load_from_numpy_arrays(self):
        """Test loading from numpy arrays."""
        arrays = {
            'ra': np.array([150.0, 151.0, 152.0]),
            'dec': np.array([2.0, 2.1, 2.2]),
            'mag': np.array([20.5, 21.0, 19.8])
        }
        
        table = self.loader.load_from_arrays(arrays)
        assert isinstance(table, Table)
        assert len(table) == 3
        assert 'ra' in table.colnames


class TestFileManager:
    """Test FileManager class functionality."""

    def setup_method(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.manager = FileManager(base_dir=self.temp_dir)

    def teardown_method(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_init(self):
        """Test initialization."""
        assert self.manager.base_dir == Path(self.temp_dir)

    def test_create_directory(self):
        """Test directory creation."""
        subdir = self.manager.create_directory('test_subdir')
        assert subdir.exists()
        assert subdir.is_dir()

    def test_get_temp_file(self):
        """Test temporary file creation."""
        temp_file = self.manager.get_temp_file(suffix='.fits')
        assert temp_file.suffix == '.fits'
        assert temp_file.parent == self.manager.base_dir

    def test_cleanup_temp_files(self):
        """Test temporary file cleanup."""
        # Create some temporary files
        temp_files = []
        for i in range(3):
            temp_file = self.manager.get_temp_file(suffix=f'.tmp{i}')
            temp_file.touch()
            temp_files.append(temp_file)
        
        # Verify they exist
        for temp_file in temp_files:
            assert temp_file.exists()
        
        # Cleanup
        self.manager.cleanup_temp_files()
        
        # Verify they're gone
        for temp_file in temp_files:
            assert not temp_file.exists()

    def test_safe_filename(self):
        """Test filename sanitization."""
        unsafe_name = "file with spaces & special/chars.fits"
        safe_name = self.manager.safe_filename(unsafe_name)
        
        assert ' ' not in safe_name
        assert '&' not in safe_name
        assert '/' not in safe_name
        assert safe_name.endswith('.fits')

    def test_get_file_info(self):
        """Test file information retrieval."""
        # Create a test file
        test_file = self.manager.base_dir / 'test.fits'
        test_file.write_text('test content')
        
        info = self.manager.get_file_info(test_file)
        
        assert info['exists'] is True
        assert info['size'] > 0
        assert 'modified_time' in info
        assert 'created_time' in info


class TestFormatConverter:
    """Test FormatConverter class functionality."""

    def setup_method(self):
        """Set up test fixtures."""
        self.converter = FormatConverter()
        self.sample_table = Table({
            'ra': [150.0, 151.0, 152.0],
            'dec': [2.0, 2.1, 2.2],
            'mag': [20.5, 21.0, 19.8]
        })

    def test_init(self):
        """Test initialization."""
        assert self.converter.supported_formats == ['csv', 'fits', 'votable', 'hdf5', 'parquet']

    def test_table_to_pandas(self):
        """Test astropy Table to pandas DataFrame conversion."""
        df = self.converter.table_to_pandas(self.sample_table)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3
        assert 'ra' in df.columns

    def test_pandas_to_table(self):
        """Test pandas DataFrame to astropy Table conversion."""
        df = self.sample_table.to_pandas()
        table = self.converter.pandas_to_table(df)
        assert isinstance(table, Table)
        assert len(table) == 3
        assert 'ra' in table.colnames

    def test_convert_format(self):
        """Test format conversion between file formats."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create input CSV file
            input_file = Path(temp_dir) / 'input.csv'
            self.sample_table.write(str(input_file), format='csv')
            
            # Convert to FITS
            output_file = Path(temp_dir) / 'output.fits'
            self.converter.convert_format(
                str(input_file), str(output_file),
                input_format='csv', output_format='fits'
            )
            
            # Verify conversion
            assert output_file.exists()
            converted_table = Table.read(str(output_file))
            assert len(converted_table) == 3

    def test_convert_format_auto_detect(self):
        """Test format conversion with automatic format detection."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create input file
            input_file = Path(temp_dir) / 'input.fits'
            self.sample_table.write(str(input_file), format='fits')
            
            # Convert to CSV (auto-detect input format)
            output_file = Path(temp_dir) / 'output.csv'
            self.converter.convert_format(
                str(input_file), str(output_file),
                output_format='csv'
            )
            
            # Verify conversion
            assert output_file.exists()
            converted_table = Table.read(str(output_file), format='csv')
            assert len(converted_table) == 3

    def test_batch_convert(self):
        """Test batch format conversion."""
        with tempfile.TemporaryDirectory() as temp_dir:
            input_dir = Path(temp_dir) / 'input'
            output_dir = Path(temp_dir) / 'output'
            input_dir.mkdir()
            output_dir.mkdir()
            
            # Create multiple input files
            input_files = []
            for i in range(3):
                input_file = input_dir / f'table_{i}.csv'
                self.sample_table.write(str(input_file), format='csv')
                input_files.append(str(input_file))
            
            # Batch convert
            result_files = self.converter.batch_convert(
                input_files, str(output_dir),
                input_format='csv', output_format='fits'
            )
            
            # Verify all files were converted
            assert len(result_files) == 3
            for result_file in result_files:
                assert Path(result_file).exists()
                assert Path(result_file).suffix == '.fits'

    def test_get_format_info(self):
        """Test format information retrieval."""
        csv_info = self.converter.get_format_info('csv')
        assert 'description' in csv_info
        assert 'extensions' in csv_info
        assert '.csv' in csv_info['extensions']

    def test_validate_format(self):
        """Test format validation."""
        assert self.converter.validate_format('csv') is True
        assert self.converter.validate_format('fits') is True
        assert self.converter.validate_format('invalid_format') is False

    def test_estimate_conversion_time(self):
        """Test conversion time estimation."""
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as f:
            self.sample_table.write(f.name, format='csv')
            
        try:
            time_estimate = self.converter.estimate_conversion_time(
                f.name, 'csv', 'fits'
            )
            assert isinstance(time_estimate, float)
            assert time_estimate > 0
        finally:
            os.unlink(f.name)
