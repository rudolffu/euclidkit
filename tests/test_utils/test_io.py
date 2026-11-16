"""Test I/O utility functionality."""

import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io.registry import IORegistryError
from astropy.table import Table

from euclidkit.utils.io import (
    DataLoader,
    FileManager,
    FormatConverter,
    load_table,
    save_table,
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
            f.write(b"not a valid table format")
            f.flush()
            
        try:
            with pytest.raises(IORegistryError):
                load_table(f.name)
        finally:
            os.unlink(f.name)

    def test_load_table_nonexistent_file(self):
        """Test error handling for non-existent files."""
        with pytest.raises(FileNotFoundError):
            load_table('/nonexistent/file.csv')



class TestDataLoader:
    """Only directory management is currently implemented by DataLoader."""

    def test_explicit_cache_and_clear(self, tmp_path):
        cache = tmp_path/'cache'
        loader = DataLoader(cache_dir=cache)
        assert loader.cache_dir == cache and cache.is_dir()
        (cache/'nested').mkdir()
        (cache/'nested'/'data').write_bytes(b'cached')
        loader.clear_cache()
        assert cache.is_dir() and not list(cache.iterdir())
        loader.clear_cache()

    def test_default_cache_is_isolated(self, tmp_path, monkeypatch):
        monkeypatch.setattr(Path, 'home', lambda: tmp_path)
        loader = DataLoader()
        assert loader.cache_dir == tmp_path/'.euclidkit'/'cache'
        assert loader.cache_dir.is_dir()


class TestFileManager:
    def test_file_info_and_path_validation(self, tmp_path):
        path = tmp_path/'data.csv'
        path.write_text('a\n1\n')
        missing = tmp_path/'missing.csv'
        info = FileManager.get_file_info(path)
        assert info['exists'] and info['is_file'] and not info['is_dir']
        assert info['size_bytes'] == path.stat().st_size
        assert info['suffix'] == '.csv' and info['name'] == 'data.csv'
        assert FileManager.get_file_info(missing) == {'exists': False}
        result = FileManager.validate_paths([path, missing])
        assert result['valid'] == [str(path)] and result['invalid'] == [str(missing)]
        assert result['total_size_mb'] == info['size_mb']

    def test_fits_info_reports_hdu_count(self, tmp_path):
        path = tmp_path/'data.fits'
        Table({'a': [1]}).write(path)
        info = FileManager.get_file_info(path)
        assert info['n_hdus'] == 2 and info['primary_header_cards'] > 0


class TestFormatConverter:
    def test_csv_fits_and_dataframe_roundtrip(self, tmp_path):
        original = Table({'ra': [150., 151.], 'object_id': [1001, 1002]})
        csv_path, fits_path = tmp_path/'sources.csv', tmp_path/'sources.fits'
        FormatConverter.table_to_csv(original, csv_path)
        FormatConverter.csv_to_fits(csv_path, fits_path)
        result = Table.read(fits_path)
        np.testing.assert_array_equal(result['ra'], original['ra'])
        np.testing.assert_array_equal(result['object_id'], original['object_id'])
        frame = FormatConverter.table_to_pandas(result)
        assert isinstance(frame, pd.DataFrame)
        converted = FormatConverter.pandas_to_table(frame)
        np.testing.assert_array_equal(converted['object_id'], original['object_id'])
