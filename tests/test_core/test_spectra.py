"""Test spectral processing functionality."""

import pytest
import numpy as np
import tempfile
import os
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

from astropy.table import Table
from astropy.io import fits

from euclidkit.core.spectra import SpectrumLoader, SpectrumProcessor, SpectrumCompiler


class TestSpectrumLoader:
    """Test cases for SpectrumLoader class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.loader = SpectrumLoader()

    def test_init(self):
        """Test initialization."""
        assert self.loader.cache_dir is None
        assert self.loader.use_cache is True

    def test_init_with_cache_dir(self):
        """Test initialization with cache directory."""
        loader = SpectrumLoader(cache_dir='/tmp/cache', use_cache=False)
        assert loader.cache_dir == Path('/tmp/cache')
        assert loader.use_cache is False

    @patch('euclidkit.core.spectra.fits.open')
    def test_load_spectrum_from_fits(self, mock_fits_open):
        """Test loading spectrum from FITS file."""
        # Mock FITS file structure
        mock_hdu = Mock()
        mock_hdu.data = np.array([
            (5000.0, 1.0, 0.1),
            (5001.0, 1.1, 0.1),
            (5002.0, 0.9, 0.1)
        ], dtype=[('wavelength', 'f8'), ('flux', 'f8'), ('error', 'f8')])
        mock_hdu.header = {'OBJECT': 'TEST_QSO', 'INSTRUME': 'NISP'}
        
        mock_hdul = Mock()
        mock_hdul.__enter__.return_value = [mock_hdu]
        mock_fits_open.return_value = mock_hdul
        
        spectrum = self.loader.load_spectrum('/path/to/spectrum.fits')
        
        assert 'wavelength' in spectrum.colnames
        assert 'flux' in spectrum.colnames
        assert 'error' in spectrum.colnames
        assert len(spectrum) == 3


class TestSpectrumProcessor:
    """Test cases for SpectrumProcessor class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.processor = SpectrumProcessor()
        
        # Sample spectrum data
        self.sample_spectrum = Table({
            'wavelength': np.linspace(4000, 8000, 1000),
            'flux': np.random.normal(1.0, 0.1, 1000),
            'error': np.full(1000, 0.1)
        })

    def test_init(self):
        """Test initialization."""
        assert self.processor.redshift_tolerance == 0.1

    def test_apply_redshift_correction(self):
        """Test redshift correction."""
        redshift = 2.0
        corrected = self.processor.apply_redshift_correction(
            self.sample_spectrum, redshift
        )
        
        # Wavelengths should be divided by (1 + z)
        expected_wave = self.sample_spectrum['wavelength'] / (1 + redshift)
        np.testing.assert_array_almost_equal(
            corrected['wavelength'], expected_wave
        )

    def test_normalize_spectrum(self):
        """Test spectrum normalization."""
        normalized = self.processor.normalize_spectrum(self.sample_spectrum)
        
        # Check that flux is normalized
        assert 'flux_normalized' in normalized.colnames
        # Mean should be close to 1 for the normalization window
        wave_mask = (normalized['wavelength'] > 5000) & (normalized['wavelength'] < 6000)
        mean_flux = np.mean(normalized['flux_normalized'][wave_mask])
        assert abs(mean_flux - 1.0) < 0.1

    def test_resample_spectrum(self):
        """Test spectrum resampling."""
        new_wavelength = np.linspace(4500, 7500, 500)
        resampled = self.processor.resample_spectrum(
            self.sample_spectrum, new_wavelength
        )
        
        assert len(resampled) == 500
        np.testing.assert_array_equal(resampled['wavelength'], new_wavelength)

    def test_mask_bad_pixels(self):
        """Test bad pixel masking."""
        # Introduce some bad pixels
        spectrum = self.sample_spectrum.copy()
        spectrum['flux'][100:110] = -999  # Bad flux values
        spectrum['error'][200:205] = 0    # Zero errors
        
        masked = self.processor.mask_bad_pixels(spectrum)
        
        assert 'mask' in masked.colnames
        # Check that bad pixels are masked
        assert masked['mask'][105] is False  # Should be masked
        assert masked['mask'][202] is False  # Should be masked
        assert masked['mask'][50] is True    # Should be good


class TestSpectrumCompiler:
    """Test cases for SpectrumCompiler class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.compiler = SpectrumCompiler(max_extensions=10)
        
        # Sample spectra table
        self.sample_spectra_table = Table({
            'object_id': [100001, 100001, 100002, 100003],
            'spectrum_id': ['spec_1', 'spec_2', 'spec_3', 'spec_4'],
            'instrument_name': ['NISP', 'VIS', 'NISP', 'VIS'],
            'file_path': ['/data/spec1.fits', '/data/spec2.fits', 
                         '/data/spec3.fits', '/data/spec4.fits'],
            'ra': [150.0, 150.0, 151.0, 152.0],
            'dec': [2.0, 2.0, 2.1, 2.2]
        })

    def test_init(self):
        """Test initialization."""
        assert self.compiler.max_extensions == 10
        assert self.compiler.loader is not None

    def test_init_default_max_extensions(self):
        """Test default max extensions."""
        compiler = SpectrumCompiler()
        assert compiler.max_extensions == 5000

    def test_calculate_chunks(self):
        """Test chunk calculation."""
        chunks = self.compiler._calculate_chunks(25, max_per_chunk=10)
        assert len(chunks) == 3
        assert chunks == [(0, 10), (10, 20), (20, 25)]

    @patch('euclidkit.core.spectra.SpectrumLoader.load_spectrum')
    @patch('astropy.io.fits.HDUList.writeto')
    def test_compile_spectra_single_file(self, mock_writeto, mock_load_spectrum):
        """Test compiling spectra into a single file."""
        # Mock loaded spectra
        mock_spectrum = Table({
            'wavelength': np.linspace(4000, 8000, 100),
            'flux': np.random.normal(1.0, 0.1, 100),
            'error': np.full(100, 0.1)
        })
        mock_load_spectrum.return_value = mock_spectrum
        
        with tempfile.TemporaryDirectory() as temp_dir:
            output_files = self.compiler.compile_spectra(
                spectra_table=self.sample_spectra_table[:5],  # 5 spectra, within limit
                output_dir=temp_dir,
                output_prefix='test_compiled'
            )
            
            assert len(output_files) == 1
            assert 'test_compiled_001.fits' in output_files[0]
            mock_writeto.assert_called_once()

    @patch('euclidkit.core.spectra.SpectrumLoader.load_spectrum')
    @patch('astropy.io.fits.HDUList.writeto')
    def test_compile_spectra_multiple_files(self, mock_writeto, mock_load_spectrum):
        """Test compiling spectra into multiple files."""
        # Create more spectra than max_extensions
        large_spectra_table = Table({
            'spectrum_id': [f'spec_{i}' for i in range(25)],
            'file_path': [f'/data/spec{i}.fits' for i in range(25)],
            'object_id': [100000 + i for i in range(25)],
            'ra': np.full(25, 150.0),
            'dec': np.full(25, 2.0)
        })
        
        mock_spectrum = Table({
            'wavelength': np.linspace(4000, 8000, 100),
            'flux': np.random.normal(1.0, 0.1, 100),
            'error': np.full(100, 0.1)
        })
        mock_load_spectrum.return_value = mock_spectrum
        
        with tempfile.TemporaryDirectory() as temp_dir:
            output_files = self.compiler.compile_spectra(
                spectra_table=large_spectra_table,
                output_dir=temp_dir,
                output_prefix='test_compiled'
            )
            
            # Should create 3 files (10, 10, 5 spectra)
            assert len(output_files) == 3
            assert mock_writeto.call_count == 3

    @patch('euclidkit.core.spectra.SpectrumLoader.load_spectrum')
    def test_compile_spectra_with_overwrite(self, mock_load_spectrum):
        """Test compiling with overwrite option."""
        mock_spectrum = Table({
            'wavelength': np.linspace(4000, 8000, 100),
            'flux': np.random.normal(1.0, 0.1, 100),
            'error': np.full(100, 0.1)
        })
        mock_load_spectrum.return_value = mock_spectrum
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create existing file
            existing_file = Path(temp_dir) / 'test_compiled_001.fits'
            existing_file.touch()
            
            with patch('astropy.io.fits.HDUList.writeto') as mock_writeto:
                output_files = self.compiler.compile_spectra(
                    spectra_table=self.sample_spectra_table[:3],
                    output_dir=temp_dir,
                    output_prefix='test_compiled',
                    overwrite=True
                )
                
                # Should call writeto with overwrite=True
                mock_writeto.assert_called_once()
                assert mock_writeto.call_args[1]['overwrite'] is True

    @patch('euclidkit.core.spectra.SpectrumLoader.load_spectrum')
    def test_compile_spectra_skip_failed_loads(self, mock_load_spectrum):
        """Test that failed spectrum loads are skipped."""
        # Mock some successful and some failed loads
        def side_effect(file_path):
            if 'spec_2' in file_path:
                raise IOError("File not found")
            return Table({
                'wavelength': np.linspace(4000, 8000, 100),
                'flux': np.random.normal(1.0, 0.1, 100),
                'error': np.full(100, 0.1)
            })
        
        mock_load_spectrum.side_effect = side_effect
        
        with tempfile.TemporaryDirectory() as temp_dir:
            with patch('astropy.io.fits.HDUList.writeto') as mock_writeto:
                output_files = self.compiler.compile_spectra(
                    spectra_table=self.sample_spectra_table,
                    output_dir=temp_dir,
                    output_prefix='test_compiled'
                )
                
                # Should still create output files despite one failed load
                assert len(output_files) == 1
                mock_writeto.assert_called_once()

    def test_create_metadata_table(self):
        """Test metadata table creation."""
        output_files = ['/path/compiled_001.fits', '/path/compiled_002.fits']
        
        with tempfile.TemporaryDirectory() as temp_dir:
            with patch('astropy.table.Table.write') as mock_write:
                metadata_file = self.compiler.create_metadata_table(
                    spectra_table=self.sample_spectra_table,
                    output_files=output_files,
                    output_dir=temp_dir,
                    output_name='metadata.fits'
                )
                
                mock_write.assert_called_once()
                assert 'metadata.fits' in metadata_file

    def test_create_spectrum_hdu(self):
        """Test HDU creation from spectrum data."""
        spectrum = Table({
            'wavelength': np.array([4000, 4001, 4002]),
            'flux': np.array([1.0, 1.1, 0.9]),
            'error': np.array([0.1, 0.1, 0.1])
        })
        
        spectrum_info = {
            'spectrum_id': 'test_spec',
            'object_id': 100001,
            'ra': 150.0,
            'dec': 2.0,
            'instrument_name': 'NISP'
        }
        
        hdu = self.compiler._create_spectrum_hdu(spectrum, spectrum_info, extension_number=1)
        
        assert isinstance(hdu, fits.BinTableHDU)
        assert hdu.header['EXTNAME'] == 'SPECTRUM'
        assert hdu.header['EXTVER'] == 1
        assert hdu.header['SPECID'] == 'test_spec'
        assert hdu.header['OBJID'] == 100001
        assert hdu.header['RA'] == 150.0
        assert hdu.header['DEC'] == 2.0
        assert hdu.header['INSTRUME'] == 'NISP'

    def test_create_primary_hdu(self):
        """Test primary HDU creation."""
        file_info = {
            'output_file': 'compiled_001.fits',
            'n_extensions': 10,
            'spectra_range': (0, 10)
        }
        
        primary_hdu = self.compiler._create_primary_hdu(file_info)
        
        assert isinstance(primary_hdu, fits.PrimaryHDU)
        assert primary_hdu.header['FILENAME'] == 'compiled_001.fits'
        assert primary_hdu.header['NEXTEND'] == 10
        assert primary_hdu.header['SPECMIN'] == 0
        assert primary_hdu.header['SPECMAX'] == 10
        assert 'CREATOR' in primary_hdu.header
        assert 'DATE' in primary_hdu.header