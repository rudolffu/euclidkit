"""Current spectrum HDU loading, validation and compilation contracts.

Removed APIs (loader caching, generic normalization/redshift transforms and
private HDU/chunk builders) are intentionally not recreated by these tests.
"""
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from euclidkit.core.spectra import SpectrumCompiler, SpectrumLoader, SpectrumProcessor


def write_spectra(path, count=3):
    hdus = [fits.PrimaryHDU()]
    for i in range(count):
        payload = np.array([(4000., float(i)), (5000., float(i + 1))],
                           dtype=[('WAVELENGTH', 'f8'), ('SIGNAL', 'f8')])
        hdus.append(fits.BinTableHDU(data=payload))
    fits.HDUList(hdus).writeto(path)
    return path


def catalog(path, count=3):
    return Table({'source_id': np.arange(count) + 100,
                  'datalabs_path': [str(path.parent)]*count,
                  'file_name': [path.name]*count,
                  'hdu_index': np.arange(count) + 1,
                  'ra_obj': np.arange(count) + 150., 'dec_obj': np.full(count, 2.)})


@pytest.mark.parametrize('astropy_only', [False, True])
def test_loader_returns_detached_hdu_payload(tmp_path, astropy_only):
    path = write_spectra(tmp_path/'spectra.fits')
    loader = SpectrumLoader()
    if astropy_only:
        with patch('euclidkit.core.spectra.fitsio', None):
            result = loader.load_spectrum(str(path), hdu_index=2)
    else:
        result = loader.load_spectrum(str(path), hdu_index=2)
    assert isinstance(result, fits.BinTableHDU)
    np.testing.assert_array_equal(result.data['WAVELENGTH'], [4000., 5000.])
    np.testing.assert_array_equal(result.data['SIGNAL'], [1., 2.])
    path.unlink()
    np.testing.assert_array_equal(result.data['SIGNAL'], [1., 2.])


def test_loader_batch_indices_and_errors(tmp_path):
    path = write_spectra(tmp_path/'spectra.fits')
    loader = SpectrumLoader()
    loaded = loader.load_spectra_batch(str(path), [1, 3, 1])
    assert sorted(loaded) == [1, 3]
    assert loader.load_spectra_batch(str(path), []) == {}
    with pytest.raises(ValueError, match='out of range'):
        loader.load_spectrum(str(path), 8)
    with pytest.raises((FileNotFoundError, OSError)):
        loader.load_spectrum(str(tmp_path/'absent.fits'))


def test_loader_verifies_local_path(tmp_path):
    loader = SpectrumLoader()
    assert loader.verify_spectrum_path(str(tmp_path), 'missing.fits') is None
    path = write_spectra(tmp_path/'spectra.fits')
    assert loader.verify_spectrum_path(str(tmp_path), path.name) == str(path)


@pytest.mark.parametrize('kind', ['valid', 'missing_columns', 'empty_hdu', 'image'])
def test_processor_validates_supported_hdus(kind):
    processor = SpectrumProcessor()
    if kind == 'valid':
        hdu = fits.BinTableHDU(Table({'WAVELENGTH': [4000., 5000.], 'SIGNAL': [-1., 3.]}))
    elif kind == 'missing_columns':
        hdu = fits.BinTableHDU(Table({'WAVELENGTH': [4000.]}))
    elif kind == 'empty_hdu':
        hdu = fits.BinTableHDU()
    else:
        hdu = fits.ImageHDU(np.ones((2, 2)))
    result = processor.validate_spectrum(hdu)
    assert result['valid'] == (kind == 'valid')
    if kind == 'valid':
        assert result['n_pixels'] == 2
        assert result['wavelength_range'] == (4000., 5000.)
        assert result['flux_range'] == (-1., 3.)
    else:
        assert result['issues']


def test_compiler_current_default_and_required_columns(tmp_path):
    compiler = SpectrumCompiler()
    assert compiler.max_extensions == 1000
    assert isinstance(compiler.loader, SpectrumLoader)
    with pytest.raises(ValueError, match='Missing required columns'):
        compiler.compile_spectra(Table({'source_id': [1]}), tmp_path)


def test_compiler_chunks_preserve_data_and_source_headers(tmp_path):
    path = write_spectra(tmp_path/'spectra.fits', count=5)
    source = catalog(path, count=5)
    compiler = SpectrumCompiler(max_extensions=2)
    outputs = compiler.compile_spectra(source, tmp_path/'out', output_prefix='compiled')
    assert [Path(p).name for p in outputs] == ['compiled_chunk_001.fits',
        'compiled_chunk_002.fits', 'compiled_chunk_003.fits']
    ids, fluxes = [], []
    for output, count in zip(outputs, [2, 2, 1]):
        with fits.open(output) as hdus:
            assert len(hdus)-1 == count
            assert hdus[0].header['NSPECTRA'] == count
            for hdu in hdus[1:]:
                ids.append(hdu.header['SOURC_ID'])
                fluxes.append(hdu.data['SIGNAL'][0])
                assert hdu.header['RA'] == source['ra_obj'][len(ids)-1]
    assert ids == [100, 101, 102, 103, 104]
    assert fluxes == [0., 1., 2., 3., 4.]
    assert 'chunk_number' not in source.colnames


def test_compiler_repeated_hdu_has_independent_source_headers(tmp_path):
    path = write_spectra(tmp_path/'spectra.fits', count=1)
    source = catalog(path, count=1)
    source = Table(np.repeat(source.as_array(), 2))
    source['source_id'] = [100, 200]
    source['ra_obj'] = [150., 151.]
    outputs = SpectrumCompiler().compile_spectra(source, tmp_path/'out')
    with fits.open(outputs[0]) as hdus:
        assert [h.header['SOURC_ID'] for h in hdus[1:]] == [100, 200]
        assert [h.header['RA'] for h in hdus[1:]] == [150., 151.]


def test_compiler_missing_file_counts_written_spectra(tmp_path):
    path = write_spectra(tmp_path/'spectra.fits')
    source = catalog(path)
    source['file_name'][1] = 'absent.fits'
    outputs = SpectrumCompiler().compile_spectra(source, tmp_path/'out')
    with fits.open(outputs[0]) as hdus:
        assert len(hdus) == 3
        assert hdus[0].header['NSPECTRA'] == len(hdus)-1
        assert [h.header['SOURC_ID'] for h in hdus[1:]] == [100, 102]


def test_compiler_resume_and_overwrite(tmp_path):
    path = write_spectra(tmp_path/'spectra.fits')
    source = catalog(path)
    compiler = SpectrumCompiler(max_extensions=2)
    out = tmp_path/'out'
    original = compiler.compile_spectra(source[:2], out)
    old_bytes = Path(original[0]).read_bytes()
    resumed = compiler.compile_spectra(source, out, overwrite=False)
    assert len(resumed) == 2
    assert Path(original[0]).read_bytes() == old_bytes
    compiler.compile_spectra(source[:1], out, overwrite=True)
    with fits.open(original[0]) as hdus:
        assert len(hdus) == 2


def test_compiler_metadata_tracks_real_chunk_sizes(tmp_path):
    path = write_spectra(tmp_path/'spectra.fits')
    source = catalog(path)
    original = source.copy()
    compiler = SpectrumCompiler(max_extensions=2)
    outputs = compiler.compile_spectra(source, tmp_path/'out')
    metadata_file = compiler.create_metadata_table(source, outputs, tmp_path/'out')
    metadata = Table.read(metadata_file)
    assert list(metadata['chunk_number']) == [1, 1, 2]
    assert list(metadata['extension_number']) == [1, 2, 1]
    assert source.colnames == original.colnames
    np.testing.assert_array_equal(source.as_array(), original.as_array())
