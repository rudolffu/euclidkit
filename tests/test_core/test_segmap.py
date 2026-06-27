"""Tests for local segmentation-map cutout compilation."""

from pathlib import Path
from typing import Optional

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from euclidkit.core.segmap import compile_segmap_cutouts


def _write_segmap(path: Path, data: Optional[np.ndarray] = None) -> None:
    if data is None:
        data = np.zeros((50, 50), dtype=np.int64)
        data[24:27, 24:27] = 42
        data[27:30, 27:30] = 99
    header = fits.Header()
    header["NAXIS"] = 2
    header["NAXIS1"] = data.shape[1]
    header["NAXIS2"] = data.shape[0]
    header["CTYPE1"] = "RA---TAN"
    header["CTYPE2"] = "DEC--TAN"
    header["CRPIX1"] = 25.0
    header["CRPIX2"] = 25.0
    header["CRVAL1"] = 150.0
    header["CRVAL2"] = 2.0
    header["CDELT1"] = -1.0 / 3600.0
    header["CDELT2"] = 1.0 / 3600.0
    fits.PrimaryHDU(data=data, header=header).writeto(path, overwrite=True)


def _center_world(path: Path, x: float = 25.0, y: float = 25.0) -> tuple[float, float]:
    with fits.open(path) as hdul:
        wcs = WCS(hdul[0].header)
        ra, dec = wcs.pixel_to_world_values(x, y)
    return float(ra), float(dec)


def _catalog(segmap_path: Path, *, use_mer_coords: bool = False) -> Table:
    ra, dec = _center_world(segmap_path)
    cols = {
        "datalabs_path": [str(segmap_path.parent)],
        "file_name": [segmap_path.name],
        "object_id": [123456789],
        "segmentation_map_id": [42000001],
        "tile_index": [101158583],
    }
    if use_mer_coords:
        cols["mer_ra"] = [ra]
        cols["mer_dec"] = [dec]
    else:
        cols["ra"] = [ra]
        cols["dec"] = [dec]
    return Table(cols)


def test_compile_segmap_cutouts_writes_raw_label_fits(tmp_path):
    """One source row should produce one raw-label FITS cutout with metadata."""
    segmap_path = tmp_path / "EUC_MER_FINAL-SEGMAP_TILE101158583-test.fits"
    _write_segmap(segmap_path)
    output_dir = tmp_path / "cutouts"

    stats = compile_segmap_cutouts(
        _catalog(segmap_path),
        output_dir,
        size_arcsec=10.0,
        show_progress=False,
    )

    assert stats.requested_rows == 1
    assert stats.written_rows == 1
    assert stats.failed_rows == 0
    assert stats.output_file_count == 1
    out_path = Path(stats.output_files[0])
    assert out_path.name == "segmap_cutout_row000000_obj123456789_seg42000001.fits"

    with fits.open(out_path) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        assert data.shape == (10, 10)
        assert 42 in data
        assert 99 in data
        assert header["OBJECTID"] == 123456789
        assert header["SEGMAPID"] == 42000001
        assert header["TILEIND"] == 101158583
        assert header["SRCFILE"] == segmap_path.name
        assert "CTYPE1" in header


def test_compile_segmap_cutouts_accepts_mer_ra_dec(tmp_path):
    """mer_ra/mer_dec should be accepted when ra/dec are absent."""
    segmap_path = tmp_path / "segmap.fits"
    _write_segmap(segmap_path)

    stats = compile_segmap_cutouts(
        _catalog(segmap_path, use_mer_coords=True),
        tmp_path / "cutouts",
        size_arcsec=10.0,
        show_progress=False,
    )

    assert stats.written_rows == 1
    with fits.open(stats.output_files[0]) as hdul:
        assert hdul[0].header["SRC_RA"] == pytest.approx(_catalog(segmap_path, use_mer_coords=True)["mer_ra"][0])


def test_compile_segmap_cutouts_skips_existing_unless_overwrite(tmp_path):
    """Existing output files are skipped by default and replaced with overwrite."""
    segmap_path = tmp_path / "segmap.fits"
    _write_segmap(segmap_path)
    output_dir = tmp_path / "cutouts"

    first = compile_segmap_cutouts(_catalog(segmap_path), output_dir, show_progress=False)
    second = compile_segmap_cutouts(_catalog(segmap_path), output_dir, show_progress=False)
    third = compile_segmap_cutouts(_catalog(segmap_path), output_dir, overwrite=True, show_progress=False)

    assert first.written_rows == 1
    assert second.written_rows == 0
    assert second.skipped_existing_rows == 1
    assert third.written_rows == 1


def test_compile_segmap_cutouts_missing_file_skip(tmp_path):
    """on_error=skip should record missing segmentation-map files as failures."""
    table = Table({
        "datalabs_path": [str(tmp_path)],
        "file_name": ["missing.fits"],
        "object_id": [1],
        "segmentation_map_id": [2],
        "ra": [150.0],
        "dec": [2.0],
    })

    stats = compile_segmap_cutouts(
        table,
        tmp_path / "cutouts",
        on_error="skip",
        show_progress=False,
    )

    assert stats.written_rows == 0
    assert stats.failed_rows == 1
    assert "not found" in stats.failures[0]["error"]


def test_compile_segmap_cutouts_missing_required_columns(tmp_path):
    """Missing required columns should raise a clear validation error."""
    table = Table({"datalabs_path": [str(tmp_path)], "file_name": ["segmap.fits"]})

    with pytest.raises(ValueError, match="segmentation_map_id"):
        compile_segmap_cutouts(table, tmp_path / "cutouts", show_progress=False)


def test_compile_segmap_cutouts_edge_source_uses_partial_mode(tmp_path):
    """Boundary sources should be padded rather than failing."""
    segmap_path = tmp_path / "segmap.fits"
    data = np.ones((50, 50), dtype=np.int64)
    _write_segmap(segmap_path, data=data)
    ra, dec = _center_world(segmap_path, x=1.0, y=1.0)
    table = _catalog(segmap_path)
    table["ra"] = [ra]
    table["dec"] = [dec]

    stats = compile_segmap_cutouts(table, tmp_path / "cutouts", size_arcsec=10.0, show_progress=False)

    assert stats.written_rows == 1
    with fits.open(stats.output_files[0]) as hdul:
        assert hdul[0].data.shape == (10, 10)
        assert 0 in hdul[0].data
