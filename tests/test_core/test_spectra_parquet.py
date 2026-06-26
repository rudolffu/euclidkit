from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.table import Table

from euclidkit.core.spectra_parquet import (
    dithers_to_parquet,
    parse_combined_extname,
    parse_dither_extname,
    spectra_to_parquet,
)


def _make_signal_hdu(
    object_id: int,
    wave: np.ndarray,
    signal: np.ndarray,
    var: np.ndarray,
    mask: np.ndarray,
    *,
    include_wavelength: bool = True,
    include_quality: bool = True,
) -> fits.BinTableHDU:
    cols = []
    if include_wavelength:
        cols.append(fits.Column(name="WAVELENGTH", array=np.asarray([wave]), format=f"{wave.size}D"))
    cols.extend(
        [
            fits.Column(name="SIGNAL", array=np.asarray([signal], dtype=np.float32), format=f"{wave.size}E"),
            fits.Column(name="VAR", array=np.asarray([var], dtype=np.float32), format=f"{wave.size}E"),
            fits.Column(name="MASK", array=np.asarray([mask], dtype=np.int32), format=f"{wave.size}J"),
            fits.Column(name="ERR", array=np.asarray([var], dtype=np.float32), format=f"{wave.size}E"),
            fits.Column(name="VALID", array=np.asarray([[True] * wave.size]), format=f"{wave.size}L"),
            fits.Column(name="IVAR", array=np.asarray([var], dtype=np.float32), format=f"{wave.size}E"),
        ]
    )
    if include_quality:
        cols.append(fits.Column(name="QUALITY", array=np.asarray([np.arange(wave.size)], dtype=np.float32), format=f"{wave.size}E"))
    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.header["EXTNAME"] = f"{object_id}_COMBINED1D_SIGNAL"
    hdu.header["OBJ_ID"] = int(object_id)
    hdu.header["WMIN"] = float(wave[0])
    hdu.header["BINWIDTH"] = float(wave[1] - wave[0])
    hdu.header["BINCOUNT"] = int(wave.size)
    return hdu


def _write_raw_test_fits(path: Path, *, lrange: str = "RGS") -> tuple[int, int]:
    wave = np.array([100.0, 101.5, 103.0, 104.5], dtype=np.float64)
    hdus = [
        fits.PrimaryHDU(),
        _make_signal_hdu(
            9001,
            wave,
            np.array([1.0, 1.5, 2.0, 2.5], dtype=np.float32),
            np.array([0.1, 0.2, 0.3, 0.4], dtype=np.float32),
            np.array([0, 1, 0, 0], dtype=np.int32),
            include_wavelength=True,
            include_quality=True,
        ),
        _make_signal_hdu(
            9002,
            wave,
            np.array([3.0, 3.5, 4.0, 4.5], dtype=np.float32),
            np.array([1.1, 1.2, 1.3, 1.4], dtype=np.float32),
            np.array([1, 0, 0, 1], dtype=np.int32),
            include_wavelength=False,
            include_quality=False,
        ),
    ]
    hdus[0].header["LRANGE"] = lrange
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return 1, 2


def _write_catalog(path: Path, datalabs_path: Path, file_name: str, hdu_a: int, hdu_b: int) -> None:
    df = pd.DataFrame(
        {
            "datalabs_path": [str(datalabs_path), str(datalabs_path)],
            "file_name": [file_name, file_name],
            "hdu_index": [hdu_b, hdu_a],
            "source_id": [9002, 9001],
            "object_id": [9002, 9001],
            "ra_obj": [12.5, 11.5],
            "dec_obj": [-0.5, 0.5],
        }
    )
    if path.suffix == ".fits":
        Table.from_pandas(df).write(path, format="fits", overwrite=True)
    else:
        df.to_parquet(path, index=False)


def _make_dither_hdu(prefix: int, dither_id: int, wave: np.ndarray, signal: np.ndarray) -> fits.BinTableHDU:
    hdu = _make_signal_hdu(
        prefix,
        wave,
        signal,
        np.ones_like(signal, dtype=np.float32),
        np.zeros_like(signal, dtype=np.int32),
        include_wavelength=False,
    )
    hdu.header["EXTNAME"] = f"{prefix}_DITH1D_{dither_id}_SIGNAL"
    hdu.header["DITH_ID"] = int(dither_id)
    hdu.header["PTGID"] = int(dither_id)
    hdu.header["GWA_POS"] = "RGS000"
    hdu.header["GWA_TILT"] = 0
    hdu.header["EXPTIME"] = 565.0
    hdu.header["DET_ID"] = 16
    hdu.header["EXT_PROF"] = "OPT_PSF"
    return hdu


def _write_dither_test_fits(path: Path, *, lrange: str = "RGS") -> int:
    wave = np.array([100.0, 101.5, 103.0, 104.5], dtype=np.float64)
    combined = _make_signal_hdu(
        995,
        wave,
        np.array([1.0, 1.5, 2.0, 2.5], dtype=np.float32),
        np.array([0.1, 0.2, 0.3, 0.4], dtype=np.float32),
        np.array([0, 1, 0, 0], dtype=np.int32),
    )
    combined.header["EXTNAME"] = "995_COMBINED1D_SIGNAL"
    hdus = [
        fits.PrimaryHDU(),
        combined,
        _make_dither_hdu(995, 3197, wave, np.array([3.0, 3.5, 4.0, 4.5], dtype=np.float32)),
        _make_dither_hdu(995, 3198, wave, np.array([5.0, 5.5, 6.0, 6.5], dtype=np.float32)),
        _make_dither_hdu(997, 3318, wave, np.array([7.0, 7.5, 8.0, 8.5], dtype=np.float32)),
    ]
    hdus[0].header["LRANGE"] = lrange
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return 1


def test_spectra_to_parquet_exports_arrays_chunks_and_manifest(tmp_path: Path):
    fits_dir = tmp_path / "fits"
    fits_dir.mkdir()
    file_name = "EUC_SIR_W-COMBSPEC_101_20250101T000000.000000Z.fits"
    hdu_a, hdu_b = _write_raw_test_fits(fits_dir / file_name)
    catalog = tmp_path / "catalog.fits"
    _write_catalog(catalog, fits_dir, file_name, hdu_a, hdu_b)

    stats = spectra_to_parquet(
        str(catalog),
        str(tmp_path / "raw_spectra"),
        chunk_size=1,
        workers=1,
        lambda_range="RGS",
    )

    assert stats.exported_rows == 2
    part1 = pd.read_parquet(tmp_path / "raw_spectra_part001.parquet")
    part2 = pd.read_parquet(tmp_path / "raw_spectra_part002.parquet")
    assert part1["object_id"].tolist() == [9002]
    assert part1["source_id"].tolist() == [9002]
    assert part1["ra"].tolist() == [12.5]
    assert part1["dec"].tolist() == [-0.5]
    for forbidden in ("file_name", "hdu_index", "err", "valid", "ivar"):
        assert forbidden not in part1.columns
    for required in ("wavelength", "signal", "var", "mask", "flux"):
        assert required in part1.columns
    np.testing.assert_allclose(np.asarray(part1["flux"].iloc[0]), [3.0, 3.5, 4.0, 4.5])
    np.testing.assert_allclose(np.asarray(part1["wavelength"].iloc[0]), [100.0, 101.5, 103.0, 104.5])
    assert "quality" not in part1.columns
    assert "quality" in part2.columns

    manifest = json.loads((tmp_path / "raw_spectra_manifest.json").read_text(encoding="utf-8"))
    assert manifest["requested_rows"] == 2
    assert manifest["exported_rows"] == 2
    assert manifest["chunk_size"] == 1
    assert manifest["workers"] == 1
    assert manifest["lambda_range"] == "RGS"


def test_spectra_to_parquet_limit_lambda_filter_overwrite_and_failures(tmp_path: Path):
    fits_dir = tmp_path / "fits"
    fits_dir.mkdir()
    file_name = "EUC_SIR_W-COMBSPEC_101_20250101T000000.000000Z.fits"
    hdu_a, hdu_b = _write_raw_test_fits(fits_dir / file_name, lrange="BGS")
    catalog = tmp_path / "catalog.parquet"
    _write_catalog(catalog, fits_dir, file_name, hdu_a, hdu_b)

    stats = spectra_to_parquet(str(catalog), str(tmp_path / "limited"), workers=1, limit=1)
    assert stats.requested_rows == 1
    assert stats.exported_rows == 1

    stats = spectra_to_parquet(str(catalog), str(tmp_path / "filtered"), workers=1, lambda_range="RGS")
    assert stats.exported_rows == 0
    assert stats.lambda_skipped_rows == 2
    assert not list(tmp_path.glob("filtered_part*.parquet"))

    with pytest.raises(FileExistsError):
        spectra_to_parquet(str(catalog), str(tmp_path / "limited"), workers=1)

    bad_catalog = tmp_path / "bad.parquet"
    pd.DataFrame(
        {
            "datalabs_path": [str(fits_dir), str(fits_dir)],
            "file_name": [file_name, file_name],
            "hdu_index": [hdu_a, 99],
            "source_id": [9001, 9999],
            "object_id": [9001, 9999],
            "ra_obj": [11.5, 99.0],
            "dec_obj": [0.5, 99.0],
        }
    ).to_parquet(bad_catalog, index=False)
    stats = spectra_to_parquet(str(bad_catalog), str(tmp_path / "skip_bad"), workers=1, on_error="skip")
    assert stats.exported_rows == 1
    assert stats.failed_rows == 1
    assert Path(stats.failures_path).exists()


def test_dithers_to_parquet_outputs_and_dithers_only(tmp_path: Path):
    fits_dir = tmp_path / "fits"
    fits_dir.mkdir()
    file_name = "EUC_SIR_W-COMBSPEC_101_20250101T000000.000000Z.fits"
    combined_hdu = _write_dither_test_fits(fits_dir / file_name)
    catalog = tmp_path / "catalog.parquet"
    pd.DataFrame(
        {
            "datalabs_path": [str(fits_dir)],
            "file_name": [file_name],
            "hdu_index": [combined_hdu],
            "source_id": [9001],
            "object_id": [9001],
            "ra_obj": [11.5],
            "dec_obj": [0.5],
        }
    ).to_parquet(catalog, index=False)

    assert parse_combined_extname("123_COMBINED1D_SIGNAL") == "123"
    assert parse_dither_extname("123_DITH1D_3197_SIGNAL") == ("123", 3197)

    stats = dithers_to_parquet(str(catalog), str(tmp_path / "raw_sir"), workers=1, lambda_range="RGS")
    assert stats.objects_exported == 1
    assert stats.combined_rows == 1
    assert stats.dither_rows == 2
    combined = pd.read_parquet(tmp_path / "raw_sir_combined_part001.parquet")
    dithers = pd.read_parquet(tmp_path / "raw_sir_dithers_part001.parquet")
    assert combined["object_id"].tolist() == [9001]
    assert combined["combined_hdu_index"].tolist() == [1]
    assert dithers["dither_id"].tolist() == [3197, 3198]
    np.testing.assert_allclose(np.asarray(dithers["wavelength"].iloc[0]), [100.0, 101.5, 103.0, 104.5])

    stats = dithers_to_parquet(
        str(catalog),
        str(tmp_path / "raw_sir_dithers_only"),
        workers=1,
        include_combined=False,
    )
    assert stats.combined_rows == 0
    assert stats.dither_rows == 2
    assert not list(tmp_path.glob("raw_sir_dithers_only_combined_part*.parquet"))
    assert (tmp_path / "raw_sir_dithers_only_dithers_part001.parquet").exists()
