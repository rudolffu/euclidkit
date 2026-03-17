"""Tests for XML-based LambdaRange annotation/filtering in SpectrumCompiler."""

from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table

from euclidkit.core.spectra import SpectrumCompiler


def _write_xml(path: Path, lambda_range: str, file_names: list[str]) -> None:
    files_xml = "\n".join([f"  <FileName>{name}</FileName>" for name in file_names])
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<Root>
  <Data>
    <LambdaRange>{lambda_range}</LambdaRange>
{files_xml}
  </Data>
</Root>
"""
    path.write_text(xml, encoding="utf-8")


def test_annotate_lambda_range_from_xml_matches_fits_gz_names(tmp_path):
    compiler = SpectrumCompiler(max_extensions=10)
    dpath = tmp_path / "dl"
    dpath.mkdir(parents=True, exist_ok=True)

    _write_xml(
        dpath / "rgs.xml",
        "RGS",
        ["EUC_SIR_W-COMBSPEC_A_2025T000000.fits.gz"],
    )
    _write_xml(
        dpath / "bgs.xml",
        "BGS",
        ["EUC_SIR_W-COMBSPEC_B_2025T000000.fits.gz"],
    )

    table = Table(
        {
            "source_id": [1, 2, 3],
            "datalabs_path": [str(dpath), str(dpath), str(dpath)],
            "file_name": [
                "EUC_SIR_W-COMBSPEC_A_2025T000000.fits",
                "EUC_SIR_W-COMBSPEC_B_2025T000000.fits",
                "EUC_SIR_W-COMBSPEC_C_2025T000000.fits",
            ],
            "hdu_index": [1, 1, 1],
        }
    )

    annotated, stats = compiler.annotate_lambda_range_from_xml(table)
    assert list(annotated["lambda_range"]) == ["RGS", "BGS", "UNKNOWN"]
    assert stats["total"] == 3
    assert stats["resolved"] == 2
    assert stats["rgs"] == 1
    assert stats["bgs"] == 1
    assert stats["unresolved"] == 1


def test_filter_table_by_lambda_range_both_keeps_resolved_only():
    compiler = SpectrumCompiler(max_extensions=10)
    table = Table(
        {
            "source_id": [1, 2, 3, 4],
            "lambda_range": ["RGS", "BGS", "UNKNOWN", "AMBIGUOUS"],
        }
    )

    selected, stats = compiler.filter_table_by_lambda_range(table, lambda_range="BOTH")
    assert len(selected) == 2
    assert set(selected["lambda_range"]) == {"RGS", "BGS"}
    assert stats["selected"] == 2
    assert stats["unresolved"] == 1
    assert stats["ambiguous"] == 1
    assert stats["skipped_other_band"] == 0


def test_deduplicate_by_source_id_keeps_first_row_order():
    compiler = SpectrumCompiler(max_extensions=10)
    table = Table(
        {
            "source_id": [11, 11, 12, 13, 12],
            "hdu_index": [1, 2, 3, 4, 5],
        }
    )
    dedup, stats = compiler.deduplicate_by_source_id(table)
    assert list(dedup["source_id"]) == [11, 12, 13]
    assert list(dedup["hdu_index"]) == [1, 3, 4]
    assert stats["input_rows"] == 5
    assert stats["unique_sources"] == 3
    assert stats["duplicate_rows_removed"] == 2


def test_deduplicate_by_source_id_falls_back_to_object_id():
    compiler = SpectrumCompiler(max_extensions=10)
    table = Table(
        {
            "object_id": [11, 11, 12, 13, 12],
            "hdu_index": [1, 2, 3, 4, 5],
        }
    )
    dedup, stats = compiler.deduplicate_by_source_id(table)
    assert list(dedup["object_id"]) == [11, 12, 13]
    assert list(dedup["hdu_index"]) == [1, 3, 4]
    assert stats["input_rows"] == 5
    assert stats["unique_sources"] == 3
    assert stats["duplicate_rows_removed"] == 2


def test_compile_spectra_datalink_all_fetches_both_arms(tmp_path):
    compiler = SpectrumCompiler(max_extensions=1000)
    table = Table(
        {
            "source_id": [101, 101, 102, 102],
            "object_id": [101, 101, 102, 102],
            "right_ascension": [10.1, 10.1, 10.2, 10.2],
            "declination": [-1.1, -1.1, -1.2, -1.2],
        }
    )

    # Minimal bintable HDU payload returned by loader for each call.
    payload = np.array([(1.0, 2.0)], dtype=[("WAVELENGTH", "f4"), ("SIGNAL", "f4")])
    hdu = fits.BinTableHDU(data=payload)

    calls = []

    def fake_load(*, euclid_client, source_id, schema, retrieval_type):
        calls.append((str(source_id), retrieval_type))
        return hdu.copy()

    compiler.loader.load_spectrum_from_datalink = fake_load  # type: ignore[assignment]
    out_files = compiler.compile_spectra_datalink(
        spectra_table=table,
        euclid_client=object(),
        output_dir=tmp_path,
        output_prefix="dl_all",
        retrieval_type="ALL",
        schema="sedm",
        overwrite=True,
    )

    assert len(out_files) == 1
    # 2 unique sources x 2 arms = 4 retrievals
    assert len(calls) == 4
    assert set(calls) == {
        ("101", "SPECTRA_RGS"),
        ("101", "SPECTRA_BGS"),
        ("102", "SPECTRA_RGS"),
        ("102", "SPECTRA_BGS"),
    }

    with fits.open(out_files[0]) as hdul:
        # Primary + 4 spectra extensions.
        assert len(hdul) == 5
        lambdas = [h.header.get("LRANGE") for h in hdul[1:]]
        assert sorted(lambdas) == ["BGS", "BGS", "RGS", "RGS"]
        # RA/DEC should be written from right_ascension/declination aliases.
        assert hdul[1].header.get("RA") is not None
        assert hdul[1].header.get("DEC") is not None


def test_compile_spectra_datalink_uses_object_id_when_source_id_missing(tmp_path):
    compiler = SpectrumCompiler(max_extensions=1000)
    table = Table(
        {
            "object_id": [101, 101, 102, 102],
            "right_ascension": [10.1, 10.1, 10.2, 10.2],
            "declination": [-1.1, -1.1, -1.2, -1.2],
        }
    )

    payload = np.array([(1.0, 2.0)], dtype=[("WAVELENGTH", "f4"), ("SIGNAL", "f4")])
    hdu = fits.BinTableHDU(data=payload)

    calls = []

    def fake_load(*, euclid_client, source_id, schema, retrieval_type):
        calls.append((str(source_id), retrieval_type))
        return hdu.copy()

    compiler.loader.load_spectrum_from_datalink = fake_load  # type: ignore[assignment]
    out_files = compiler.compile_spectra_datalink(
        spectra_table=table,
        euclid_client=object(),
        output_dir=tmp_path,
        output_prefix="dl_oid",
        retrieval_type="SPECTRA_RGS",
        schema="sedm",
        overwrite=True,
    )

    assert len(out_files) == 1
    assert len(calls) == 2
    assert set(calls) == {
        ("101", "SPECTRA_RGS"),
        ("102", "SPECTRA_RGS"),
    }

    with fits.open(out_files[0]) as hdul:
        assert len(hdul) == 3
        assert hdul[1].header["SOURC_ID"] == "101"
        assert hdul[2].header["SOURC_ID"] == "102"


def test_resolve_radec_columns_aliases():
    compiler = SpectrumCompiler(max_extensions=10)
    table = Table({"right_ascension": [1.0], "declination": [2.0]})
    ra_col, dec_col = compiler._resolve_radec_columns(table)
    assert ra_col == "right_ascension"
    assert dec_col == "declination"
