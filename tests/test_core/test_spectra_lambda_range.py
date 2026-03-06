"""Tests for XML-based LambdaRange annotation/filtering in SpectrumCompiler."""

from pathlib import Path

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
