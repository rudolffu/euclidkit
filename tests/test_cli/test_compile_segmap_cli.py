"""Tests for the compile-segmap CLI command."""

from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from euclidkit.cli.crossmatch_cli import compile_segmap
from euclidkit.core.segmap import SegmapCutoutStats


def test_compile_segmap_cli_forwards_options_and_prints_summary(tmp_path):
    """compile-segmap should call the core compiler and print count summaries."""
    input_path = tmp_path / "segmaps.fits"
    input_path.write_bytes(b"placeholder")
    output_dir = tmp_path / "cutouts"
    stats = SegmapCutoutStats(
        requested_rows=3,
        written_rows=2,
        failed_rows=1,
        skipped_existing_rows=0,
        output_files=[str(output_dir / "a.fits"), str(output_dir / "b.fits")],
    )

    with patch("euclidkit.core.segmap.compile_segmap_cutouts", return_value=stats) as mock_compile:
        result = CliRunner().invoke(compile_segmap, [
            "--input", str(input_path),
            "--output-dir", str(output_dir),
            "--size-arcsec", "12.5",
            "--overwrite",
            "--on-error", "skip",
            "--no-progress",
        ])

    assert result.exit_code == 0
    mock_compile.assert_called_once_with(
        input_table=str(input_path),
        output_dir=str(output_dir),
        size_arcsec=12.5,
        overwrite=True,
        on_error="skip",
        show_progress=False,
    )
    assert "Segmentation-map cutouts completed." in result.output
    assert "Requested rows: 3" in result.output
    assert "Written cutouts: 2" in result.output
    assert "Failed rows: 1" in result.output
    assert "Output files: 2" in result.output
