"""Tests for compile-spectra datalink mode."""

import os
import tempfile
from unittest.mock import Mock, patch

from astropy.table import Table
from click.testing import CliRunner

from euclidkit.cli.crossmatch_cli import compile_spectra


def _make_temp_spectra_file(table: Table) -> str:
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.fits', delete=False)
    table.write(tmp.name, format='fits', overwrite=True)
    return tmp.name


def test_compile_spectra_use_datalink_mode():
    """compile-spectra should use datalink path and archive login when requested."""
    runner = CliRunner()
    spectra_table = Table({
        'source_id': [100001, 100002],
        'object_id': [100001, 100002],
    })
    spectra_file = _make_temp_spectra_file(spectra_table)

    try:
        with tempfile.TemporaryDirectory() as output_dir:
            mock_compiler = Mock()
            mock_compiler.deduplicate_by_source_id.return_value = (
                spectra_table,
                {'input_rows': 2, 'unique_sources': 2, 'duplicate_rows_removed': 0},
            )
            mock_compiler.compile_spectra_datalink.return_value = [f"{output_dir}/dl_compiled_chunk_001.fits"]
            mock_compiler.create_metadata_table.return_value = f"{output_dir}/dl_compiled_metadata.fits"

            mock_archive = Mock()
            mock_archive.euclid = Mock()

            with patch('euclidkit.core.spectra.SpectrumCompiler', return_value=mock_compiler):
                with patch('euclidkit.core.data_access.EuclidArchive', return_value=mock_archive):
                    with patch('euclidkit.utils.io.load_table', return_value=spectra_table):
                        result = runner.invoke(compile_spectra, [
                            '--spectra-table', spectra_file,
                            '--output-dir', output_dir,
                        '--output-format', 'fits',
                            '--prefix', 'dl_compiled',
                            '--use-datalink',
                            '--environment', 'IDR',
                            '--retrieval-type', 'SPECTRA_BGS',
                            '--schema', 'dr1',
                            '-L', 'BGS',
                        ])

            assert result.exit_code == 0
            mock_archive.login.assert_called_once()
            mock_archive.logout.assert_called_once()
            mock_compiler.compile_spectra_datalink.assert_called_once()
            kwargs = mock_compiler.compile_spectra_datalink.call_args.kwargs
            assert kwargs['retrieval_type'] == 'SPECTRA_BGS'
            assert kwargs['schema'] == 'dr1'
            assert kwargs['euclid_client'] is mock_archive.euclid
            assert "Datalink arm selection:" in result.output
    finally:
        if os.path.exists(spectra_file):
            os.unlink(spectra_file)


def test_compile_spectra_non_empty_dir_uses_resume_mode_by_default():
    """If output dir is non-empty, CLI should keep path and run in resume mode by default."""
    runner = CliRunner()
    spectra_table = Table({
        'source_id': [100001, 100002],
        'datalabs_path': ['/tmp', '/tmp'],
        'file_name': ['a.fits', 'b.fits'],
        'hdu_index': [1, 1],
    })
    spectra_file = _make_temp_spectra_file(spectra_table)

    try:
        with tempfile.TemporaryDirectory() as output_dir:
            with open(os.path.join(output_dir, 'existing.txt'), 'w', encoding='utf-8') as f:
                f.write('occupied')

            mock_compiler = Mock()
            mock_compiler.compile_spectra.return_value = [f"{output_dir}_1/compiled_chunk_001.fits"]
            mock_compiler.create_metadata_table.return_value = f"{output_dir}_1/compiled_metadata.fits"

            with patch('euclidkit.core.spectra.SpectrumCompiler', return_value=mock_compiler):
                with patch('euclidkit.utils.io.load_table', return_value=spectra_table):
                    result = runner.invoke(compile_spectra, [
                        '--spectra-table', spectra_file,
                        '--output-dir', output_dir,
                        '--output-format', 'fits',
                        '--prefix', 'compiled',
                    ])

            assert result.exit_code == 0
            mock_compiler.compile_spectra.assert_called_once()
            kwargs = mock_compiler.compile_spectra.call_args.kwargs
            assert kwargs['output_dir'] == output_dir
            assert kwargs['overwrite'] is False
            assert "Resume mode: existing chunk files will be kept." in result.output
    finally:
        if os.path.exists(spectra_file):
            os.unlink(spectra_file)


def test_compile_spectra_workers_passed_to_canonical_compile():
    """--workers should be forwarded to canonical compile_spectra."""
    runner = CliRunner()
    spectra_table = Table({
        'source_id': [100001, 100002],
        'datalabs_path': ['/tmp', '/tmp'],
        'file_name': ['a.fits', 'b.fits'],
        'hdu_index': [1, 1],
    })
    spectra_file = _make_temp_spectra_file(spectra_table)

    try:
        with tempfile.TemporaryDirectory() as output_dir:
            mock_compiler = Mock()
            mock_compiler.compile_spectra.return_value = [f"{output_dir}/compiled_chunk_001.fits"]
            mock_compiler.create_metadata_table.return_value = f"{output_dir}/compiled_metadata.fits"

            with patch('euclidkit.core.spectra.SpectrumCompiler', return_value=mock_compiler):
                with patch('euclidkit.utils.io.load_table', return_value=spectra_table):
                    result = runner.invoke(compile_spectra, [
                        '--spectra-table', spectra_file,
                        '--output-dir', output_dir,
                        '--output-format', 'fits',
                        '--workers', '2',
                    ])

            assert result.exit_code == 0
            mock_compiler.compile_spectra.assert_called_once()
            kwargs = mock_compiler.compile_spectra.call_args.kwargs
            assert kwargs['workers'] == 2
    finally:
        if os.path.exists(spectra_file):
            os.unlink(spectra_file)


def test_compile_spectra_non_deep_lambda_fallback_warns():
    """Non-DEEP canonical mode should warn and fall back to RGS when BGS/BOTH requested."""
    runner = CliRunner()
    spectra_table = Table({
        'source_id': [100001, 100002],
        'datalabs_path': ['/tmp', '/tmp'],
        'file_name': ['a.fits', 'b.fits'],
        'hdu_index': [1, 1],
    })
    spectra_file = _make_temp_spectra_file(spectra_table)

    try:
        with tempfile.TemporaryDirectory() as output_dir:
            mock_compiler = Mock()
            mock_compiler.compile_spectra.return_value = [f"{output_dir}/compiled_chunk_001.fits"]
            mock_compiler.create_metadata_table.return_value = f"{output_dir}/compiled_metadata.fits"

            with patch('euclidkit.core.spectra.SpectrumCompiler', return_value=mock_compiler):
                with patch('euclidkit.utils.io.load_table', return_value=spectra_table):
                    result = runner.invoke(compile_spectra, [
                        '--spectra-table', spectra_file,
                        '--output-dir', output_dir,
                        '--output-format', 'fits',
                        '--environment', 'PDR',
                        '-L', 'BOTH',
                    ])

            assert result.exit_code == 0
            assert "falling back to RGS" in result.output
            mock_compiler.compile_spectra.assert_called_once()
            assert not mock_compiler.annotate_lambda_range_from_xml.called
    finally:
        if os.path.exists(spectra_file):
            os.unlink(spectra_file)


def test_compile_spectra_idr_deep_both_splits_outputs():
    """IDR DEEP + -L BOTH should compile RGS and BGS separately."""
    runner = CliRunner()
    spectra_table = Table({
        'source_id': [100001, 100002, 100003],
        'datalabs_path': ['/tmp', '/tmp', '/tmp'],
        'file_name': ['a.fits', 'b.fits', 'c.fits'],
        'hdu_index': [1, 1, 1],
    })
    spectra_file = _make_temp_spectra_file(spectra_table)

    rgs_table = spectra_table[:2].copy()
    rgs_table['lambda_range'] = ['RGS', 'RGS']
    bgs_table = spectra_table[2:].copy()
    bgs_table['lambda_range'] = ['BGS']

    try:
        with tempfile.TemporaryDirectory() as output_dir:
            mock_compiler = Mock()
            annotated = spectra_table.copy()
            annotated['lambda_range'] = ['RGS', 'RGS', 'BGS']
            mock_compiler.annotate_lambda_range_from_xml.return_value = (
                annotated,
                {'total': 3, 'resolved': 3, 'rgs': 2, 'bgs': 1, 'unresolved': 0, 'ambiguous': 0},
            )
            mock_compiler.filter_table_by_lambda_range.side_effect = [
                (rgs_table, {'selected': 2, 'unresolved': 0, 'ambiguous': 0, 'skipped_other_band': 1}),
                (bgs_table, {'selected': 1, 'unresolved': 0, 'ambiguous': 0, 'skipped_other_band': 2}),
            ]
            mock_compiler.compile_spectra.side_effect = [
                [f"{output_dir}/myprefix_rgs_chunk_001.fits"],
                [f"{output_dir}/myprefix_bgs_chunk_001.fits"],
            ]
            mock_compiler.create_metadata_table.side_effect = [
                f"{output_dir}/myprefix_rgs_metadata.fits",
                f"{output_dir}/myprefix_bgs_metadata.fits",
            ]

            with patch('euclidkit.core.spectra.SpectrumCompiler', return_value=mock_compiler):
                with patch('euclidkit.utils.io.load_table', return_value=spectra_table):
                    result = runner.invoke(compile_spectra, [
                        '--spectra-table', spectra_file,
                        '--output-dir', output_dir,
                        '--output-format', 'fits',
                        '--prefix', 'myprefix',
                        '--environment', 'IDR',
                        '--idr-field', 'DEEP',
                        '-L', 'BOTH',
                    ])

            assert result.exit_code == 0
            assert "LambdaRange XML summary:" in result.output
            assert "LambdaRange selection summary (BOTH):" in result.output
            assert mock_compiler.compile_spectra.call_count == 2
            compile_kwargs_1 = mock_compiler.compile_spectra.call_args_list[0].kwargs
            compile_kwargs_2 = mock_compiler.compile_spectra.call_args_list[1].kwargs
            assert compile_kwargs_1['output_prefix'] == 'myprefix_rgs'
            assert compile_kwargs_2['output_prefix'] == 'myprefix_bgs'
    finally:
        if os.path.exists(spectra_file):
            os.unlink(spectra_file)


def test_compile_spectra_datalink_deduplicates_source_rows():
    """Datalink mode should deduplicate duplicate source_id rows before compile/metadata."""
    runner = CliRunner()
    spectra_table = Table({
        'source_id': [100001, 100001, 100002],
        'object_id': [100001, 100001, 100002],
    })
    spectra_file = _make_temp_spectra_file(spectra_table)

    try:
        with tempfile.TemporaryDirectory() as output_dir:
            dedup_table = Table({
                'source_id': [100001, 100002],
                'object_id': [100001, 100002],
            })
            mock_compiler = Mock()
            mock_compiler.deduplicate_by_source_id.return_value = (
                dedup_table,
                {'input_rows': 3, 'unique_sources': 2, 'duplicate_rows_removed': 1},
            )
            mock_compiler.compile_spectra_datalink.return_value = [f"{output_dir}/dl_compiled_chunk_001.fits"]
            mock_compiler.create_metadata_table.return_value = f"{output_dir}/dl_compiled_metadata.fits"

            mock_archive = Mock()
            mock_archive.euclid = Mock()

            with patch('euclidkit.core.spectra.SpectrumCompiler', return_value=mock_compiler):
                with patch('euclidkit.core.data_access.EuclidArchive', return_value=mock_archive):
                    with patch('euclidkit.utils.io.load_table', return_value=spectra_table):
                        result = runner.invoke(compile_spectra, [
                            '--spectra-table', spectra_file,
                            '--output-dir', output_dir,
                        '--output-format', 'fits',
                            '--use-datalink',
                        ])

            assert result.exit_code == 0
            assert "Datalink deduplication:" in result.output
            compile_kwargs = mock_compiler.compile_spectra_datalink.call_args.kwargs
            metadata_kwargs = mock_compiler.create_metadata_table.call_args.kwargs
            assert compile_kwargs['spectra_table'] is dedup_table
            assert metadata_kwargs['spectra_table'] is dedup_table
    finally:
        if os.path.exists(spectra_file):
            os.unlink(spectra_file)


def test_compile_spectra_datalink_lambda_range_maps_to_retrieval_type_all():
    """In datalink mode, -L BOTH should compile separate RGS and BGS outputs."""
    runner = CliRunner()
    spectra_table = Table({
        'source_id': [100001, 100002],
        'object_id': [100001, 100002],
    })
    spectra_file = _make_temp_spectra_file(spectra_table)

    try:
        with tempfile.TemporaryDirectory() as output_dir:
            mock_compiler = Mock()
            mock_compiler.deduplicate_by_source_id.return_value = (
                spectra_table,
                {'input_rows': 2, 'unique_sources': 2, 'duplicate_rows_removed': 0},
            )
            mock_compiler.compile_spectra_datalink.side_effect = [
                [f"{output_dir}/dl_compiled_rgs_chunk_001.fits"],
                [f"{output_dir}/dl_compiled_bgs_chunk_001.fits"],
            ]
            mock_compiler.create_metadata_table.side_effect = [
                f"{output_dir}/dl_compiled_rgs_metadata.fits",
                f"{output_dir}/dl_compiled_bgs_metadata.fits",
            ]
            mock_archive = Mock()
            mock_archive.euclid = Mock()

            with patch('euclidkit.core.spectra.SpectrumCompiler', return_value=mock_compiler):
                with patch('euclidkit.core.data_access.EuclidArchive', return_value=mock_archive):
                    with patch('euclidkit.utils.io.load_table', return_value=spectra_table):
                        result = runner.invoke(compile_spectra, [
                            '--spectra-table', spectra_file,
                            '--output-dir', output_dir,
                        '--output-format', 'fits',
                            '--use-datalink',
                            '-L', 'BOTH',
                        ])

            assert result.exit_code == 0
            assert mock_compiler.compile_spectra_datalink.call_count == 2
            kwargs_rgs = mock_compiler.compile_spectra_datalink.call_args_list[0].kwargs
            kwargs_bgs = mock_compiler.compile_spectra_datalink.call_args_list[1].kwargs
            assert kwargs_rgs['retrieval_type'] == 'SPECTRA_RGS'
            assert kwargs_bgs['retrieval_type'] == 'SPECTRA_BGS'
            assert "Datalink BOTH mode: writing separate RGS and BGS compiled files." in result.output
    finally:
        if os.path.exists(spectra_file):
            os.unlink(spectra_file)


def test_compile_spectra_datalink_conflict_lambda_vs_retrieval_prefers_lambda():
    """If both are explicit and conflict, lambda-range should take precedence."""
    runner = CliRunner()
    spectra_table = Table({
        'source_id': [100001],
        'object_id': [100001],
    })
    spectra_file = _make_temp_spectra_file(spectra_table)

    try:
        with tempfile.TemporaryDirectory() as output_dir:
            mock_compiler = Mock()
            mock_compiler.deduplicate_by_source_id.return_value = (
                spectra_table,
                {'input_rows': 1, 'unique_sources': 1, 'duplicate_rows_removed': 0},
            )
            mock_compiler.compile_spectra_datalink.return_value = [f"{output_dir}/dl_compiled_chunk_001.fits"]
            mock_compiler.create_metadata_table.return_value = f"{output_dir}/dl_compiled_metadata.fits"
            mock_archive = Mock()
            mock_archive.euclid = Mock()

            with patch('euclidkit.core.spectra.SpectrumCompiler', return_value=mock_compiler):
                with patch('euclidkit.core.data_access.EuclidArchive', return_value=mock_archive):
                    with patch('euclidkit.utils.io.load_table', return_value=spectra_table):
                        result = runner.invoke(compile_spectra, [
                            '--spectra-table', spectra_file,
                            '--output-dir', output_dir,
                        '--output-format', 'fits',
                            '--use-datalink',
                            '-L', 'BGS',
                            '--retrieval-type', 'SPECTRA_RGS',
                        ])

            assert result.exit_code == 0
            kwargs = mock_compiler.compile_spectra_datalink.call_args.kwargs
            assert kwargs['retrieval_type'] == 'SPECTRA_BGS'
            assert "disagree in datalink mode" in result.output
    finally:
        if os.path.exists(spectra_file):
            os.unlink(spectra_file)


def test_compile_spectra_default_uses_parquet_export(tmp_path):
    """Non-datalink compile-spectra should default to parquet with auto workers."""
    from types import SimpleNamespace

    runner = CliRunner()
    catalog = tmp_path / "catalog.fits"
    Table({'object_id': [1]}).write(catalog, format='fits', overwrite=True)

    stats = SimpleNamespace(
        requested_rows=1,
        exported_rows=1,
        skipped_rows=0,
        failed_rows=0,
        file_missing_rows=0,
        output_files=[str(tmp_path / 'compiled_part001.parquet')],
        manifest_path=str(tmp_path / 'compiled_manifest.json'),
        failures_path='',
    )

    with patch('euclidkit.core.spectra_parquet.default_raw_parquet_workers', return_value=7):
        with patch('euclidkit.core.spectra_parquet.spectra_to_parquet', return_value=stats) as mock_export:
            result = runner.invoke(compile_spectra, [
                '--spectra-table', str(catalog),
                '--output-dir', str(tmp_path),
                '--prefix', 'compiled',
            ])

    assert result.exit_code == 0
    assert "Parquet export completed successfully!" in result.output
    mock_export.assert_called_once()
    kwargs = mock_export.call_args.kwargs
    assert kwargs['catalog_table'] == str(catalog)
    assert kwargs['output_prefix'] == str(tmp_path / 'compiled')
    assert kwargs['chunk_size'] == 2000
    assert kwargs['workers'] == 7
    assert kwargs['lambda_range'] == 'RGS'
    assert kwargs['on_error'] == 'fail'


def test_compile_spectra_parquet_both_runs_two_arm_exports(tmp_path):
    """Parquet -L BOTH should write separate RGS/BGS prefix families."""
    from types import SimpleNamespace

    runner = CliRunner()
    catalog = tmp_path / "catalog.fits"
    raw_frame = tmp_path / "raw_frame.vot"
    Table({'object_id': [1]}).write(catalog, format='fits', overwrite=True)
    Table(
        {
            'pointing_id': [1],
            'technique': ['SPECTROIMAGE'],
            'obs_time_mjd': [60000.0],
            'obs_time_utc': ['2025-01-01T00:00:00Z'],
            'pa': [42.0],
        }
    ).write(raw_frame, format='votable', overwrite=True)
    stats = SimpleNamespace(
        requested_rows=1,
        exported_rows=1,
        skipped_rows=0,
        failed_rows=0,
        file_missing_rows=0,
        output_files=[],
        manifest_path=str(tmp_path / 'manifest.json'),
        failures_path='',
    )

    with patch('euclidkit.core.spectra_parquet.spectra_to_parquet', return_value=stats) as mock_export:
        result = runner.invoke(compile_spectra, [
            '--spectra-table', str(catalog),
            '--output-dir', str(tmp_path),
            '--prefix', 'compiled',
            '-L', 'BOTH',
            '--workers', '3',
            '--chunk-size', '10',
            '--on-error', 'skip',
        ])

    assert result.exit_code == 0
    assert "Parquet BOTH mode" in result.output
    assert mock_export.call_count == 2
    first = mock_export.call_args_list[0].kwargs
    second = mock_export.call_args_list[1].kwargs
    assert first['output_prefix'] == str(tmp_path / 'compiled_rgs')
    assert second['output_prefix'] == str(tmp_path / 'compiled_bgs')
    assert [first['lambda_range'], second['lambda_range']] == ['RGS', 'BGS']
    assert first['workers'] == second['workers'] == 3
    assert first['chunk_size'] == second['chunk_size'] == 10
    assert first['on_error'] == second['on_error'] == 'skip'


def test_dithers_to_parquet_main_cli_forwards_options(tmp_path):
    """The main euclidkit CLI should expose dithers-to-parquet."""
    from types import SimpleNamespace
    from euclidkit.cli.main import main

    runner = CliRunner()
    catalog = tmp_path / "catalog.fits"
    raw_frame = tmp_path / "raw_frame.vot"
    Table({'object_id': [1]}).write(catalog, format='fits', overwrite=True)
    Table(
        {
            'pointing_id': [1],
            'technique': ['SPECTROIMAGE'],
            'obs_time_mjd': [60000.0],
            'obs_time_utc': ['2025-01-01T00:00:00Z'],
            'pa': [42.0],
        }
    ).write(raw_frame, format='votable', overwrite=True)
    stats = SimpleNamespace(
        objects_requested=1,
        objects_exported=1,
        skipped_rows=0,
        failed_rows=0,
        file_missing_rows=0,
        combined_rows=0,
        dither_rows=2,
        combined_output_files=[],
        dither_output_files=[str(tmp_path / 'raw_dithers_part001.parquet')],
        manifest_path=str(tmp_path / 'raw_manifest.json'),
        failures_path='',
        dither_rows_with_raw_frame_metadata=1,
        dither_rows_missing_raw_frame_metadata=0,
    )

    with patch('euclidkit.core.spectra_parquet.dithers_to_parquet', return_value=stats) as mock_export:
        result = runner.invoke(main, [
            'dithers-to-parquet',
            '--catalog-table', str(catalog),
            '--output-prefix', str(tmp_path / 'raw'),
            '--chunk-size', '11',
            '--workers', '2',
            '--lambda-range', 'RGS',
            '--raw-frame-table', str(raw_frame),
            '--dithers-only',
            '--overwrite',
            '--on-error', 'skip',
            '--no-progress',
        ])

    assert result.exit_code == 0
    assert "Dither parquet export completed successfully!" in result.output
    mock_export.assert_called_once()
    kwargs = mock_export.call_args.kwargs
    assert kwargs['catalog_table'] == str(catalog)
    assert kwargs['output_prefix'] == str(tmp_path / 'raw')
    assert kwargs['chunk_size'] == 11
    assert kwargs['workers'] == 2
    assert kwargs['lambda_range'] == 'RGS'
    assert kwargs['raw_frame_table'] == str(raw_frame)
    assert kwargs['include_combined'] is False
    assert kwargs['overwrite'] is True
    assert kwargs['on_error'] == 'skip'
    assert kwargs['show_progress'] is False
