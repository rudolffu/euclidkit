"""Tests for Cutana segmentation-map shortcut logic in cutouts module."""

import os
import tempfile
from unittest.mock import Mock

import pandas as pd
from astropy.table import Table

from euclidkit.core.cutouts import CutoutGenerator
from euclidkit.core.data_access import EuclidArchive


def _build_generator_with_archive():
    archive = Mock()
    archive._logged_in = True
    archive.environment = 'IDR'
    archive._get_mer_table_name.return_value = 'catalogue.mer_catalogue_deep_survey'
    archive._get_mer_table_names.return_value = ['catalogue.mer_catalogue_deep_survey']
    archive.login.return_value = None
    archive._ensure_client.return_value = None
    archive._combine_mer_results.side_effect = lambda tables: tables[0] if tables else Table()
    archive.euclid = Mock()
    return CutoutGenerator(archive=archive), archive


def test_lookup_mer_source_info_uses_mer_table_selection():
    """MER lookup should use upload+join table selection and include segmentation_map_id."""
    generator, archive = _build_generator_with_archive()

    input_df = pd.DataFrame({'object_id': [100001, 100002]})
    mer_result = Table(
        {
            'object_id': ['100001', '100002'],
            'right_ascension': [150.0, 151.0],
            'declination': [2.0, 2.1],
            'det_quality_flag': [0, 1],
            'segmentation_map_id': [1234567890000, 9876543210000],
        }
    )
    mock_job = Mock()
    mock_job.get_results.return_value = mer_result
    archive.euclid.launch_job_async.return_value = mock_job

    resolved = generator._lookup_mer_source_info(input_df, idr_field='DEEP')

    archive._get_mer_table_names.assert_called_once_with(
        idr_field='DEEP',
        idr_deep_partition='survey',
    )
    query = archive.euclid.launch_job_async.call_args[0][0]
    assert 'FROM TAP_UPLOAD.' in query
    assert 'JOIN catalogue.mer_catalogue_deep_survey AS m ON u.object_id = m.object_id' in query
    assert 'segmentation_map_id' in query
    for field in [
        'det_quality_flag',
        'parent_id',
        'spurious_flag',
        'vis_det',
        'flag_vis',
        'flag_y',
        'flag_j',
        'flag_h',
    ]:
        assert f'm.{field}' in query
    assert 'det_quality_flag' in resolved.columns
    assert 'segmentation_map_id' in resolved.columns
    assert len(resolved) == 2


def test_lookup_mer_source_info_can_query_both_deep_partitions():
    """IDR DEEP Cutana lookup should support survey+mode partition selection."""
    generator, archive = _build_generator_with_archive()
    archive._get_mer_table_names.return_value = [
        'catalogue.mer_catalogue_deep_survey',
        'catalogue.mer_catalogue_deep_mode',
    ]
    archive._combine_mer_results.side_effect = EuclidArchive._combine_mer_results

    input_df = pd.DataFrame({'object_id': [100001, 100002]})
    survey_result = Table(
        {
            'object_id': ['100001'],
            'right_ascension': [150.0],
            'declination': [2.0],
            'segmentation_map_id': [1234567890000],
        }
    )
    mode_result = Table(
        {
            'object_id': ['100002'],
            'right_ascension': [151.0],
            'declination': [2.1],
            'segmentation_map_id': [9876543210000],
        }
    )
    survey_job = Mock()
    survey_job.get_results.return_value = survey_result
    mode_job = Mock()
    mode_job.get_results.return_value = mode_result
    archive.euclid.launch_job_async.side_effect = [survey_job, mode_job]

    resolved = generator._lookup_mer_source_info(
        input_df,
        idr_field='DEEP',
        idr_deep_partition='both',
    )

    archive._get_mer_table_names.assert_called_once_with(
        idr_field='DEEP',
        idr_deep_partition='both',
    )
    queries = [call.args[0] for call in archive.euclid.launch_job_async.call_args_list]
    assert 'JOIN catalogue.mer_catalogue_deep_survey AS m' in queries[0]
    assert 'JOIN catalogue.mer_catalogue_deep_mode AS m' in queries[1]
    assert len(resolved) == 2
    assert set(resolved['object_id']) == {'100001', '100002'}


def test_lookup_mer_source_info_queries_mode_only_for_unmatched_wide_ids():
    """IDR WIDE Cutana lookup should use wide_mode only for IDs missing in wide_survey."""
    generator, archive = _build_generator_with_archive()
    archive._get_mer_table_names.return_value = [
        'catalogue.mer_catalogue_wide_survey',
        'catalogue.mer_catalogue_wide_mode',
    ]
    archive._combine_mer_results.side_effect = EuclidArchive._combine_mer_results

    input_df = pd.DataFrame({'object_id': [100001, 100002]})
    survey_result = Table(
        {
            'object_id': ['100001'],
            'right_ascension': [150.0],
            'declination': [2.0],
            'segmentation_map_id': [1234567890000],
        }
    )
    mode_result = Table(
        {
            'object_id': ['100002'],
            'right_ascension': [151.0],
            'declination': [2.1],
            'segmentation_map_id': [9876543210000],
        }
    )
    survey_job = Mock()
    survey_job.get_results.return_value = survey_result
    mode_job = Mock()
    mode_job.get_results.return_value = mode_result
    archive.euclid.launch_job_async.side_effect = [survey_job, mode_job]

    resolved = generator._lookup_mer_source_info(input_df, idr_field='WIDE')

    assert archive.euclid.launch_job_async.call_count == 2
    queries = [call.args[0] for call in archive.euclid.launch_job_async.call_args_list]
    assert 'JOIN catalogue.mer_catalogue_wide_survey AS m' in queries[0]
    assert 'JOIN catalogue.mer_catalogue_wide_mode AS m' in queries[1]
    assert len(resolved) == 2
    assert set(resolved['object_id']) == {'100001', '100002'}


def test_generate_cutana_input_prefers_segmentation_shortcut():
    """Cutana input should use segmentation shortcut when segmentation_map_id is present."""
    generator, _ = _build_generator_with_archive()
    generator._resolve_coordinates = Mock(
        return_value=pd.DataFrame(
            {
                'object_id': [100001],
                'right_ascension': [150.0],
                'declination': [2.0],
                'segmentation_map_id': [1234567890000],
            }
        )
    )
    generator._get_files_mosaic_from_segmentation = Mock(
        return_value=pd.DataFrame(
            {
                'object_id': [100001],
                'right_ascension': [150.0],
                'declination': [2.0],
                'segmentation_map_id': [1234567890000],
                'datalabs_path': ['/data/euclid_q1'],
                'file_name': ['mock.fits'],
                'filter_name': ['VIS'],
                'dist': [0.0],
            }
        )
    )
    generator._get_files_mosaic = Mock(side_effect=AssertionError("positional path should not be used"))

    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as output_file:
        try:
            result = generator.generate_cutana_input(
                sources=pd.DataFrame({'object_id': [100001]}),
                output_file=output_file.name,
            )
        finally:
            os.unlink(output_file.name)

    generator._get_files_mosaic_from_segmentation.assert_called_once()
    assert len(result) == 1
    assert 'diameter_arcsec' in result.columns


def test_generate_cutana_input_raises_when_shortcut_empty():
    """Cutana input should fail when segmentation shortcut has no mosaic matches."""
    generator, _ = _build_generator_with_archive()
    coords_df = pd.DataFrame(
        {
            'object_id': [100001],
            'right_ascension': [150.0],
            'declination': [2.0],
            'segmentation_map_id': [1234567890000],
        }
    )
    generator._resolve_coordinates = Mock(return_value=coords_df)
    generator._get_files_mosaic_from_segmentation = Mock(return_value=pd.DataFrame())

    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as output_file:
        try:
            try:
                generator.generate_cutana_input(
                    sources=pd.DataFrame({'object_id': [100001]}),
                    output_file=output_file.name,
                )
                assert False, "Expected ValueError when no segmentation-based mosaic matches are found"
            except ValueError as exc:
                assert "No mosaic files matched" in str(exc)
        finally:
            os.unlink(output_file.name)

    generator._get_files_mosaic_from_segmentation.assert_called_once()


def test_lookup_mer_source_info_by_position_uses_crossmatch():
    """Coordinate-only inputs should resolve segmentation_map_id via crossmatch_sources."""
    generator, archive = _build_generator_with_archive()
    input_df = pd.DataFrame({'ra': [150.0], 'dec': [2.0]})
    crossmatch_result = Table(
        {
            'ra': [150.0],
            'dec': [2.0],
            'object_id': [100001],
            'mer_ra': [150.0001],
            'mer_dec': [2.0001],
            'segmentation_map_id': [1234567890000],
            'separation_arcsec': [0.2],
        }
    )
    archive.crossmatch_sources.return_value = crossmatch_result

    resolved = generator._lookup_mer_source_info_by_position(input_df, idr_field='DEEP')

    archive.crossmatch_sources.assert_called_once()
    kwargs = archive.crossmatch_sources.call_args.kwargs
    assert kwargs['use_object_id'] is False
    assert kwargs['idr_field'] == 'DEEP'
    assert kwargs['idr_deep_partition'] == 'survey'
    assert 'segmentation_map_id' in resolved.columns
    assert len(resolved) == 1


def test_get_files_mosaic_from_segmentation_uses_computed_tile_index():
    """Segmentation shortcut should upload floor(id / 1e6) tile_index and join on it."""
    generator, archive = _build_generator_with_archive()
    captured = {}
    job = Mock()
    job.get_results.return_value = Table({
        'object_id': [100001],
        'right_ascension': [150.0],
        'declination': [2.0],
        'segmentation_map_id': [1234567890000],
        'file_name': ['mock.fits'],
        'file_path': ['/data/euclid/mock.fits'],
        'datalabs_path': ['/data/euclid'],
        'mosaic_product_oid': [1],
        'tile_index': [1234567],
        'instrument_name': ['VIS'],
        'filter_name': ['VIS'],
        'image_ra': [150.0],
        'image_dec': [2.0],
    })

    def fake_launch_job_async(query, upload_resource=None, upload_table_name=None):
        captured['query'] = query
        captured['upload_table'] = Table.read(upload_resource, format='votable')
        captured['upload_table_name'] = upload_table_name
        return job

    archive.euclid.launch_job_async.side_effect = fake_launch_job_async
    sources = pd.DataFrame({
        'object_id': [100001, 100002, 100003],
        'right_ascension': [150.0, 151.0, 152.0],
        'declination': [2.0, 2.1, 2.2],
        'segmentation_map_id': [1234567890000, None, 'bad'],
    })

    result = generator._get_files_mosaic_from_segmentation(
        sources,
        instrument_name='VIS',
        nisp_filters=None,
        radec_colnames={'ra_colname': 'right_ascension', 'dec_colname': 'declination'},
    )

    query = captured['query']
    assert 'mosaics.tile_index = u.tile_index' in query
    assert 'SUBSTRING' not in query
    assert list(captured['upload_table']['tile_index']) == [1234567]
    assert [str(value) for value in captured['upload_table']['segmentation_map_id']] == [
        '1234567890000'
    ]
    assert len(result) == 1
