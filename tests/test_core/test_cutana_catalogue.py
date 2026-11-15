"""Cutana preparation contracts without network access or a Cutana dependency."""
import json
from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest

from euclidkit.core.cutana_catalogue import (
    download_alternatives, validate_size, write_catalogue,
)
from euclidkit.core.cutouts import CutoutGenerator


def inputs():
    coords = pd.DataFrame({'object_id': [1, 2, 3], 'ra': [359.99, 20., 30.],
                           'dec': [0., 1., 2.], 'segmentation_map_id': [1, 2, 3]})
    files = pd.DataFrame({
        'object_id': [1, 1, 1, 1, 2], 'ra': [359.99]*4 + [20.],
        'dec': [0.]*4 + [1.], 'filter_name': ['NIR_H', 'NIR_H', 'NIR_H', 'NIR_Y', 'NIR_Y'],
        'file_name': ['far_H.fits', 'near_H.fits', 'near_H.fits', 'one_Y.fits', 'two_Y.fits'],
        'datalabs_path': ['/data']*5, 'file_path': ['/archive']*5,
        'mosaic_product_oid': [1, 2, 2, 3, 4],
        'image_ra': [359., 0., 0., 359.99, 20.], 'image_dec': [0., 0., 0., 0., 1.],
    })
    return files, coords


def prepare(tmp_path, files=None, coords=None, keep=False):
    f, c = inputs()
    return write_catalogue(f if files is None else files, c if coords is None else coords,
                           tmp_path/'sources.csv', 'ra', 'dec', ['NIR_Y', 'NIR_H'],
                           'arcsec', 15., keep)


def test_selection_coverage_and_group_catalogues(tmp_path):
    result = prepare(tmp_path)
    assert result.SourceID.tolist() == ['OBJ_1', 'OBJ_2']
    assert json.loads(result.iloc[0].fits_file_paths) == ['/data/one_Y.fits', '/data/near_H.fits']
    assert json.loads(result.iloc[1].fits_file_paths) == ['/data/two_Y.fits']
    manifest = json.loads((tmp_path/'sources.manifest.json').read_text())
    assert manifest['sources_without_mosaics'] == 1
    assert manifest['sources_with_alternatives'] == 1
    assert manifest['combined_catalogue_requires_separate_group_runs']
    assert sorted(g['bands'] for g in manifest['groups']) == [['NIR_Y'], ['NIR_Y', 'NIR_H']]
    for group in manifest['groups']:
        table = pd.read_csv(group['catalogue'])
        assert len(table) == group['source_count']
        assert all(len(json.loads(v)) == len(group['bands']) for v in table.fits_file_paths)
    alternatives = pd.read_csv(manifest['alternatives'])
    assert len(alternatives) == 2  # exact duplicate records do not count
    assert alternatives.loc[alternatives.selected, 'file_name'].tolist() == ['near_H.fits']
    assert {'object_id', 'mosaic_product_oid', 'file_path', 'candidate_rank'} <= set(alternatives)
    coverage = pd.read_csv(manifest['coverage'])
    assert len(coverage) == 3
    assert not coverage.loc[coverage.object_id.eq(3), 'has_mosaic'].iloc[0]


def test_candidate_shuffle_does_not_change_selection(tmp_path):
    files, coords = inputs()
    a = prepare(tmp_path/'a', files, coords)
    b = prepare(tmp_path/'b', files.sample(frac=1, random_state=3), coords)
    pd.testing.assert_frame_equal(a, b)  # attrs do not participate


def test_ties_and_missing_centres_use_path_order(tmp_path):
    files, coords = inputs()
    files = files.drop(columns=['image_ra', 'image_dec'])
    result = prepare(tmp_path, files, coords)
    assert '/data/far_H.fits' in json.loads(result.iloc[0].fits_file_paths)


def test_duplicate_ids_positions_and_metadata(tmp_path):
    files, coords = inputs()
    # Same object at different requested positions must retain distinct IDs.
    files.loc[files.object_id.eq(2), 'object_id'] = 1
    coords.loc[coords.object_id.eq(2), 'object_id'] = 1
    coords['label'] = ['first', 'second', 'unmatched']
    coords = pd.concat([coords, coords.iloc[:1]], ignore_index=True)
    result = prepare(tmp_path, files, coords, keep=True)
    assert len(result) == 2 and result.SourceID.is_unique
    assert result.SourceID.str.startswith('OBJ_1_').all()
    assert set(result.label) == {'first', 'second'}


def test_coordinate_ids_do_not_remove_decimal_points(tmp_path):
    files, coords = inputs()
    files, coords = files.drop(columns='object_id'), coords.drop(columns='object_id')
    result = prepare(tmp_path, files, coords)
    assert result.SourceID.is_unique


@pytest.mark.parametrize('unit,value', [('arcsec', -1), ('arcsec', 0), ('arcsec', 3601),
    ('arcsec', np.nan), ('arcsec', np.inf), ('pixel', 1.5), ('pixel', 10001)])
def test_invalid_size_rejected_before_archive_query(unit, value, tmp_path):
    generator = CutoutGenerator(archive=Mock())
    generator._resolve_coordinates = Mock()
    with pytest.raises(ValueError):
        generator.generate_cutana_input({'object_id': 1}, tmp_path/'bad.csv',
                                        cutout_size=unit, cutout_size_value=value)
    generator._resolve_coordinates.assert_not_called()
    assert not (tmp_path/'bad.csv').exists()


def test_size_boundaries():
    assert validate_size('pixel', 10000.) == 10000
    assert isinstance(validate_size('pixel', 15.), int)
    assert validate_size('arcsec', 3600.) == 3600.


def test_single_band_set_uses_original_catalogue(tmp_path):
    files, coords = inputs()
    result = prepare(tmp_path, files[files.object_id.eq(1)], coords.iloc[:1])
    groups = result.attrs['cutana_groups']
    assert len(groups) == 1
    assert Path(groups[0]['catalogue']).name == 'sources.csv'


def test_invalid_coordinates_do_not_write_output(tmp_path):
    files, coords = inputs()
    coords.loc[0, 'ra'] = np.nan
    with pytest.raises(ValueError, match='Invalid ra'):
        prepare(tmp_path, files, coords)
    assert not (tmp_path/'sources.csv').exists()


def test_download_alternatives_deduplicates_and_resumes(tmp_path):
    prepare(tmp_path)
    archive = Mock()
    def get_product(**kwargs):
        path = Path(kwargs['output_file']); path.write_bytes(b'FITS content'); return [str(path)]
    archive.euclid.get_product.side_effect = get_product
    table = pd.read_csv(tmp_path/'sources.alternatives.csv')
    pd.concat([table, table]).to_csv(tmp_path/'duplicates.csv', index=False)
    paths = download_alternatives(archive, tmp_path/'duplicates.csv', tmp_path/'download', 'dr1')
    assert len(paths) == 1 and Path(paths[0]).name == 'far_H.fits'
    assert archive.euclid.get_product.call_count == 1
    assert archive.euclid.get_product.call_args.kwargs['schema'] == 'dr1'
    download_alternatives(archive, tmp_path/'duplicates.csv', tmp_path/'download')
    assert archive.euclid.get_product.call_count == 1


def test_failed_download_does_not_leave_resumable_partial_file(tmp_path):
    prepare(tmp_path)
    archive = Mock()
    def fail(**kwargs):
        Path(kwargs['output_file']).write_bytes(b'partial')
        return None
    archive.euclid.get_product.side_effect = fail
    with pytest.raises(RuntimeError):
        download_alternatives(archive, tmp_path/'sources.alternatives.csv', tmp_path/'download')
    assert not list((tmp_path/'download').iterdir())


def test_segmentation_query_checks_source_centre_footprint():
    generator = CutoutGenerator(archive=Mock())
    generator.archive._logged_in = True
    generator.archive.euclid.launch_job_async.return_value.get_results.return_value.to_pandas.return_value = pd.DataFrame()
    generator._get_files_mosaic_from_segmentation(
        pd.DataFrame({'ra': [1.], 'dec': [2.], 'segmentation_map_id': [1000000]}),
        'VIS', None, {'ra_colname': 'ra', 'dec_colname': 'dec'})
    query = generator.archive.euclid.launch_job_async.call_args.args[0]
    assert "CONTAINS(POINT('ICRS', u.ra, u.dec), mosaics.fov)=1" in query
    assert 'mosaics.fov IS NOT NULL' in query


def test_cutana_accepts_generated_groups_and_rejects_missing_configured_channel(tmp_path):
    cutana = pytest.importorskip('cutana')
    from cutana.fits_paths import parse_fits_file_paths
    from cutana.catalogue_preprocessor import validate_catalogue_columns, validate_coordinate_ranges
    from cutana.image_processor import combine_channels
    result = prepare(tmp_path)
    for group in result.attrs['cutana_groups']:
        frame = pd.read_csv(group['catalogue'])
        assert not validate_catalogue_columns(frame)
        assert not validate_coordinate_ranges(frame)
        assert all(len(parse_fits_file_paths(value)) == len(group['bands'])
                   for value in frame.fits_file_paths)
    with pytest.raises(ValueError):
        combine_channels(np.ones((1, 2, 2, 1), dtype='f4'),
                         {'NIR-Y': [1.], 'NIR-H': [1.]}, channel_names=['NIR-Y'])


def test_empty_alternatives_download_is_a_noop(tmp_path):
    files, coords = inputs()
    prepare(tmp_path, files[files.object_id.eq(2)], coords.iloc[1:2])
    archive = Mock()
    assert download_alternatives(archive, tmp_path/'sources.alternatives.csv', tmp_path/'download') == []
    archive.euclid.get_product.assert_not_called()


def test_download_rejects_parent_path(tmp_path):
    table = tmp_path/'alternatives.csv'
    pd.DataFrame({'selected': [False], 'file_name': ['../escape.fits']}).to_csv(table, index=False)
    archive = Mock()
    with pytest.raises(ValueError, match='basename'):
        download_alternatives(archive, table, tmp_path/'download')
    archive.euclid.get_product.assert_not_called()


def test_metadata_does_not_multiply_objects_at_same_position(tmp_path):
    files, coords = inputs()
    files.loc[files.object_id.eq(2), ['ra', 'dec']] = [359.99, 0.]
    coords.loc[coords.object_id.eq(2), ['ra', 'dec']] = [359.99, 0.]
    coords['label'] = ['one', 'two', 'three']
    result = prepare(tmp_path, files, coords, keep=True)
    assert len(result) == 2
    assert dict(zip(result.object_id, result.label)) == {1: 'one', 2: 'two'}


def test_large_catalogue_keeps_partial_sources_without_per_source_queries(tmp_path):
    count = 2000
    coords = pd.DataFrame({'object_id': np.arange(count), 'ra': np.arange(count)/100.,
                           'dec': np.zeros(count)})
    files = pd.concat([coords.assign(filter_name='NIR_Y', file_name='Y.fits'),
                       coords.iloc[::2].assign(filter_name='NIR_H', file_name='H.fits')],
                      ignore_index=True).assign(datalabs_path='/data')
    result = prepare(tmp_path, files.sample(frac=1, random_state=42), coords)
    assert len(result) == count and result.SourceID.is_unique
    assert sorted(g['source_count'] for g in result.attrs['cutana_groups']) == [1000, 1000]
