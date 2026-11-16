"""Public releases and internal services must never be conflated."""
import json
from types import SimpleNamespace
from unittest.mock import Mock, patch

import pandas as pd
import pytest
from astropy.table import Table

from euclidkit.core.archive_service import product_schema_for_client
from euclidkit.core.cutana_catalogue import download_alternatives
from euclidkit.core.cutouts import CutoutGenerator
from euclidkit.core.data_access import EuclidArchive
from euclidkit.core.spectra_parquet import _raw_frame_table_name


@pytest.mark.parametrize('release', [None, 'DR1', 'dr1'])
def test_public_defaults_to_dr1(release):
    archive = EuclidArchive(data_release=release)
    assert archive.environment == 'PDR' and archive.data_release == 'DR1'
    assert archive._get_mer_table_names() == ['dr1.mer_catalogue_wide_survey', 'dr1.mer_catalogue_wide_mode']
    assert archive._get_spectra_source_table_name('DEEP') == 'dr1.spectra_source_deep'
    assert archive._get_segmentation_map_table_name() == 'dr1.mer_segmentation_map'
    assert archive._get_product_table_name('mosaic_product') == 'dr1.mosaic_product'
    assert _raw_frame_table_name('PDR', release) == 'dr1.raw_frame'


def test_q1_has_no_dr1_routes():
    archive = EuclidArchive(data_release='Q1')
    names = [*archive._get_mer_table_names(), archive._get_spectra_source_table_name(),
             archive._get_segmentation_map_table_name(), archive._get_product_table_name('mosaic_product'),
             _raw_frame_table_name('PDR', 'Q1')]
    assert all(name.startswith('q1.') for name in names)
    for operation in [lambda: archive._get_zspe_candidate_table_names('qso'),
                      lambda: archive.service.morphology_tables(),
                      lambda: archive._get_mer_table_names('DEEP')]:
        with pytest.raises(ValueError):
            operation()


@pytest.mark.parametrize('environment', ['IDR', 'OTF', 'REG'])
@pytest.mark.parametrize('release', ['DR1', 'Q1'])
def test_public_release_on_other_service_is_rejected(environment, release):
    with pytest.raises(ValueError, match='require environment=PDR'):
        EuclidArchive(environment, data_release=release)


@pytest.mark.parametrize('environment', ['PDR', 'IDR', 'OTF', 'REG'])
def test_client_stays_on_selected_service(environment):
    archive = EuclidArchive(environment)
    with patch('astroquery.esa.euclid.core.EuclidClass') as client:
        archive._ensure_client()
    client.assert_called_once_with(environment=environment)
    assert archive.data_release == ('DR1' if environment == 'PDR' else None)


@pytest.mark.parametrize('partition,expected', [('survey', ['survey']), ('mode', ['mode']), ('both', ['survey', 'mode'])])
@pytest.mark.parametrize('environment,prefix', [('PDR', 'dr1'), ('IDR', 'catalogue')])
def test_deep_partitions(environment, prefix, partition, expected):
    archive = EuclidArchive(environment)
    assert archive._get_mer_table_names('DEEP', partition) == [f'{prefix}.mer_catalogue_deep_{p}' for p in expected]
    assert archive.service.morphology_tables('DEEP', partition) == [f'{prefix}.mer_morphology_deep_{p}' for p in expected]
    assert archive._get_zspe_candidate_table_names('qso', 'DEEP') == [(f'{prefix}.spectro_z_spe_qso_candidates_deep', 'deep')]


def test_public_morphology_uses_dedicated_tables():
    from tests.test_core.test_morphology import jobs, result
    archive = EuclidArchive()
    archive._logged_in = True
    archive.euclid = Mock()
    calls = jobs(archive, [result([1]), result([2])])
    matched = archive.query_morphology(Table({'object_id': [1, 2]}))
    assert len(matched) == 2
    assert 'dr1.mer_morphology_wide_survey' in calls[0][0]
    assert 'dr1.mer_morphology_wide_mode' in calls[1][0]
    assert list(calls[1][1]['object_id']) == [2]


def test_public_morphology_recovery_records_release(tmp_path):
    from tests.test_core.test_morphology import jobs
    archive = EuclidArchive()
    archive._logged_in = True
    archive.euclid = Mock()
    jobs(archive, [RuntimeError('download interrupted')], asynchronous=True)
    output = tmp_path/'recovery.fits'
    metadata = archive.query_morphology(Table({'object_id': [1]}), output, full_async=True)
    assert metadata['environment'] == 'PDR' and metadata['data_release'] == 'DR1'
    record = json.loads((tmp_path/'recovery.fits.job.json').read_text())
    assert record['data_release'] == 'DR1'
    assert all(name.startswith('dr1.') for name in record['tables'])


def test_public_crossmatch_retains_survey_first_fallback():
    archive = EuclidArchive()
    archive._logged_in = True
    archive.euclid = Mock()
    first = Table({'object_id': [1], 'object_id_user': [1]})
    second = Table({'object_id': [2], 'object_id_user': [2]})
    with patch.object(archive, '_crossmatch_batch', side_effect=[first, second]) as query:
        result = archive.crossmatch_sources(Table({'object_id': [1, 2]}), use_object_id=True)
    assert len(result) == 2
    assert query.call_args_list[0].args[4] == 'dr1.mer_catalogue_wide_survey'
    assert query.call_args_list[1].args[4] == 'dr1.mer_catalogue_wide_mode'
    assert list(query.call_args_list[1].args[0]['object_id']) == [2]


def test_custom_mer_override_remains_single_table():
    archive = EuclidArchive()
    archive._logged_in = True
    archive.euclid = Mock()
    with patch.object(archive, '_crossmatch_batch', return_value=Table({'object_id': [1]})) as query:
        archive.crossmatch_sources(Table({'object_id': [1]}), mer_table='custom.my_mer', use_object_id=True)
    assert query.call_count == 1 and query.call_args.args[4] == 'custom.my_mer'


@pytest.mark.parametrize('release,prefix', [(None, 'dr1'), ('Q1', 'q1')])
def test_cutana_query_and_manifest_carry_release(tmp_path, release, prefix):
    archive = EuclidArchive(data_release=release)
    archive.euclid = Mock()
    archive._logged_in = True
    source = pd.DataFrame({'object_id': [1], 'ra': [10.], 'dec': [2.], 'segmentation_map_id': [1000000]})
    files = source.assign(datalabs_path='/data', file_name='VIS.fits', filter_name='VIS')
    archive.euclid.launch_job_async.return_value.get_results.return_value = Table.from_pandas(files)
    generator = CutoutGenerator(archive=archive)
    with patch.object(generator, '_resolve_coordinates', return_value=source):
        generator.generate_cutana_input(source, tmp_path/'sources.csv')
    query = archive.euclid.launch_job_async.call_args.args[0]
    assert f'JOIN {prefix}.mosaic_product' in query
    manifest = json.loads((tmp_path/'sources.manifest.json').read_text())
    assert manifest['archive_environment'] == 'PDR' and manifest['data_release'] == prefix.upper()


def test_alternative_download_rejects_wrong_release(tmp_path):
    path = tmp_path/'alternatives.csv'
    pd.DataFrame({'selected': [False], 'file_name': ['a.fits'], 'data_release': ['Q1']}).to_csv(path, index=False)
    archive = EuclidArchive()
    archive.euclid = Mock()
    with pytest.raises(ValueError, match='matching service/release'):
        download_alternatives(archive, path, tmp_path/'files')
    archive.euclid.get_product.assert_not_called()


def test_product_schema_defaults_and_override():
    assert product_schema_for_client(SimpleNamespace(environment='PDR')) == 'dr1'
    assert product_schema_for_client(SimpleNamespace(environment='PDR'), data_release='Q1') == 'q1'
    assert product_schema_for_client(SimpleNamespace(environment='IDR')) == 'sedm'
    assert product_schema_for_client(SimpleNamespace(environment='PDR'), schema='custom') == 'custom'


@pytest.mark.parametrize('environment,release', [('PDR', None), ('PDR', 'Q1'), ('IDR', None)])
def test_client_environment_is_never_redirected(monkeypatch, environment, release):
    import sys
    from unittest.mock import Mock
    factory = Mock()
    monkeypatch.setitem(sys.modules, 'astroquery.esa.euclid.core', SimpleNamespace(EuclidClass=factory))
    archive = EuclidArchive(environment=environment, data_release=release)
    assert archive._ensure_client() is factory.return_value
    factory.assert_called_once_with(environment=environment)


@pytest.mark.parametrize('release,prefix', [(None, 'dr1'), ('Q1', 'q1')])
def test_segmentation_and_raw_frame_queries_use_public_release(release, prefix):
    from unittest.mock import Mock

    from euclidkit.core.spectra_parquet import _query_raw_frame_lookup
    archive = EuclidArchive(data_release=release)
    archive._logged_in = True
    archive.euclid = Mock()
    archive.euclid.launch_job.return_value.get_results.return_value = Table()
    archive.query_segmentation_maps(Table({'object_id': [1], 'segmentation_map_id': [1000001],
                                           'ra': [1.], 'dec': [2.]}))
    query = archive.euclid.launch_job.call_args.args[0]
    assert f'JOIN {prefix}.mer_segmentation_map' in query
    archive.euclid.launch_job.return_value.get_results.return_value = Table({
        'pointing_id': [123], 'grism_wheel_pos': ['RGS000'], 'obs_time_mjd': [60000.],
        'obs_time_utc': ['2023-02-25T00:00:00'], 'pa': [0.]})
    _query_raw_frame_lookup([{'ptgid': 123, 'gwa_pos': 'RGS000'}], archive=archive,
                            environment='PDR', data_release=release)
    query = archive.euclid.launch_job.call_args.args[0]
    assert f'FROM {prefix}.raw_frame' in query
    assert 'r.grism_wheel_pos = p.gwa_pos' in query


def test_public_spe_queries_use_dr1_and_only_fallback_unmatched():
    archive = EuclidArchive()
    archive._logged_in = True
    archive.euclid = Mock()
    captures = []

    def launch(query, upload_resource, upload_table_name):
        uploaded = Table.read(upload_resource)
        captures.append((query, list(uploaded['object_id'])))
        oid = len(captures)
        response = Table({'object_id': [oid], 'spe_rank': [1], 'spe_z': [1.2],
                          'spe_z_err': [0.1], 'source_table': ['wide_survey' if oid == 1 else 'wide']})
        job = Mock()
        job.get_results.return_value = response
        return job

    archive.euclid.launch_job.side_effect = launch
    result = archive.query_zspe_candidates(Table({'object_id': [1, 2]}))
    assert len(result) == 2
    assert 'dr1.spectro_z_spe_qso_candidates_wide_survey' in captures[0][0]
    assert 'dr1.spectro_z_spe_qso_candidates_wide ' in captures[1][0]
    assert captures[0][1] == [1, 2] and captures[1][1] == [2]


def test_q1_spe_python_api_rejects_before_authentication():
    archive = EuclidArchive(data_release='Q1')
    archive.login = Mock()
    with pytest.raises(ValueError, match='Q1 is not supported'):
        archive.query_zspe_candidates(Table({'object_id': [1]}))
    archive.login.assert_not_called()
