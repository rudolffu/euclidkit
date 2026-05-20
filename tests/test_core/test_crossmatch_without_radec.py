"""Ensure object-id mode does not require RA/Dec columns."""

import json
from unittest.mock import Mock, patch

from astropy.table import Table

from euclidkit.core.data_access import EuclidArchive


def test_idr_wide_mer_table_names_include_survey_and_mode():
    arch = EuclidArchive(environment='IDR')

    assert arch._get_mer_table_names(idr_field='WIDE') == [
        'catalogue.mer_catalogue_wide_survey',
        'catalogue.mer_catalogue_wide_mode',
    ]
    assert arch._get_mer_table_name(idr_field='WIDE') == 'catalogue.mer_catalogue_wide_survey'
    assert arch._get_mer_table_names(idr_field='DEEP') == ['catalogue.mer_catalogue_deep']


def test_crossmatch_sources_object_id_only_table(monkeypatch):
    # Table lacks 'ra'/'dec' but has object_id_euclid and euclid RA/Dec columns
    user_table = Table({
        'object_id_euclid': [111, 222],
        'right_ascension_euclid': [150.0, 151.0],
        'declination_euclid': [2.0, 2.1],
    })

    arch = EuclidArchive(environment='REG')
    arch.euclid = Mock()

    # Patch _crossmatch_batch to avoid ADQL and just echo input
    with patch.object(arch, '_crossmatch_batch', return_value=Table({'object_id': [111, 222]})) as mock_batch:
        out = arch.crossmatch_sources(
            user_table=user_table,
            radius=1.0,
            ra_col='ra',  # intentionally missing
            dec_col='dec',  # intentionally missing
            use_object_id=True,
        )

        assert len(out) == 2
        # Verify _crossmatch_batch was called with use_object_id=True
        assert mock_batch.called
        args, _ = mock_batch.call_args
        assert args[5] is True


def test_crossmatch_sources_idr_field_selection(monkeypatch):
    """Ensure IDR field selection picks the correct MER catalogue."""
    user_table = Table({
        'ra': [150.0],
        'dec': [2.0],
    })

    arch = EuclidArchive(environment='IDR')
    arch.euclid = Mock()
    arch._logged_in = True

    with patch.object(arch, '_crossmatch_batch', return_value=Table({'object_id': [123]})) as mock_batch:
        arch.crossmatch_sources(user_table=user_table, radius=1.0)
        arch.crossmatch_sources(user_table=user_table, radius=1.0, idr_field='DEEP')

    # IDR WIDE queries survey+mode; IDR DEEP queries a single MER table.
    assert mock_batch.call_count == 3

    first_call_args, _ = mock_batch.call_args_list[0]
    second_call_args, _ = mock_batch.call_args_list[1]
    third_call_args, _ = mock_batch.call_args_list[2]

    assert first_call_args[4] == 'catalogue.mer_catalogue_wide_survey'
    assert second_call_args[4] == 'catalogue.mer_catalogue_wide_mode'
    assert third_call_args[4] == 'catalogue.mer_catalogue_deep'


def test_crossmatch_sources_idr_wide_deduplicates_primary_mer_table():
    """IDR WIDE should keep wide_survey rows when wide_mode duplicates object_id."""
    user_table = Table({'source_id': [123], 'ra': [150.0], 'dec': [2.0]})
    arch = EuclidArchive(environment='IDR')
    arch.euclid = Mock()
    arch._logged_in = True

    def fake_batch(batch, ra_col, dec_col, radius, mer_table, use_object_id, force_async=False):
        if mer_table.endswith('_survey'):
            return Table({'object_id': [123], 'mer_ra': [150.0], 'source_table': ['survey']})
        return Table({'object_id': [123], 'mer_ra': [151.0], 'source_table': ['mode']})

    with patch.object(arch, '_crossmatch_batch', side_effect=fake_batch):
        out = arch.crossmatch_sources(
            user_table=user_table,
            radius=1.0,
            use_object_id=True,
            idr_field='WIDE',
        )

    assert len(out) == 1
    assert out['source_table'][0] == 'survey'


def test_spatial_crossmatch_contains_expression(monkeypatch):
    """Spatial joins should use CONTAINS predicate and compute separations locally."""
    user_table = Table({
        'RACAT': [150.0],
        'DECCAT': [2.0],
    })

    arch = EuclidArchive(environment='PDR')
    arch.euclid = Mock()
    arch._logged_in = True

    captured = {}

    def fake_launch_job(query, upload_resource=None, upload_table_name=None):
        captured['query'] = query
        mock_job = Mock()
        mock_job.get_results.return_value = Table({
            'racat': [150.0],
            'deccat': [2.0],
            'mer_ra': [150.0001],
            'mer_dec': [2.0001],
            'object_id': [123],
        })
        return mock_job

    arch.euclid.launch_job.side_effect = fake_launch_job
    arch.euclid.launch_job_async.side_effect = fake_launch_job

    result = arch.crossmatch_sources(
        user_table=user_table,
        ra_col='RACAT',
        dec_col='DECCAT',
        radius=1.0,
        use_object_id=False,
    )

    query = captured['query']
    assert "DISTANCE(" in query
    assert 'separation_arcsec' in result.colnames


def test_crossmatch_sources_full_async(tmp_path):
    """full_async should submit a single async job and return metadata."""
    user_table = Table({'ra': [150.0, 151.0], 'dec': [2.0, 2.1]})
    arch = EuclidArchive(environment='PDR')
    arch.euclid = Mock()
    arch._logged_in = True

    fake_job = Mock()
    fake_job.jobid = 'job-123'
    fake_job.phase = 'QUEUED'
    fake_job.url = 'https://tap/job/job-123'
    fake_job.remote_results_location = 'https://tap/job/job-123/results'
    fake_job.get_results.return_value = Table({'object_id': [123], 'mer_ra': [150.0], 'mer_dec': [2.0]})
    submission = {
        'job': fake_job,
        'query': 'SELECT 1',
        'upload_table_name': 'user_batch_12345',
        'row_count': len(user_table),
    }

    output_path = tmp_path / 'job_info.fits'

    with patch.object(arch, '_crossmatch_batch', return_value=submission) as mock_batch:
        job_info = arch.crossmatch_sources(
            user_table=user_table,
            radius=1.0,
            full_async=True,
            use_object_id=False,
            output_file=output_path,
        )

    mock_batch.assert_called_once()
    _, kwargs = mock_batch.call_args
    assert kwargs['force_async'] is True
    assert kwargs['fetch_results'] is False
    assert job_info['job_id'] == 'job-123'
    assert job_info['results_downloaded'] is True
    assert job_info['result_row_count'] == 1

    saved = Table.read(output_path, format='fits')
    assert len(saved) == 1
    arch.euclid.remove_jobs.assert_called_once_with(['job-123'])


def test_crossmatch_sources_full_async_chunked_persists_and_merges(tmp_path):
    """Chunked full_async should save chunk files, delete jobs, and merge output."""
    user_table = Table({
        'source_id': [1, 2, 3, 4, 5],
        'ra': [150.0, 151.0, 152.0, 153.0, 154.0],
        'dec': [2.0, 2.1, 2.2, 2.3, 2.4],
    })
    arch = EuclidArchive(environment='PDR')
    arch.euclid = Mock()
    arch._logged_in = True

    def make_submission(start_value):
        fake_job = Mock()
        fake_job.jobid = f'job-{start_value}'
        fake_job.get_results.return_value = Table({
            'object_id': [start_value, start_value + 1],
            'mer_ra': [150.0, 151.0],
            'mer_dec': [2.0, 2.1],
        })
        return {
            'job': fake_job,
            'query': f'SELECT {start_value}',
            'upload_table_name': f'user_batch_{start_value}',
        }

    submissions = [
        make_submission(10),
        make_submission(20),
        {
            'job': Mock(
                jobid='job-30',
                get_results=Mock(return_value=Table({
                    'object_id': [30],
                    'mer_ra': [152.0],
                    'mer_dec': [2.2],
                })),
            ),
            'query': 'SELECT 30',
            'upload_table_name': 'user_batch_30',
        },
    ]

    output_path = tmp_path / 'merged.fits'

    with patch.object(arch, '_crossmatch_batch', side_effect=submissions) as mock_batch:
        job_info = arch.crossmatch_sources(
            user_table=user_table,
            radius=1.0,
            full_async=True,
            use_object_id=True,
            output_file=output_path,
            async_chunk_size=2,
        )

    assert mock_batch.call_count == 3
    assert output_path.exists()
    manifest_path = tmp_path / 'merged.fits.manifest.json'
    assert manifest_path.exists()
    part_files = sorted(tmp_path.glob('merged_part_*.fits'))
    assert len(part_files) == 3

    merged = Table.read(output_path, format='fits')
    assert len(merged) == 5
    assert job_info['result_row_count'] == 5
    assert job_info['chunk_count'] == 3
    assert job_info['manifest_file'] == str(manifest_path)

    manifest = json.loads(manifest_path.read_text())
    assert len(manifest['chunks']) == 3
    assert [c['status'] for c in manifest['chunks']] == ['COMPLETED', 'COMPLETED', 'COMPLETED']
    assert [c['mer_table'] for c in manifest['chunks']] == ['catalogue.mer_catalogue'] * 3

    arch.euclid.remove_jobs.assert_any_call(['job-10'])
    arch.euclid.remove_jobs.assert_any_call(['job-20'])
    arch.euclid.remove_jobs.assert_any_call(['job-30'])
