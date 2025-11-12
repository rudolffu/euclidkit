"""Ensure object-id mode does not require RA/Dec columns."""

import json
from unittest.mock import Mock, patch

from astropy.table import Table

from euclidqso.core.data_access import EuclidArchive


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
        _, kwargs = mock_batch.call_args
        assert kwargs['use_object_id'] is True


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

    # Two batches processed (one per invocation)
    assert mock_batch.call_count == 2

    first_call_args, _ = mock_batch.call_args_list[0]
    second_call_args, _ = mock_batch.call_args_list[1]

    assert first_call_args[4] == 'catalogue.mer_catalogue_wide'
    assert second_call_args[4] == 'catalogue.mer_catalogue_deep'


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
            'RACAT': [150.0],
            'DECCAT': [2.0],
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
    assert "CONTAINS(" in query
    assert "DISTANCE" not in query
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
    submission = {
        'job': fake_job,
        'query': 'SELECT 1',
        'upload_table_name': 'user_batch_12345',
        'row_count': len(user_table),
    }

    output_path = tmp_path / 'job_info.json'

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

    saved = json.loads(output_path.read_text())
    assert saved['job_id'] == 'job-123'
    assert saved['query'] == 'SELECT 1'
