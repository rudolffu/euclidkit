"""Tests for object_id-based join path in _crossmatch_batch."""

from unittest.mock import Mock

import numpy as np
from astropy.table import MaskedColumn, Table

from euclidkit.core.data_access import EuclidArchive


def test_crossmatch_batch_object_id_join_query_build(monkeypatch):
    archive = EuclidArchive(environment='REG')
    archive.euclid = Mock()

    # Prepare user batch with object_id and a column colliding with MER (e.g., mu_max)
    batch = Table({
        'object_id': [1234567890],
        'ra': [150.0],
        'dec': [2.0],
        'mu_max': [25.0],
    })

    # Mock job result (minimal) to allow method to complete
    result_table = Table({
        'object_id': [1234567890],
        'mer_mu_max': [26.0],
    })
    mock_job = Mock()
    mock_job.get_results.return_value = result_table

    captured = {}

    def fake_launch_job(query, upload_resource=None, upload_table_name=None):
        captured['query'] = query
        return mock_job

    archive.euclid.launch_job.side_effect = fake_launch_job
    archive.euclid.launch_job_async.side_effect = fake_launch_job

    out = archive._crossmatch_batch(
        batch=batch,
        ra_col='ra',
        dec_col='dec',
        radius=1.0,
        mer_table='catalogue.mer_catalogue'
    )

    # Validate that the query uses equality join and avoids u.*
    q = captured['query']
    assert 'JOIN catalogue.mer_catalogue AS m ON u.object_id = m.object_id' in q
    assert 'u.*' not in q
    # User object_id is aliased to avoid collision
    assert 'u.object_id AS object_id_user' in q
    # MER columns that collide are prefixed
    assert 'm.mu_max AS mer_mu_max' in q
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
        assert f'm.{field} AS {field}' in q

    # Output is whatever job returned (no separation in this path)
    assert 'separation_arcsec' not in out.colnames
    assert 'object_id' in out.colnames


def test_sanitize_upload_table_columns_makes_masked_object_column_votable_safe(tmp_path):
    archive = EuclidArchive(environment='REG')
    table = Table(
        {
            'object_id': [1, 2],
            'addon_sampling_group': MaskedColumn(
                data=np.array([None, None], dtype=object),
                mask=[True, True],
            ),
        }
    )

    sanitized, _ = archive._sanitize_upload_table_columns(table)

    assert sanitized['addon_sampling_group'].dtype.kind in {'U', 'S'}
    sanitized.write(tmp_path / 'upload.vot', format='votable')
