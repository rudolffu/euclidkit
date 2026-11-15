"""Morphology partition, upload, and output contracts without live TAP access."""

import json
from unittest.mock import Mock

import numpy as np
import pytest
from astropy.table import Table, MaskedColumn

from euclidkit.core.data_access import EuclidArchive
from euclidkit.core.morphology import MORPHOLOGY_COLUMNS


@pytest.fixture
def archive():
    archive = EuclidArchive(environment='IDR')
    archive._logged_in = True
    archive.euclid = Mock()
    archive._ensure_client = Mock()
    archive._delete_async_job = Mock()
    return archive


def result(ids, priorities=None):
    table = Table({'object_id': np.asarray(ids, dtype='i8')})
    for name in MORPHOLOGY_COLUMNS:
        table[name] = np.arange(len(ids), dtype='f4')
    table['partition_priority'] = priorities if priorities is not None else [0] * len(ids)
    return table


def jobs(archive, responses, asynchronous=False):
    captured = []

    def launch(query, upload_resource, upload_table_name):
        captured.append((query, Table.read(upload_resource)))
        response = responses[len(captured) - 1]
        job = Mock(jobid=f'job-{len(captured)}')
        if isinstance(response, Exception):
            job.get_results.side_effect = response
        else:
            job.get_results.return_value = response
        return job

    method = archive.euclid.launch_job_async if asynchronous else archive.euclid.launch_job
    method.side_effect = launch
    return captured


def test_wide_sync_join_precision_and_nulls(archive, tmp_path):
    big = 9007199254740993
    source = Table({'OBJECT_ID': [big, 2, big, 3], 'label': ['a', 'b', 'c', 'd'],
                    'gini': [99., 99., 99., 99.]})
    first = result([big])
    first['gini'] = MaskedColumn([1.], mask=[True], dtype='f4')
    calls = jobs(archive, [first, result([2])])
    output = tmp_path / 'morph.fits'
    matched = archive.query_morphology(source, output)
    assert len(matched) == 3
    assert list(calls[0][1]['object_id']) == [big, 2, 3]
    assert list(calls[1][1]['object_id']) == [2, 3]
    assert 'mer_morphology_wide_survey' in calls[0][0]
    assert 'mer_morphology_wide_mode' in calls[1][0]
    assert all(f'm.{name}' in calls[0][0] for name in MORPHOLOGY_COLUMNS)
    assert 'partition_priority' not in matched.colnames
    assert matched['gini'].mask.sum() == 2
    assert np.count_nonzero(Table.read(output)['object_id'] == big) == 2
    assert 'OBJECT_ID' in source.colnames and np.all(source['gini'] == 99)


@pytest.mark.parametrize('partition,suffixes', [
    ('survey', ['deep_survey']), ('mode', ['deep_mode']),
    ('both', ['deep_survey', 'deep_mode']),
])
def test_deep_selection(archive, partition, suffixes):
    calls = jobs(archive, [result([1])] * len(suffixes))
    matched = archive.query_morphology(Table({'object_id': [1]}), idr_field='DEEP',
                                       idr_deep_partition=partition)
    assert len(matched) == 1
    assert len(calls) == len(suffixes)
    for call, suffix in zip(calls, suffixes):
        assert f'catalogue.mer_morphology_{suffix}' in call[0]
        assert len(call[1]) == 1


@pytest.mark.parametrize('ids', [[None], [1.5], [np.nan], ['abc'], [True], [2**63]])
def test_invalid_ids_before_submission(archive, ids):
    with pytest.raises(ValueError, match='object_id'):
        archive.query_morphology(Table({'object_id': ids}))
    archive.euclid.launch_job.assert_not_called()


def test_missing_id_and_environment(archive):
    with pytest.raises(ValueError, match='object_id'):
        archive.query_morphology(Table({'ra': [1]}))
    archive.environment = 'PDR'
    with pytest.raises(ValueError, match='IDR'):
        archive.query_morphology(Table({'object_id': [1]}))


def test_empty_schema(archive, tmp_path):
    jobs(archive, [Table(), Table()])
    output = tmp_path / 'empty.fits'
    matched = archive.query_morphology(Table({'object_id': [1], 'ra': [2.]}), output)
    assert len(matched) == 0
    assert matched.colnames == ['object_id', 'ra', *MORPHOLOGY_COLUMNS]
    assert Table.read(output).colnames == matched.colnames


@pytest.mark.parametrize('field,partition', [('WIDE', 'survey'), ('DEEP', 'both'), ('DEEP', 'survey')])
def test_async_query_and_precedence(archive, tmp_path, field, partition):
    raw = result([1, 1], [1, 0])
    calls = jobs(archive, [raw], asynchronous=True)
    output = tmp_path / 'async.fits'
    metadata = archive.query_morphology(Table({'object_id': [1]}), output,
        idr_field=field, idr_deep_partition=partition, full_async=True)
    query = calls[0][0]
    assert ('NOT EXISTS' in query) == (field == 'WIDE')
    assert ('UNION ALL' in query) == (field == 'WIDE' or partition == 'both')
    assert 'ORDER BY' not in query
    assert metadata['results_downloaded']
    assert Table.read(output)['gini'][0] == 1  # Survey row arrived second.
    archive._delete_async_job.assert_called_once_with('job-1')


@pytest.mark.parametrize('failure', [False, True])
def test_chunked_and_failed_downloads(archive, tmp_path, failure):
    jobs(archive, [result([1]), RuntimeError('offline') if failure else result([2])], True)
    output = tmp_path / 'chunks.fits'
    metadata = archive.query_morphology(Table({'object_id': [1, 2, 1]}), output,
                                       full_async=True, async_chunk_size=1)
    assert metadata['results_downloaded'] == (not failure)
    assert output.exists() == (not failure)
    manifest = json.loads((tmp_path / 'chunks.fits.manifest.json').read_text())
    assert len(manifest['chunks']) == 2
    assert 'partition_priority' in Table.read(tmp_path / 'chunks_part_0001.fits').colnames
    if not failure:
        assert len(Table.read(output)) == 3


def test_single_download_failure(archive, tmp_path):
    jobs(archive, [RuntimeError('offline')], True)
    output = tmp_path / 'single.fits'
    metadata = archive.query_morphology(Table({'object_id': [1]}), output, full_async=True)
    assert not metadata['results_downloaded'] and not output.exists()
    assert json.loads((tmp_path / 'single.fits.job.json').read_text())['job_id'] == 'job-1'
    archive._delete_async_job.assert_not_called()


def test_sync_batch_size(archive):
    calls = jobs(archive, [result(range(1000)), result([1000])])
    matched = archive.query_morphology(Table({'object_id': np.arange(1001)}))
    assert len(matched) == 1001
    assert [len(call[1]) for call in calls] == [1000, 1]
