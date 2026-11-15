"""IDR MER morphology queries and input-table enrichment."""

import logging
import re
import tempfile
from pathlib import Path
from uuid import uuid4

import numpy as np
from astropy.table import Table, join, vstack

from euclidkit.utils.io import save_table

logger = logging.getLogger(__name__)

MORPHOLOGY_COLUMNS = tuple(
    column
    for measurement in (
        'concentration', 'asymmetry', 'smoothness', 'moment_20', 'gini',
        'sersic_sersic_vis_radius', 'sersic_sersic_vis_axis_ratio',
        'sersic_sersic_vis_index', 'sersic_angle', 'flux_radius',
    )
    for column in (measurement, measurement + '_err')
)


def _prepare_input(source_table):
    names = [name for name in source_table.colnames if name.lower() == 'object_id']
    if len(names) != 1:
        raise ValueError('Input must contain exactly one object_id column (case-insensitive).')
    result = source_table.copy()
    if names[0] != 'object_id':
        result.rename_column(names[0], 'object_id')
    ids = []
    for value in result['object_id']:
        # Never route identifiers through floating point, which can round int64 IDs.
        if (np.ma.is_masked(value) or isinstance(value, (bool, np.bool_))
                or not isinstance(value, (int, np.integer, str))
                or not re.fullmatch(r'[+-]?\d+', str(value).strip())):
            raise ValueError(f'Invalid object_id: {value!r}; expected a non-null int64 ID.')
        integer = int(value)
        if not np.iinfo(np.int64).min <= integer <= np.iinfo(np.int64).max:
            raise ValueError(f'object_id outside int64 range: {value!r}')
        ids.append(integer)
    result['object_id'] = np.asarray(ids, dtype=np.int64)
    return result, Table({'object_id': np.asarray(list(dict.fromkeys(ids)), dtype=np.int64)})


def _empty_result():
    return Table(names=['object_id', *MORPHOLOGY_COLUMNS, 'partition_priority'],
                 dtype=['i8', *(['f4'] * len(MORPHOLOGY_COLUMNS)), 'i4'])


def _query(upload, tables, wide_fallback):
    branches = []
    for priority, table in enumerate(tables):
        fields = ', '.join(['u.object_id', *[f'm.{c}' for c in MORPHOLOGY_COLUMNS],
                            f'{priority} AS partition_priority'])
        branch = (f'SELECT {fields} FROM TAP_UPLOAD.{upload} AS u '
                  f'JOIN {table} AS m ON u.object_id = m.object_id')
        if priority and wide_fallback:
            branch += (f' WHERE NOT EXISTS (SELECT 1 FROM {tables[0]} AS s '
                       'WHERE s.object_id = u.object_id)')
        branches.append(branch)
    return '\nUNION ALL\n'.join(branches)


def _submit(archive, batch, tables, wide_fallback, asynchronous):
    upload = 'user_morph_' + uuid4().hex
    query = _query(upload, tables, wide_fallback)
    with tempfile.TemporaryDirectory(prefix='euclidkit_morph_') as directory:
        path = Path(directory) / 'upload.vot'
        batch.write(path, format='votable')
        launch = archive.euclid.launch_job_async if asynchronous else archive.euclid.launch_job
        job = launch(query, upload_resource=str(path), upload_table_name=upload)
    return job, query


def _join_input(source, results):
    left = source.copy()
    for name in MORPHOLOGY_COLUMNS:
        if name in left.colnames:
            left.remove_column(name)
    right = results.copy()
    if len(right):
        right.sort('partition_priority')
        _, first = np.unique(right['object_id'], return_index=True)
        right = right[np.sort(first)]
    right = right[['object_id', *MORPHOLOGY_COLUMNS]]
    if not len(right) or not len(left):
        empty = left[:0].copy()
        for name in MORPHOLOGY_COLUMNS:
            empty.add_column(right[name][:0])
        return empty
    return join(left, right, keys='object_id', join_type='inner', metadata_conflicts='silent')


def query_morphology(archive, source_table, output_file=None, idr_field='WIDE',
                     idr_deep_partition='survey', full_async=False,
                     async_chunk_size=500000):
    """Query unique IDs and inner join compact morphology measurements to input rows."""
    if archive.environment != 'IDR':
        raise ValueError('MER morphology queries currently support only the IDR environment.')
    source, upload = _prepare_input(source_table)
    mer_tables = archive._get_mer_table_names(idr_field, idr_deep_partition)
    tables = [name.replace('mer_catalogue', 'mer_morphology') for name in mer_tables]
    wide = archive._uses_idr_wide_fallback(mer_tables)
    if full_async and (output_file is None or not isinstance(async_chunk_size, int)
                       or isinstance(async_chunk_size, bool) or async_chunk_size <= 0):
        raise ValueError('Full async requires output_file and a positive integer async_chunk_size.')
    if len(upload):
        if not archive._logged_in:
            archive.login()
        archive._ensure_client()
    logger.info('Querying morphology for %d unique IDs using %s', len(upload), tables)

    if not full_async:
        results = []
        for start in range(0, len(upload), 1000):
            batch = upload[start:start + 1000]
            for priority, table in enumerate(tables):
                if not len(batch):
                    break
                job, _ = _submit(archive, batch, [table], False, False)
                result = job.get_results()
                if result is not None and len(result):
                    result['partition_priority'] = priority
                    results.append(result)
                    if wide:
                        batch = batch[~np.isin(batch['object_id'], result['object_id'])]
        raw = vstack(results, metadata_conflicts='silent') if results else _empty_result()
        final = _join_input(source, raw)
        if output_file is not None:
            save_table(final, output_file)
        return final

    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    chunked = len(upload) > async_chunk_size
    metadata = {
        'type': 'euclidkit_morph_async_chunked' if chunked else 'euclidkit_morph_async_job',
        'environment': 'IDR', 'tables': tables, 'row_count': len(upload),
        'output_file': str(output), 'results_downloaded': False,
        'result_row_count': 0, 'chunk_size': async_chunk_size, 'chunks': [],
    }
    metadata_path = output.with_name(output.name + ('.manifest.json' if chunked else '.job.json'))
    metadata['manifest_file' if chunked else 'job_info_file'] = str(metadata_path)
    results = []
    failed = False
    for index, start in enumerate(range(0, len(upload), async_chunk_size), 1):
        job, query = _submit(archive, upload[start:start + async_chunk_size], tables, wide, True)
        job_id = getattr(job, 'jobid', None)
        chunk = {'index': index, 'job_id': job_id, 'query': query,
                 'tables': tables, 'status': 'SUBMITTED'}
        metadata['chunks'].append(chunk)
        metadata['query_count'] = len(metadata['chunks'])
        if not chunked:
            metadata.update(job_id=job_id, query=query)
        # Persist recovery information before attempting result retrieval.
        archive._write_job_info(metadata, metadata_path)
        try:
            result = job.get_results()
            result = result if result is not None and len(result) else _empty_result()
            if chunked:
                part = output.with_name(f'{output.stem}_part_{index:04d}{output.suffix or ".fits"}')
                save_table(result, part)
                chunk['file'] = str(part)
            results.append(result)
            chunk.update(status='COMPLETED', rows=len(result))
            if chunked:
                archive._delete_async_job(job_id)
        except Exception as exc:
            failed = True
            chunk.update(status='FAILED', error=str(exc))
            metadata['download_error'] = str(exc)
        archive._write_job_info(metadata, metadata_path)
    metadata['chunk_count'] = len(metadata['chunks'])
    metadata['job_ids'] = [c['job_id'] for c in metadata['chunks'] if c['job_id']]
    if not failed:
        raw = vstack(results, metadata_conflicts='silent') if results else _empty_result()
        final = _join_input(source, raw)
        save_table(final, output)
        metadata.update(results_downloaded=True, result_row_count=len(final))
        if not chunked and metadata.get('job_id'):
            archive._delete_async_job(metadata['job_id'])
    archive._write_job_info(metadata, metadata_path)
    return metadata
