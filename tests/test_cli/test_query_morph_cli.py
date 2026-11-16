from unittest.mock import patch

import pytest
from astropy.table import Table
from click.testing import CliRunner

from euclidkit.cli.main import main


@pytest.mark.parametrize('asynchronous', [False, True])
def test_query_morph_forwarding(tmp_path, asynchronous):
    source = tmp_path / 'input.fits'
    Table({'OBJECT_ID': [1]}).write(source)
    args = ['query-morph', '-i', str(source), '-o', str(tmp_path / 'out.fits'),
            '--idr-field', 'DEEP', '--idr-deep-partition', 'both']
    if asynchronous:
        args += ['--full-async', '--async-chunk-size', '10']
    with patch('euclidkit.core.data_access.EuclidArchive') as factory:
        archive = factory.return_value
        archive.query_morphology.return_value = (
            {'results_downloaded': True, 'result_row_count': 1, 'chunk_count': 1,
             'job_ids': ['test-job']} if asynchronous else Table({'object_id': [1]})
        )
        response = CliRunner().invoke(main, args)
        assert response.exit_code == 0, response.output
        assert 'Matched rows: 1' in response.output
        kwargs = archive.query_morphology.call_args.kwargs
        assert kwargs['full_async'] == asynchronous
        assert kwargs['async_chunk_size'] == (10 if asynchronous else 500000)
        assert kwargs['idr_deep_partition'] == 'both'
        factory.assert_called_once_with(environment='PDR', data_release=None)
        archive.login.assert_called_once()
        archive.logout.assert_called_once()


def test_query_morph_help():
    response = CliRunner().invoke(main, ['--help'])
    assert 'Query MER morphology by object ID.' in response.output


def test_query_morph_error(tmp_path):
    source = tmp_path / 'input.fits'
    Table({'ra': [1]}).write(source)
    with patch('euclidkit.core.data_access.EuclidArchive') as factory:
        factory.return_value.query_morphology.side_effect = ValueError('Missing object_id')
        response = CliRunner().invoke(main, ['query-morph', '-i', str(source), '-o', 'out.fits'])
        assert response.exit_code == 1 and 'Missing object_id' in response.output
        factory.return_value.logout.assert_called_once()
