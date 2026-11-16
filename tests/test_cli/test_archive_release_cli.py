"""Release-aware CLI defaults and legacy field aliases."""
from unittest.mock import patch

import pytest
from astropy.table import Table
from click.testing import CliRunner

from euclidkit.cli.main import main


@pytest.mark.parametrize('field,partition', [('--field', '--deep-partition'),
                                           ('--idr-field', '--idr-deep-partition')])
@pytest.mark.parametrize('environment,release', [('PDR', 'DR1'), ('IDR', None)])
def test_morphology_routes_and_aliases(tmp_path, field, partition, environment, release):
    source = tmp_path/'sources.fits'
    Table({'object_id': [1]}).write(source)
    args = ['query-morph', '-i', str(source), '-o', str(tmp_path/'out.fits'),
            '--environment', environment, field, 'DEEP', partition, 'both']
    if release:
        args += ['--release', release]
    with patch('euclidkit.core.data_access.EuclidArchive') as factory:
        factory.return_value.query_morphology.return_value = Table({'object_id': [1]})
        response = CliRunner().invoke(main, args)
    assert response.exit_code == 0, response.output
    factory.assert_called_once_with(environment=environment, data_release=release)
    assert factory.return_value.query_morphology.call_args.kwargs['idr_deep_partition'] == 'both'


@pytest.mark.parametrize('command', ['query-morph', 'query-zspe'])
def test_q1_rejects_unavailable_commands_before_login(tmp_path, command):
    source = tmp_path/'sources.fits'
    Table({'object_id': [1]}).write(source)
    flag = '-i' if command == 'query-morph' else '--crossmatch'
    with patch('euclidkit.core.data_access.EuclidArchive') as factory:
        response = CliRunner().invoke(main, [command, flag, str(source), '-o',
            str(tmp_path/'out.fits'), '--release', 'Q1'])
    assert response.exit_code != 0 and 'Q1 is not supported' in response.output
    factory.assert_not_called()


def test_idr_public_release_rejected_before_login(tmp_path):
    source = tmp_path/'sources.fits'
    Table({'object_id': [1]}).write(source)
    with patch('euclidkit.core.data_access.EuclidArchive') as factory:
        response = CliRunner().invoke(main, ['query-morph', '-i', str(source), '-o',
            str(tmp_path/'out.fits'), '--environment', 'IDR', '--release', 'DR1'])
    assert response.exit_code != 0 and 'require environment=PDR' in response.output
    factory.assert_not_called()


@pytest.mark.parametrize('release', [None, 'Q1'])
def test_spectra_selection_forwarded(tmp_path, release):
    source = tmp_path/'sources.fits'
    Table({'object_id': [1]}).write(source)
    args = ['query-spectra', '--crossmatch', str(source), '-o', str(tmp_path/'out.fits')]
    if release:
        args += ['--release', release]
    with patch('euclidkit.core.data_access.EuclidArchive') as factory:
        factory.return_value.query_spectra_sources.return_value = Table({'source_id': [1]})
        response = CliRunner().invoke(main, args)
    assert response.exit_code == 0, response.output
    factory.assert_called_once_with(environment='PDR', data_release=release)
    field = factory.return_value.query_spectra_sources.call_args.kwargs['idr_field']
    assert field == ('WIDE' if release is None else None)


def test_public_q1_deep_selection_rejected(tmp_path):
    source = tmp_path/'sources.fits'
    Table({'object_id': [1]}).write(source)
    response = CliRunner().invoke(main, ['query-spectra', '--crossmatch', str(source), '-o',
        str(tmp_path/'out.fits'), '--release', 'Q1', '--field', 'DEEP'])
    assert response.exit_code != 0 and 'not available for Q1' in response.output


@pytest.mark.parametrize('command', ['crossmatch', 'query-spectra', 'query-cutana', 'query-segmap',
    'query-morph', 'query-zspe', 'cutouts', 'compile-spectra', 'dithers-to-parquet',
    'download-cutana-alternatives', 'upload-table'])
def test_archive_commands_expose_release(command):
    response = CliRunner().invoke(main, [command, '--help'])
    assert response.exit_code == 0, response.output
    assert '--release [DR1|Q1]' in response.output
