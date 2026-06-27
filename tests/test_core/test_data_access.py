"""Test data access functionality."""

import pytest
import numpy as np
import pandas as pd
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path
import tempfile
import os

from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord

from euclidkit.core.data_access import EuclidArchive


class TestEuclidArchive:
    """Test cases for EuclidArchive class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.archive = EuclidArchive(environment='REG')
        
        # Mock the Euclid client
        self.mock_client = Mock()
        self.archive.euclid = self.mock_client
        
        # Sample test data
        self.sample_coordinates = Table({
            'ra': [150.0, 151.0, 152.0],
            'dec': [2.0, 2.1, 2.2],
            'source_id': [1, 2, 3]
        })
        
        self.sample_crossmatch_result = Table({
            'object_id': [100001, 100002],
            'ra': [150.001, 151.001],
            'dec': [2.001, 2.101],
            'separation_arcsec': [0.1, 0.2],
            'user_source_id': [1, 2]
        })

    def test_init_default_environment(self):
        """Test default initialization."""
        archive = EuclidArchive()
        assert archive.environment == 'PDR'
        assert archive.verbose is False
        
    def test_init_custom_environment(self):
        """Test initialization with custom environment."""
        archive = EuclidArchive(environment='OTF', verbose=True)
        assert archive.environment == 'OTF'
        assert archive.verbose is True

    def test_get_mer_table_name(self):
        """Test MER table name generation for different environments."""
        test_cases = [
            ('PDR', 'catalogue.mer_catalogue'),
            ('IDR', 'catalogue.mer_catalogue'),
            ('OTF', 'catalogue.mer_catalogue'),
            ('REG', 'test.mer_catalogue')
        ]
        
        for env, expected in test_cases:
            archive = EuclidArchive(environment=env)
            assert archive.get_mer_table_name() == expected

    def test_get_spectra_source_table_name(self):
        """Spectra source table should map to the expected release prefix by environment."""
        test_cases = [
            ('PDR', 'q1.spectra_source'),
            ('OTF', 'q1.spectra_source'),
            ('REG', 'dr1.spectra_source'),
            ('IDR', 'dr1.spectra_source'),
        ]

        for env, expected in test_cases:
            archive = EuclidArchive(environment=env)
            assert archive._get_spectra_source_table_name() == expected

    def test_get_segmentation_map_table_name(self):
        """Segmentation-map table should map to the expected release by environment."""
        test_cases = [
            ('PDR', 'q1.mer_segmentation_map'),
            ('IDR', 'dr1.mer_segmentation_map'),
            ('OTF', 'sedm.mer_segmentation_map'),
            ('REG', 'sedm.mer_segmentation_map'),
        ]

        for env, expected in test_cases:
            archive = EuclidArchive(environment=env)
            assert archive._get_segmentation_map_table_name() == expected

    def test_query_segmentation_maps_requires_segmentation_id(self):
        """query-segmap should ask users to crossmatch when SEGMENTATION_MAP_ID is absent."""
        archive = EuclidArchive(environment='PDR')
        archive._logged_in = True

        with pytest.raises(ValueError, match="Run euclidkit crossmatch first"):
            archive.query_segmentation_maps(Table({
                'object_id': [1],
                'ra': [150.0],
                'dec': [2.0],
            }))

    def test_query_segmentation_maps_requires_output_columns(self):
        """query-segmap should validate object_id, ra, and dec input columns."""
        archive = EuclidArchive(environment='PDR')
        archive._logged_in = True

        with pytest.raises(ValueError, match="object_id"):
            archive.query_segmentation_maps(Table({
                'SEGMENTATION_MAP_ID': [1234567890000],
                'ra': [150.0],
                'dec': [2.0],
            }))

    def test_query_segmentation_maps_builds_tile_join_and_saves(self, tmp_path):
        """query-segmap should compute tile_index and join TAP_UPLOAD to mer_segmentation_map."""
        archive = EuclidArchive(environment='IDR')
        archive._logged_in = True
        archive.euclid = Mock()
        captured = {}

        result_table = Table({
            'object_id': [100001],
            'segmentation_map_id': [1234567890000],
            'ra': [150.0],
            'dec': [2.0],
            'datalabs_path': ['/data/euclid'],
            'file_path': ['/data/euclid/mock.fits'],
            'file_name': ['mock.fits'],
            'crpix1': [1.0],
            'crpix2': [2.0],
            'crval1': [150.0],
            'crval2': [2.0],
            'seg_dec': [2.0],
            'seg_ra': [150.0],
            'data_set_release': ['dr1'],
            'environment': ['IDR'],
            'tile_index': [1234567],
            'processing_mode': ['wide'],
        })
        job = Mock()
        job.get_results.return_value = result_table

        def fake_launch_job(query, upload_resource=None, upload_table_name=None):
            captured['query'] = query
            captured['upload_table'] = Table.read(upload_resource, format='votable')
            captured['upload_table_name'] = upload_table_name
            return job

        archive.euclid.launch_job.side_effect = fake_launch_job
        source = Table({
            'OBJECT_ID': [100001, 100002, 100003, 100004],
            'SEGMENTATION_MAP_ID': [1234567890000, None, '9876543210000', 'bad'],
            'RA': [150.0, 151.0, 152.0, 153.0],
            'DEC': [2.0, 2.1, 2.2, 2.3],
        })

        with patch('euclidkit.core.data_access.save_table') as mock_save:
            result = archive.query_segmentation_maps(
                source_table=source,
                output_file=tmp_path / 'segmap.fits',
            )

        assert result is result_table
        query = captured['query']
        assert 'TAP_UPLOAD' in query
        assert 'dr1.mer_segmentation_map' in query
        assert 's.dec AS seg_dec' in query
        assert 's.ra AS seg_ra' in query
        assert 'u.tile_index = s.tile_index' in query
        assert 'SUBSTRING' not in query
        upload_table = captured['upload_table']
        assert list(upload_table['tile_index']) == [1234567, 9876543]
        assert archive._last_segmap_valid_tile_count == 2
        mock_save.assert_called_once()

    def test_query_segmentation_maps_accepts_mer_ra_dec_fallback(self):
        """query-segmap should use mer_ra/mer_dec when ra/dec are absent."""
        archive = EuclidArchive(environment='PDR')
        source = Table({
            'object_id': [100001],
            'segmentation_map_id': [1234567890000],
            'mer_ra': [150.0001],
            'mer_dec': [2.0001],
        })

        upload_table = archive._prepare_segmentation_upload_table(source)

        assert 'ra' in upload_table.colnames
        assert 'dec' in upload_table.colnames
        assert list(upload_table['ra']) == [150.0001]
        assert list(upload_table['dec']) == [2.0001]
        assert list(upload_table['tile_index']) == [1234567]

    def test_get_zspe_candidate_table_names(self):
        """SPE redshift candidate tables should resolve by object type and IDR field."""
        archive = EuclidArchive(environment='IDR')

        assert archive._get_zspe_candidate_table_names('qso', 'WIDE') == [
            ('catalogue.spectro_z_spe_qso_candidates_wide_survey', 'wide_survey'),
            ('catalogue.spectro_z_spe_qso_candidates_wide', 'wide'),
        ]
        assert archive._get_zspe_candidate_table_names('galaxy', 'DEEP') == [
            ('catalogue.spectro_z_spe_galaxy_candidates_deep', 'deep'),
        ]

    def test_query_zspe_candidates_wide_fallback(self):
        """WIDE zspe queries should pass only survey-unmatched objects to wide."""
        archive = EuclidArchive(environment='IDR')
        archive._logged_in = True
        cross_tab = Table({'object_id': [1, 2, 3], 'label': ['a', 'b', 'c']})
        survey_result = Table({
            'object_id': [1, 1, 2],
            'spe_rank': [1, 2, 1],
            'spe_z': [2.1, 2.2, 0.7],
            'spe_z_err': [0.01, 0.02, 0.03],
            'source_table': ['wide_survey', 'wide_survey', 'wide_survey'],
        })
        wide_result = Table({
            'object_id': [3],
            'spe_rank': [1],
            'spe_z': [1.1],
            'spe_z_err': [0.04],
            'source_table': ['wide'],
        })
        archive._query_zspe_batch = Mock(side_effect=[survey_result, wide_result])

        result = archive.query_zspe_candidates(
            crossmatch_table=cross_tab,
            object_type='qso',
            idr_field='WIDE',
        )

        assert len(result) == 4
        assert 'label' in result.colnames
        assert 'spe_rank' in result.colnames
        assert 'spe_z' in result.colnames
        assert 'spe_z_err' in result.colnames
        assert 'source_table' not in result.colnames
        assert archive._query_zspe_batch.call_count == 2
        first_batch = archive._query_zspe_batch.call_args_list[0].args[0]
        second_batch = archive._query_zspe_batch.call_args_list[1].args[0]
        assert list(first_batch['object_id']) == [1, 2, 3]
        assert list(second_batch['object_id']) == [3]

    def test_query_zspe_candidates_requires_object_id(self):
        """SPE redshift candidate queries require object_id input."""
        archive = EuclidArchive(environment='IDR')
        archive._logged_in = True

        with pytest.raises(ValueError, match="object_id"):
            archive.query_zspe_candidates(crossmatch_table=Table({'source_id': [1]}))

    def test_query_zspe_batch_has_no_order_by(self):
        """Synchronous SPE redshift batch query should not use server-side ordering."""
        archive = EuclidArchive(environment='IDR')
        archive.euclid = Mock()
        job = Mock()
        job.get_results.return_value = Table()
        archive.euclid.launch_job.return_value = job

        archive._query_zspe_batch(
            Table({'object_id': [1]}),
            'catalogue.spectro_z_spe_qso_candidates_deep',
            'deep',
        )

        query = archive.euclid.launch_job.call_args.args[0]
        assert 'ORDER BY' not in query

    def test_query_zspe_candidates_full_async_wide_union_query(self, tmp_path):
        """WIDE full_async should use one server-side union fallback query."""
        archive = EuclidArchive(environment='IDR')
        archive._logged_in = True
        archive.euclid = Mock()
        job = Mock()
        job.jobid = 'wide-job'
        job.get_results.return_value = Table({
            'object_id': [1, 2],
            'spe_rank': [1, 1],
            'spe_z': [2.0, 1.0],
            'spe_z_err': [0.1, 0.2],
            'source_table': ['wide_survey', 'wide'],
        })
        archive.euclid.launch_job_async.return_value = job

        meta = archive.query_zspe_candidates(
            crossmatch_table=Table({'object_id': [1, 2], 'label': ['a', 'b']}),
            output_file=tmp_path / 'zs.fits',
            object_type='qso',
            idr_field='WIDE',
            full_async=True,
        )

        query = archive.euclid.launch_job_async.call_args.args[0]
        assert 'catalogue.spectro_z_spe_qso_candidates_wide_survey' in query
        assert 'catalogue.spectro_z_spe_qso_candidates_wide' in query
        assert 'UNION ALL' in query
        assert 'NOT EXISTS' in query
        assert 'ORDER BY' not in query
        assert meta['results_downloaded'] is True
        assert meta['result_row_count'] == 2
        saved = Table.read(tmp_path / 'zs.fits', format='fits')
        assert 'label' in saved.colnames
        assert 'spe_rank' in saved.colnames
        assert 'spe_z' in saved.colnames
        assert 'spe_z_err' in saved.colnames
        assert 'source_table' not in saved.colnames

    def test_query_zspe_candidates_full_async_deep_single_table(self, tmp_path):
        """DEEP full_async should query only the deep candidate table."""
        archive = EuclidArchive(environment='IDR')
        archive._logged_in = True
        archive.euclid = Mock()
        job = Mock()
        job.jobid = 'deep-job'
        job.get_results.return_value = Table({
            'object_id': [1],
            'spe_rank': [1],
            'spe_z': [0.5],
            'spe_z_err': [0.1],
            'source_table': ['deep'],
        })
        archive.euclid.launch_job_async.return_value = job

        archive.query_zspe_candidates(
            crossmatch_table=Table({'object_id': [1], 'label': ['a']}),
            output_file=tmp_path / 'deep.fits',
            object_type='galaxy',
            idr_field='DEEP',
            full_async=True,
        )

        query = archive.euclid.launch_job_async.call_args.args[0]
        assert 'catalogue.spectro_z_spe_galaxy_candidates_deep' in query
        assert 'wide_survey' not in query
        assert 'UNION ALL' not in query
        assert 'ORDER BY' not in query
        saved = Table.read(tmp_path / 'deep.fits', format='fits')
        assert 'label' in saved.colnames
        assert 'source_table' not in saved.colnames

    def test_query_zspe_candidates_full_async_download_failure_writes_job_info(self, tmp_path):
        """Full async should preserve job metadata when immediate download fails."""
        archive = EuclidArchive(environment='IDR')
        archive._logged_in = True
        archive.euclid = Mock()
        job = Mock()
        job.jobid = 'failed-job'
        job.get_results.side_effect = RuntimeError('not ready')
        archive.euclid.launch_job_async.return_value = job
        output_path = tmp_path / 'failed.fits'

        meta = archive.query_zspe_candidates(
            crossmatch_table=Table({'object_id': [1]}),
            output_file=output_path,
            full_async=True,
        )

        assert meta['results_downloaded'] is False
        assert meta['download_error'] == 'not ready'
        assert Path(meta['job_info_file']).exists()

    def test_query_zspe_candidates_full_async_chunked_merges_manifest(self, tmp_path):
        """Chunked zspe full_async should save parts, manifest, and merged output."""
        archive = EuclidArchive(environment='IDR')
        archive._logged_in = True
        archive.euclid = Mock()
        jobs = []
        for idx in range(2):
            job = Mock()
            job.jobid = f'chunk-{idx}'
            job.get_results.return_value = Table({
                'object_id': [idx + 1],
                'spe_rank': [1],
                'spe_z': [idx + 0.5],
                'spe_z_err': [0.1],
                'source_table': ['wide_survey'],
            })
            jobs.append(job)
        archive.euclid.launch_job_async.side_effect = jobs
        output_path = tmp_path / 'chunked.fits'

        meta = archive.query_zspe_candidates(
            crossmatch_table=Table({'object_id': [1, 2, 3], 'label': ['a', 'b', 'c']}),
            output_file=output_path,
            full_async=True,
            async_chunk_size=2,
        )

        assert meta['results_downloaded'] is True
        assert meta['chunk_count'] == 2
        assert meta['result_row_count'] == 2
        assert Path(meta['manifest_file']).exists()
        assert output_path.exists()
        saved = Table.read(output_path, format='fits')
        assert 'label' in saved.colnames
        assert 'source_table' not in saved.colnames

    @patch('euclidkit.core.data_access.Euclid')
    def test_login_default_credentials(self, mock_euclid_class):
        """Test login with default credentials."""
        mock_instance = Mock()
        mock_euclid_class.return_value = mock_instance
        
        archive = EuclidArchive()
        archive.login()
        
        mock_euclid_class.assert_called_once()
        mock_instance.login.assert_called_once_with('/media/user/cred.txt')

    @patch('euclidkit.core.data_access.Euclid')
    def test_login_custom_credentials(self, mock_euclid_class):
        """Test login with custom credentials file."""
        mock_instance = Mock()
        mock_euclid_class.return_value = mock_instance
        
        archive = EuclidArchive()
        archive.login(credentials_file='/path/to/cred.txt')
        
        mock_instance.login.assert_called_once_with('/path/to/cred.txt')

    def test_logout(self):
        """Test logout functionality."""
        self.archive.logout()
        self.mock_client.logout.assert_called_once()

    def test_crossmatch_sources_astropy_table(self):
        """Test crossmatching with astropy Table input."""
        # Mock the query result
        mock_result = self.sample_crossmatch_result
        self.mock_client.query_object_async.return_value.get_results.return_value = mock_result
        
        result = self.archive.crossmatch_sources(
            user_table=self.sample_coordinates,
            radius=1.0
        )
        
        assert isinstance(result, Table)
        assert len(result) == 2
        assert 'separation_arcsec' in result.colnames

    def test_crossmatch_sources_pandas_dataframe(self):
        """Test crossmatching with pandas DataFrame input."""
        df = self.sample_coordinates.to_pandas()
        mock_result = self.sample_crossmatch_result
        self.mock_client.query_object_async.return_value.get_results.return_value = mock_result
        
        result = self.archive.crossmatch_sources(
            user_table=df,
            radius=1.0
        )
        
        assert isinstance(result, Table)
        assert len(result) == 2

    def test_crossmatch_sources_custom_columns(self):
        """Test crossmatching with custom RA/Dec column names."""
        custom_table = Table({
            'right_ascension': [150.0, 151.0],
            'declination': [2.0, 2.1],
            'id': [1, 2]
        })
        
        mock_result = self.sample_crossmatch_result
        self.mock_client.query_object_async.return_value.get_results.return_value = mock_result
        
        result = self.archive.crossmatch_sources(
            user_table=custom_table,
            radius=1.0,
            ra_col='right_ascension',
            dec_col='declination'
        )
        
        assert isinstance(result, Table)

    def test_crossmatch_sources_file_input(self):
        """Test crossmatching with file path input."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            self.sample_coordinates.write(f.name, format='csv', overwrite=True)
            
        try:
            mock_result = self.sample_crossmatch_result
            self.mock_client.query_object_async.return_value.get_results.return_value = mock_result
            
            with patch('euclidkit.core.data_access.load_table', return_value=self.sample_coordinates):
                result = self.archive.crossmatch_sources(
                    user_table=f.name,
                    radius=1.0
                )
            
            assert isinstance(result, Table)
        finally:
            os.unlink(f.name)

    def test_crossmatch_sources_output_file(self):
        """Test crossmatching with output file."""
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as f:
            output_path = f.name
            
        try:
            mock_result = self.sample_crossmatch_result
            self.mock_client.query_object_async.return_value.get_results.return_value = mock_result
            
            with patch('astropy.table.Table.write') as mock_write:
                result = self.archive.crossmatch_sources(
                    user_table=self.sample_coordinates,
                    radius=1.0,
                    output_file=output_path
                )
                
                mock_write.assert_called_once()
                assert isinstance(result, Table)
        finally:
            if os.path.exists(output_path):
                os.unlink(output_path)

    def test_crossmatch_sources_max_sources_limit(self):
        """Test crossmatching with max_sources limit."""
        large_table = Table({
            'ra': np.random.uniform(150, 152, 100),
            'dec': np.random.uniform(2, 3, 100),
            'source_id': range(100)
        })
        
        mock_result = self.sample_crossmatch_result
        self.mock_client.query_object_async.return_value.get_results.return_value = mock_result
        
        result = self.archive.crossmatch_sources(
            user_table=large_table,
            radius=1.0,
            max_sources=10
        )
        
        # Should have processed only 10 sources
        assert isinstance(result, Table)

    def test_query_spectra_sources_with_crossmatch_table(self):
        """Test spectral source querying with crossmatch table."""
        mock_spectra = Table({
            'object_id': [100001, 100002],
            'instrument_name': ['NISP', 'VIS'],
            'spectrum_id': ['spec_1', 'spec_2'],
            'file_path': ['/path/1.fits', '/path/2.fits']
        })
        
        self.mock_client.query_object_async.return_value.get_results.return_value = mock_spectra
        
        result = self.archive.query_spectra_sources(
            crossmatch_table=self.sample_crossmatch_result
        )
        
        assert isinstance(result, Table)
        assert 'spectrum_id' in result.colnames

    def test_query_spectra_sources_with_crossmatch_ids(self):
        """Test spectral source querying using object IDs from crossmatch table."""
        cross_tab = Table({
            'object_id': [100001, 100002, 100003],
        })
        mock_spectra = Table({
            'object_id': [100001, 100002, 100003],
            'instrument_name': ['NISP', 'VIS', 'NISP'],
            'spectrum_id': ['spec_1', 'spec_2', 'spec_3']
        })
        self.mock_client.launch_job.return_value.get_results.return_value = mock_spectra
        result = self.archive.query_spectra_sources(crossmatch_table=cross_tab)
        assert isinstance(result, Table)
        assert len(result) == 3

    def test_query_spectra_batch_uses_environment_table_name(self):
        """Spectra batch query should join the environment-specific spectra table."""
        batch = Table({'object_id': [100001]})
        self.archive.environment = 'IDR'
        self.mock_client.launch_job.return_value.get_results.return_value = Table({
            'source_id': [100001],
            'object_id': [100001],
        })

        self.archive._query_spectra_batch(batch)

        args, _ = self.mock_client.launch_job.call_args
        query = args[0]
        assert 'JOIN dr1.spectra_source AS s' in query

    def test_query_spectra_sources_with_output_file(self):
        """Test spectral source querying with output file."""
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as f:
            output_path = f.name
            
        try:
            mock_spectra = Table({'object_id': [100001], 'spectrum_id': ['spec_1']})
            self.mock_client.query_object_async.return_value.get_results.return_value = mock_spectra
            
            with patch('astropy.table.Table.write') as mock_write:
                cross_tab = Table({'object_id': [100001]})
                result = self.archive.query_spectra_sources(
                    crossmatch_table=cross_tab,
                    output_file=output_path
                )
                
                mock_write.assert_called_once()
        finally:
            if os.path.exists(output_path):
                os.unlink(output_path)

    def test_query_spectra_sources_error_no_input(self):
        """Test error when no input provided to query_spectra_sources."""
        with pytest.raises(ValueError, match="Must provide crossmatch_table"):
            self.archive.query_spectra_sources()

    def test_batch_processing(self):
        """Test batch processing for large datasets."""
        # Test the _batch_process method
        large_data = list(range(2500))  # > batch_size of 1000
        
        def mock_process_func(batch):
            return [x * 2 for x in batch]
        
        results = self.archive._batch_process(large_data, mock_process_func, batch_size=1000)
        
        assert len(results) == 2500
        assert results[0] == 0
        assert results[-1] == 4998

    def test_crossmatch_sources_forces_async_batches_over_2000_rows(self):
        """Crossmatch should force async batch jobs when input row count exceeds 2000."""
        self.archive._logged_in = True
        large_table = Table({
            'ra': np.linspace(150.0, 151.0, 2501),
            'dec': np.linspace(2.0, 3.0, 2501),
        })

        with patch.object(self.archive, '_crossmatch_batch', return_value=Table()) as mock_batch:
            self.archive.crossmatch_sources(user_table=large_table, radius=1.0)

        assert mock_batch.call_count == 3
        for call in mock_batch.call_args_list:
            assert call.kwargs['force_async'] is True

    def test_upload_user_table_from_file(self, tmp_path):
        """Uploading from file should infer format and allow overwrite."""
        table_path = tmp_path / 'sample.csv'
        table_path.write_text('ra,dec\n10.0,2.0\n')

        fake_job = Mock()
        fake_job.jobid = 'job-42'
        fake_job.phase = 'EXECUTING'
        self.mock_client.upload_table.return_value = fake_job

        result = self.archive.upload_user_table(
            table=table_path,
            table_name='user_table',
            description='demo table',
            overwrite=True,
            verbose=True,
        )

        self.mock_client.delete_user_table.assert_called_once_with(
            table_name='user_table',
            force_removal=True,
            verbose=True,
        )
        self.mock_client.upload_table.assert_called_once_with(
            upload_resource=str(table_path),
            table_name='user_table',
            table_description='demo table',
            format='csv',
            verbose=True,
        )
        assert result['job_id'] == 'job-42'
        assert result['format'] == 'csv'

    def test_upload_user_table_astropy_table(self):
        """Uploading an astropy Table defaults to VOTable format."""
        tbl = Table({'ra': [10.0], 'dec': [2.0]})
        self.mock_client.upload_table.return_value = None

        result = self.archive.upload_user_table(
            table=tbl,
            table_name='table_from_astropy',
            description=None,
        )

        args, kwargs = self.mock_client.upload_table.call_args
        assert kwargs['upload_resource'] is tbl
        assert kwargs['format'] == 'votable'
        assert result['resource_type'] == 'table'
