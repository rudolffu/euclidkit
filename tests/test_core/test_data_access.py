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
