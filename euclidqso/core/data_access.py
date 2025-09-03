"""
Data access module for euclidqso package.

Provides interfaces to Euclid science archive and data volumes.
"""

import os
import logging
from pathlib import Path
from typing import Union, Optional, List, Dict, Any, Tuple
import tempfile

import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.esa.euclid.core import EuclidClass

from euclidqso.config import config
from euclidqso.utils.io import load_table, save_table

logger = logging.getLogger(__name__)


class EuclidArchive:
    """Main interface to Euclid science archive."""
    
    def __init__(self, environment: str = 'PDR'):
        """
        Initialize Euclid archive client.
        
        Parameters
        ----------
        environment : str
            Archive environment: 'PDR', 'OTF', 'REG', or 'IDR'
        """
        self.environment = environment
        self.euclid = EuclidClass(environment=environment)
        self._logged_in = False
        
        logger.info(f"Initialized EuclidArchive for {environment} environment")
    
    def login(self, credentials_file: Optional[str] = None, user: Optional[str] = None):
        """
        Login to Euclid archive.
        
        Parameters
        ----------
        credentials_file : str, optional
            Path to credentials file
        user : str, optional
            Username for interactive login
        """
        if credentials_file is None:
            credentials_file = config.get('data.credentials_file')
        
        try:
            if credentials_file and Path(credentials_file).exists():
                self.euclid.login(credentials_file=credentials_file)
                logger.info("Successfully logged in with credentials file")
            elif user:
                self.euclid.login(user=user)
                logger.info(f"Successfully logged in as {user}")
            else:
                logger.warning("No credentials provided - some functionality may be limited")
                return
            
            self._logged_in = True
            
        except Exception as e:
            logger.error(f"Login failed: {e}")
            raise
    
    def logout(self):
        """Logout from Euclid archive."""
        if self._logged_in:
            self.euclid.logout()
            self._logged_in = False
            logger.info("Logged out from Euclid archive")
    
    def crossmatch_sources(
        self,
        user_table: Union[str, Path, Table, pd.DataFrame],
        radius: float = 1.0,
        mer_table: str = None,
        output_file: Optional[Union[str, Path]] = None,
        ra_col: str = 'ra',
        dec_col: str = 'dec',
        max_sources: Optional[int] = None
    ) -> Table:
        """
        Crossmatch user table with Euclid MER catalogue.
        
        Parameters
        ----------
        user_table : str, Path, Table, or DataFrame
            Input table with source positions
        radius : float, default 1.0
            Search radius in arcseconds
        mer_table : str, optional
            MER catalogue table name (environment-dependent)
        output_file : str or Path, optional
            Output file path to save results
        ra_col : str, default 'ra'
            RA column name in user table
        dec_col : str, default 'dec'
            Dec column name in user table
        max_sources : int, optional
            Maximum number of sources to process
            
        Returns
        -------
        Table
            Crossmatched results
        """
        if not self._logged_in:
            logger.warning("Not logged in - attempting login with default credentials")
            self.login()
        
        # Load user table
        if isinstance(user_table, (str, Path)):
            user_data = load_table(user_table)
            logger.info(f"Loaded user table from {user_table}")
        elif isinstance(user_table, pd.DataFrame):
            user_data = Table.from_pandas(user_table)
            logger.info("Converted pandas DataFrame to astropy Table")
        elif isinstance(user_table, Table):
            user_data = user_table.copy()
            logger.info("Using provided astropy Table")
        else:
            raise ValueError("user_table must be file path, astropy Table, or pandas DataFrame")
        
        if max_sources:
            user_data = user_data[:max_sources]
            logger.info(f"Limited to {max_sources} sources")
        
        # Determine MER table name based on environment
        if mer_table is None:
            mer_table = self._get_mer_table_name()
        
        logger.info(f"Using MER table: {mer_table}")
        logger.info(f"Crossmatching {len(user_data)} sources with radius {radius}\"")
        
        # Check column names
        if ra_col not in user_data.colnames:
            # Try alternative column names
            for alt_ra in ['right_ascension', 'RA', 'ra_deg']:
                if alt_ra in user_data.colnames:
                    ra_col = alt_ra
                    break
            else:
                raise ValueError(f"RA column '{ra_col}' not found in user table")
        
        if dec_col not in user_data.colnames:
            # Try alternative column names
            for alt_dec in ['declination', 'DEC', 'dec_deg']:
                if alt_dec in user_data.colnames:
                    dec_col = alt_dec
                    break
            else:
                raise ValueError(f"Dec column '{dec_col}' not found in user table")
        
        # Perform crossmatch using ADQL
        crossmatch_results = []
        batch_size = 1000  # Process in batches
        
        for i in range(0, len(user_data), batch_size):
            batch = user_data[i:i+batch_size]
            logger.info(f"Processing batch {i//batch_size + 1}/{(len(user_data)-1)//batch_size + 1}")
            
            batch_result = self._crossmatch_batch(
                batch, ra_col, dec_col, radius, mer_table
            )
            
            if len(batch_result) > 0:
                crossmatch_results.append(batch_result)
        
        if crossmatch_results:
            final_result = Table(np.concatenate([r for r in crossmatch_results]))
        else:
            logger.warning("No crossmatches found")
            return Table()
        
        logger.info(f"Found {len(final_result)} crossmatches")
        
        # Save results if requested
        if output_file:
            save_table(final_result, output_file)
            logger.info(f"Saved crossmatch results to {output_file}")
        
        return final_result
    
    def _get_mer_table_name(self) -> str:
        """Get MER catalogue table name for current environment."""
        table_names = {
            'PDR': 'catalogue.mer_catalogue',
            'IDR': 'catalogue.mer_catalogue', 
            'OTF': 'catalogue.mer_catalogue',
            'REG': 'catalogue.mer_final_catalog_fits_file_regreproc1_r2'
        }
        return table_names.get(self.environment, 'catalogue.mer_catalogue')
    
    def _crossmatch_batch(
        self, 
        batch: Table, 
        ra_col: str, 
        dec_col: str, 
        radius: float, 
        mer_table: str
    ) -> Table:
        """Crossmatch a batch of sources."""
        
        # Upload user table to archive for crossmatching
        with tempfile.NamedTemporaryFile(suffix='.vot', delete=False) as tmp_file:
            batch.write(tmp_file.name, format='votable', overwrite=True)
            tmp_name = tmp_file.name
        
        try:
            # Use temporary upload with launch_job
            upload_name = f"user_batch_{np.random.randint(10000, 99999)}"
            
            # Construct crossmatch query
            radius_deg = radius / 3600.0  # Convert to degrees
            
            query = f"""
            SELECT u.*, 
                   m.object_id, 
                   m.right_ascension AS mer_ra, 
                   m.declination AS mer_dec,
                   m.vis_det, m.det_quality_flag, m.spurious_flag,
                   m.flux_detection_total, m.flux_vis_sersic, 
                   m.flux_y_sersic, m.flux_j_sersic, m.flux_h_sersic,
                   m.segmentation_map_id, m.segmentation_area, m.kron_radius,
                   DISTANCE(u.{ra_col}, u.{dec_col}, m.right_ascension, m.declination) AS separation_deg
            FROM TAP_UPLOAD.{upload_name} AS u
            JOIN {mer_table} AS m 
            ON DISTANCE(u.{ra_col}, u.{dec_col}, m.right_ascension, m.declination) < {radius_deg}
            ORDER BY u.{ra_col}, separation_deg
            """
            
            # Execute query with temporary table upload
            if len(batch) < 2000:
                job = self.euclid.launch_job(query, upload_resource=tmp_name, 
                                           upload_table_name=upload_name)
            else:
                job = self.euclid.launch_job_async(query, upload_resource=tmp_name,
                                                 upload_table_name=upload_name)
            
            result = job.get_results()
            
            # Convert separation from degrees to arcseconds
            if len(result) > 0 and 'separation_deg' in result.colnames:
                result['separation_arcsec'] = result['separation_deg'] * 3600.0
                # Remove the degree column
                result.remove_column('separation_deg')
            
            return result
            
        finally:
            # Clean up temporary file
            if os.path.exists(tmp_name):
                os.unlink(tmp_name)
    
    def query_spectra_sources(
        self, 
        object_ids: Optional[List[int]] = None,
        crossmatch_table: Optional[Table] = None,
        output_file: Optional[Union[str, Path]] = None
    ) -> Table:
        """
        Query spectral sources from Euclid archive.
        
        Parameters
        ----------
        object_ids : list of int, optional
            List of object IDs to query
        crossmatch_table : Table, optional
            Table with object_id column from crossmatch results
        output_file : str or Path, optional
            Output file path to save results
            
        Returns
        -------
        Table
            Spectral source information
        """
        if not self._logged_in:
            logger.warning("Not logged in - attempting login with default credentials")
            self.login()
        
        # Determine object IDs to query
        if object_ids is None and crossmatch_table is not None:
            if 'object_id' not in crossmatch_table.colnames:
                raise ValueError("crossmatch_table must contain 'object_id' column")
            object_ids = list(set(crossmatch_table['object_id']))
        elif object_ids is None:
            raise ValueError("Must provide either object_ids or crossmatch_table")
        
        logger.info(f"Querying spectra for {len(object_ids)} unique objects")
        
        # Query in batches to avoid query limits
        batch_size = 1000
        all_results = []
        
        for i in range(0, len(object_ids), batch_size):
            batch_ids = object_ids[i:i+batch_size]
            logger.info(f"Querying batch {i//batch_size + 1}/{(len(object_ids)-1)//batch_size + 1}")
            
            # Create IN clause for object IDs
            ids_str = ','.join(map(str, batch_ids))
            
            query = f"""
            SELECT s.source_id, s.ra_obj, s.dec_obj, s.datalabs_path, 
                   s.file_name, s.hdu_index, s.instrument_name,
                   m.object_id, m.right_ascension, m.declination
            FROM sedm.spectra_source AS s
            JOIN {self._get_mer_table_name()} AS m
            ON s.source_id = m.object_id
            WHERE m.object_id IN ({ids_str})
            AND s.datalabs_path IS NOT NULL
            ORDER BY s.source_id
            """
            
            try:
                if len(batch_ids) < 2000:
                    job = self.euclid.launch_job(query)
                else:
                    job = self.euclid.launch_job_async(query)
                
                batch_result = job.get_results()
                
                if len(batch_result) > 0:
                    all_results.append(batch_result)
                    logger.info(f"Found {len(batch_result)} spectra in batch")
                
            except Exception as e:
                logger.error(f"Error querying batch: {e}")
                continue
        
        if all_results:
            final_result = Table(np.concatenate([r for r in all_results]))
        else:
            logger.warning("No spectra found")
            return Table()
        
        logger.info(f"Found {len(final_result)} total spectra")
        
        # Save results if requested
        if output_file:
            save_table(final_result, output_file)
            logger.info(f"Saved spectral source results to {output_file}")
        
        return final_result
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.logout()