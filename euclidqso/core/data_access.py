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

try:
    import keyring  # type: ignore
except Exception:  # pragma: no cover
    keyring = None

try:
    from getpass import getpass
except Exception:  # pragma: no cover
    getpass = None


def _flux_to_abmag(flux: np.ndarray, zeropoint: float = 23.9) -> np.ndarray:
    """
    Convert flux to AB magnitude.

    Assumes flux is in microJanskys (uJy) for the default zeropoint of 23.9.
    Non-positive flux values yield NaN magnitudes.

    Parameters
    ----------
    flux : array-like
        Flux values (expected in uJy for zeropoint=23.9).
    zeropoint : float, default 23.9
        AB zeropoint corresponding to flux units.

    Returns
    -------
    np.ndarray
        Magnitudes with NaN where flux <= 0.
    """
    arr = np.asarray(flux, dtype=float)
    mag = np.full(arr.shape, np.nan, dtype=float)
    pos = arr > 0
    if np.any(pos):
        mag[pos] = -2.5 * np.log10(arr[pos]) + zeropoint
    return mag


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
    
    def login(
        self,
        credentials_file: Optional[str] = None,
        user: Optional[str] = None,
        password: Optional[str] = None,
        use_keyring: bool = True,
        prompt: bool = False,
        store_in_keyring: bool = False,
    ):
        """
        Login to Euclid archive.
        
        Parameters
        ----------
        credentials_file : str, optional
            Path to credentials file
        user : str, optional
            Username for interactive login
        password : str, optional
            Password to use (avoid if possible; prefer keyring or prompt)
        use_keyring : bool, default True
            Try OS keyring to retrieve password if available
        prompt : bool, default False
            Prompt for password if not provided and running interactively
        store_in_keyring : bool, default False
            If prompting, store password in keyring for future use
        """
        if credentials_file is None:
            credentials_file = os.environ.get('EUCLIDQSO_CREDENTIALS_FILE') or config.get('data.credentials_file')
        
        try:
            # 1) Credentials file (on-disk). Use only if explicitly provided or configured.
            if credentials_file and Path(credentials_file).exists():
                self.euclid.login(credentials_file=credentials_file)
                logger.info("Successfully logged in with credentials file")
            else:
                # 2) Environment variables
                env_user = os.environ.get('EUCLID_USER')
                env_pass = os.environ.get('EUCLID_PASSWORD')
                use_user = user or env_user

                if use_user and (password or env_pass):
                    self.euclid.login(user=use_user, password=password or env_pass)
                    logger.info(f"Successfully logged in as {use_user} via environment/user credentials")
                else:
                    # 3) OS keyring
                    if use_user and use_keyring and keyring is not None:
                        try:
                            kr_pass = keyring.get_password('euclidqso', use_user)
                        except Exception:
                            kr_pass = None
                        if kr_pass:
                            self.euclid.login(user=use_user, password=kr_pass)
                            logger.info(f"Successfully logged in as {use_user} via keyring")
                            self._logged_in = True
                            return
                    # 4) Prompt if allowed
                    if prompt and use_user and getpass is not None:
                        pw = getpass(f"Password for {use_user}: ")
                        if pw:
                            self.euclid.login(user=use_user, password=pw)
                            logger.info(f"Successfully logged in as {use_user} via prompt")
                            if store_in_keyring and keyring is not None:
                                try:
                                    keyring.set_password('euclidqso', use_user, pw)
                                except Exception:
                                    pass
                        else:
                            logger.warning("Empty password entered; skipping login")
                            return
                    else:
                        # 5) Fall back to interactive flow provided by astroquery (may open browser or prompt)
                        if use_user:
                            self.euclid.login(user=use_user)
                            logger.info(f"Successfully logged in as {use_user}")
                        else:
                            logger.warning("No credentials provided - set EUCLID_USER/EUCLID_PASSWORD or use keyring/prompt")
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
        mer_table: Optional[str] = None,
        output_file: Optional[Union[str, Path]] = None,
        ra_col: str = 'ra',
        dec_col: str = 'dec',
        max_sources: Optional[int] = None,
        use_object_id: Optional[bool] = None,
        idr_field: Optional[str] = None,
        full_async: bool = False,
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
        use_object_id : {None, True, False}, optional
            Matching preference. None (default) auto-detects and uses object_id join
            when a suitable column exists; True forces equality join (warns if missing);
            False forces spatial crossmatch even if object_id is present.
        idr_field : {'WIDE', 'DEEP'}, optional
            IDR field selector. Only used when environment='IDR' and mer_table is not
            explicitly provided. Defaults to 'WIDE'.
        full_async : bool, optional
            If True, submit the entire user table as a single asynchronous job instead of
            processing in 1000-source batches.
            
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
            mer_table = self._get_mer_table_name(idr_field=idr_field)
        
        logger.info(f"Using MER table: {mer_table}")
        logger.info(f"Crossmatching {len(user_data)} sources with radius {radius}\"")

        # Decide matching mode
        user_cols = list(user_data.colnames)
        has_oid = 'object_id' in user_cols
        has_oid_alt = 'object_id_euclid' in user_cols
        want_oid = (
            (use_object_id is True and (has_oid or has_oid_alt))
            or (use_object_id is None and (has_oid or has_oid_alt))
        )
        if use_object_id is True and not (has_oid or has_oid_alt):
            logger.warning("use_object_id=True but no object_id/object_id_euclid found; using spatial match")

        # Check/resolve RA/Dec columns only if spatial matching
        if not want_oid:
            if ra_col not in user_data.colnames:
                # Try alternative column names
                for alt_ra in ['right_ascension', 'RA', 'ra_deg', 'RAJ2000', 'right_ascension_euclid', 'ra_euclid']:
                    if alt_ra in user_data.colnames:
                        ra_col = alt_ra
                        break
                else:
                    raise ValueError(f"RA column '{ra_col}' not found in user table")
            
            if dec_col not in user_data.colnames:
                # Try alternative column names
                for alt_dec in ['declination', 'DEC', 'dec_deg', 'DEJ2000', 'declination_euclid', 'dec_euclid']:
                    if alt_dec in user_data.colnames:
                        dec_col = alt_dec
                        break
                else:
                    raise ValueError(f"Dec column '{dec_col}' not found in user table")
        
        # Perform crossmatch using ADQL
        crossmatch_results = []
        batch_size = 1000  # Process in batches
        
        if full_async:
            batches = [(0, len(user_data))]
        else:
            batches = [(i, min(i + batch_size, len(user_data))) for i in range(0, len(user_data), batch_size)]
        
        total_batches = len(batches)
        if total_batches == 0:
            return Table()
        
        for idx, (start, end) in enumerate(batches, 1):
            batch = user_data[start:end]
            logger.info(f"Processing batch {idx}/{total_batches}")
            
            batch_result = self._crossmatch_batch(
                batch, ra_col, dec_col, radius, mer_table, use_object_id, force_async=full_async
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
    
    def _get_mer_table_name(self, idr_field: Optional[str] = None) -> str:
        """Get MER catalogue table name for current environment."""
        if self.environment == 'IDR':
            field = (idr_field or 'WIDE').upper()
            valid_fields = {'WIDE', 'DEEP'}
            if field not in valid_fields:
                raise ValueError(
                    f"Invalid IDR field '{idr_field}'. Expected one of {sorted(valid_fields)}"
                )
            return (
                'catalogue.mer_catalogue_wide'
                if field == 'WIDE'
                else 'catalogue.mer_catalogue_deep'
            )
        
        table_names = {
            'PDR': 'catalogue.mer_catalogue',
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
        mer_table: str,
        use_object_id: Optional[bool] = None,
        force_async: bool = False,
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
            user_cols = list(batch.colnames)
            has_oid = 'object_id' in user_cols
            has_oid_alt = 'object_id_euclid' in user_cols
            # Determine if we want to use object_id join
            want_oid = (use_object_id is True) or (use_object_id is None and (has_oid or has_oid_alt))

            # Columns we commonly want from MER
            mer_columns = [
                ('object_id', 'object_id'),
                ('right_ascension', 'mer_ra'),
                ('declination', 'mer_dec'),
                ('mu_max', 'mu_max'),
                ('mumax_minus_mag', 'mumax_minus_mag'),
                ('kron_radius', 'kron_radius'),
                ('kron_radius_err', 'kron_radius_err'),
                ('gaia_id', 'gaia_id'),
                ('gaia_match_quality', 'gaia_match_quality'),
                ('flux_y_templfit', 'flux_y_templfit'),
                ('flux_h_templfit', 'flux_h_templfit'),
                ('flux_j_templfit', 'flux_j_templfit'),
                ('flux_u_ext_decam_templfit', 'flux_u_ext_decam_templfit'),
                ('flux_g_ext_decam_templfit', 'flux_g_ext_decam_templfit'),
                ('flux_r_ext_decam_templfit', 'flux_r_ext_decam_templfit'),
                ('flux_i_ext_decam_templfit', 'flux_i_ext_decam_templfit'),
                ('flux_z_ext_decam_templfit', 'flux_z_ext_decam_templfit'),
                ('flux_u_ext_megacam_templfit', 'flux_u_ext_megacam_templfit'),
                ('flux_r_ext_megacam_templfit', 'flux_r_ext_megacam_templfit'),
                ('flux_g_ext_jpcam_templfit', 'flux_g_ext_jpcam_templfit'),
                ('flux_i_ext_panstarrs_templfit', 'flux_i_ext_panstarrs_templfit'),
                ('flux_z_ext_panstarrs_templfit', 'flux_z_ext_panstarrs_templfit'),
                ('flux_g_ext_hsc_templfit', 'flux_g_ext_hsc_templfit'),
                ('flux_z_ext_hsc_templfit', 'flux_z_ext_hsc_templfit'),
                ('fluxerr_y_templfit', 'fluxerr_y_templfit'),
                ('fluxerr_j_templfit', 'fluxerr_j_templfit'),
                ('fluxerr_h_templfit', 'fluxerr_h_templfit'),
                ('fluxerr_r_ext_decam_templfit', 'fluxerr_r_ext_decam_templfit'),
                ('fluxerr_i_ext_decam_templfit', 'fluxerr_i_ext_decam_templfit'),
                ('fluxerr_z_ext_decam_templfit', 'fluxerr_z_ext_decam_templfit'),
                ('fluxerr_u_ext_megacam_templfit', 'fluxerr_u_ext_megacam_templfit'),
                ('fluxerr_r_ext_megacam_templfit', 'fluxerr_r_ext_megacam_templfit'),
                ('fluxerr_g_ext_jpcam_templfit', 'fluxerr_g_ext_jpcam_templfit'),
                ('fluxerr_i_ext_panstarrs_templfit', 'fluxerr_i_ext_panstarrs_templfit'),
                ('fluxerr_z_ext_panstarrs_templfit', 'fluxerr_z_ext_panstarrs_templfit'),
                ('fluxerr_g_ext_hsc_templfit', 'fluxerr_g_ext_hsc_templfit'),
                ('fluxerr_z_ext_hsc_templfit', 'fluxerr_z_ext_hsc_templfit'),
                ('fluxerr_u_ext_decam_templfit', 'fluxerr_u_ext_decam_templfit'),
                ('fluxerr_g_ext_decam_templfit', 'fluxerr_g_ext_decam_templfit'),
                ('flux_vis_psf', 'flux_vis_psf'),
                ('fluxerr_vis_psf', 'fluxerr_vis_psf'),
                ('segmentation_map_id', 'segmentation_map_id'),
                ('segmentation_area', 'segmentation_area'),
            ]

            if want_oid and (has_oid or has_oid_alt):
                # Equality join on object_id
                user_id_col = 'object_id' if has_oid else 'object_id_euclid'

                # Build user select list, avoiding name collision with MER object_id
                user_aliases = set()
                user_select_parts = []
                for col in user_cols:
                    if col == 'object_id':
                        alias = 'object_id_user'
                    elif col == 'object_id_euclid':
                        alias = 'object_id_euclid'  # keep as-is
                    else:
                        alias = col
                    user_aliases.add(alias)
                    user_select_parts.append(f"u.{col} AS {alias}")

                # Build MER select list, prefix alias if it would duplicate a user alias
                mer_select_parts = []
                for mcol, alias in mer_columns:
                    out_alias = alias
                    if out_alias in user_aliases:
                        out_alias = f"mer_{out_alias}"
                    mer_select_parts.append(f"m.{mcol} AS {out_alias}")

                select_list = user_select_parts + mer_select_parts
                query = (
                    f"SELECT {', '.join(select_list)}\n"
                    f"FROM TAP_UPLOAD.{upload_name} AS u\n"
                    f"JOIN {mer_table} AS m ON u.{user_id_col} = m.object_id\n"
                    f"ORDER BY u.{user_id_col}"
                )
            else:
                if use_object_id is True and not (has_oid or has_oid_alt):
                    logger.warning("use_object_id=True but no object_id column found; falling back to spatial match")
                # Spatial crossmatch via distance predicate
                radius_deg = radius / 3600.0  # Convert to degrees
                distance_expr = f"DISTANCE(u.{ra_col}, u.{dec_col}, m.right_ascension, m.declination)"
                query = f"""
                SELECT u.*, 
                       m.object_id, 
                       m.right_ascension AS mer_ra, 
                       m.declination AS mer_dec,
                       mu_max, mumax_minus_mag, kron_radius, kron_radius_err, gaia_id, gaia_match_quality, 
                       flux_y_templfit, flux_h_templfit, flux_j_templfit, flux_u_ext_decam_templfit, flux_g_ext_decam_templfit, 
                       flux_r_ext_decam_templfit, flux_i_ext_decam_templfit, flux_z_ext_decam_templfit, 
                       flux_u_ext_megacam_templfit, flux_r_ext_megacam_templfit, flux_g_ext_jpcam_templfit, 
                       flux_i_ext_panstarrs_templfit, flux_z_ext_panstarrs_templfit, flux_g_ext_hsc_templfit, 
                       flux_z_ext_hsc_templfit, fluxerr_y_templfit, fluxerr_j_templfit, fluxerr_h_templfit, 
                       fluxerr_r_ext_decam_templfit, fluxerr_i_ext_decam_templfit, fluxerr_z_ext_decam_templfit, 
                       fluxerr_u_ext_megacam_templfit, fluxerr_r_ext_megacam_templfit, fluxerr_g_ext_jpcam_templfit, 
                       fluxerr_i_ext_panstarrs_templfit, fluxerr_z_ext_panstarrs_templfit, fluxerr_g_ext_hsc_templfit, 
                       fluxerr_z_ext_hsc_templfit, fluxerr_u_ext_decam_templfit, fluxerr_g_ext_decam_templfit, 
                       flux_vis_psf, fluxerr_vis_psf, m.segmentation_map_id, m.segmentation_area,
                       {distance_expr} AS separation_deg
                FROM TAP_UPLOAD.{upload_name} AS u
                JOIN {mer_table} AS m 
                ON {distance_expr} < {radius_deg}
                ORDER BY u.{ra_col}, separation_deg
                """
            
            # Execute query with temporary table upload
            if force_async or len(batch) >= 2000:
                job = self.euclid.launch_job_async(
                    query, upload_resource=tmp_name, upload_table_name=upload_name
                )
            else:
                job = self.euclid.launch_job(
                    query, upload_resource=tmp_name, upload_table_name=upload_name
                )
            
            result = job.get_results()
            
            # Convert separation from degrees to arcseconds
            if len(result) > 0 and 'separation_deg' in result.colnames:
                result['separation_arcsec'] = result['separation_deg'] * 3600.0
                # Remove the degree column
                result.remove_column('separation_deg')

            # Add AB magnitudes for selected flux columns when present
            if len(result) > 0:
                flux_to_mag_map = {
                    'flux_y_templfit': 'mag_y_templfit',
                    'flux_h_templfit': 'mag_h_templfit',
                    'flux_j_templfit': 'mag_j_templfit',
                    'flux_vis_psf': 'mag_vis_psf',
                }
                for flux_col, mag_col in flux_to_mag_map.items():
                    if flux_col in result.colnames and mag_col not in result.colnames:
                        try:
                            result[mag_col] = _flux_to_abmag(result[flux_col])
                        except Exception as e:
                            logger.warning(f"Failed to compute {mag_col} from {flux_col}: {e}")
            
            return result
            
        finally:
            # Clean up temporary file
            if os.path.exists(tmp_name):
                os.unlink(tmp_name)
    
    def _query_spectra_batch(self, batch: Table) -> Table:
        """Query spectra for a batch of object IDs."""
        
        # Upload user table to archive for querying
        with tempfile.NamedTemporaryFile(suffix='.vot', delete=False) as tmp_file:
            batch.write(tmp_file.name, format='votable', overwrite=True)
            tmp_name = tmp_file.name
        
        try:
            # Use temporary upload with launch_job
            upload_name = f"user_spectra_batch_{np.random.randint(10000, 99999)}"
            
            # Construct spectra query
            query = f"""
            SELECT s.source_id, s.ra_obj, s.dec_obj, s.datalabs_path, 
                   s.file_name, s.hdu_index,
                   u.object_id
            FROM TAP_UPLOAD.{upload_name} AS u
            JOIN sedm.spectra_source AS s
            ON s.source_id = u.object_id
            WHERE s.datalabs_path IS NOT NULL
            ORDER BY s.source_id
            """
            
            # Execute query with temporary table upload
            if len(batch) < 2000:
                job = self.euclid.launch_job(query, upload_resource=tmp_name, 
                                           upload_table_name=upload_name)
            else:
                job = self.euclid.launch_job_async(query, upload_resource=tmp_name,
                                                 upload_table_name=upload_name)
            
            result = job.get_results()
            return result
            
        except Exception as e:
            logger.error(f"Error querying spectra batch: {e}")
            return Table()
            
        finally:
            # Clean up temporary file
            if os.path.exists(tmp_name):
                os.unlink(tmp_name)
    
    def query_spectra_sources(
        self,
        crossmatch_table: Optional[Table] = None,
        output_file: Optional[Union[str, Path]] = None
    ) -> Table:
        """
        Query spectral sources from Euclid archive.
        
        Parameters
        ----------
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
        if crossmatch_table is None:
            raise ValueError("Must provide crossmatch_table with 'object_id' or 'object_id_euclid'")
        # Accept either object_id or object_id_euclid
        if 'object_id' in crossmatch_table.colnames:
            object_ids = list(set([int(x) for x in crossmatch_table['object_id']]))
        elif 'object_id_euclid' in crossmatch_table.colnames:
            object_ids = list(set([int(x) for x in crossmatch_table['object_id_euclid']]))
        else:
            raise ValueError("crossmatch_table must contain 'object_id' or 'object_id_euclid'")
        
        logger.info(f"Querying spectra for {len(object_ids)} unique objects")
        
        # Create user table with object IDs
        user_table = Table()
        user_table['object_id'] = object_ids
        
        # Query in batches to avoid query limits  
        batch_size = 1000
        all_results = []
        
        for i in range(0, len(user_table), batch_size):
            batch = user_table[i:i+batch_size]
            logger.info(f"Querying batch {i//batch_size + 1}/{(len(user_table)-1)//batch_size + 1}")
            
            batch_result = self._query_spectra_batch(batch)
            
            if len(batch_result) > 0:
                all_results.append(batch_result)
                logger.info(f"Found {len(batch_result)} spectra in batch")
        
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
    
    def get_individual_spectrum(
        self,
        datalabs_path: str,
        file_name: str,
        hdu_index: int
    ):
        """
        Get an individual spectrum HDU from a combined spectra file.
        
        This is a convenience wrapper around SpectrumLoader.
        
        Parameters
        ----------
        datalabs_path : str
            Path to the data directory on ESA Datalabs
        file_name : str
            Name of the combined spectra FITS file
        hdu_index : int
            HDU index of the specific spectrum
            
        Returns
        -------
        fits.HDU
            Individual spectrum HDU
        """
        from euclidqso.core.spectra import SpectrumLoader
        
        loader = SpectrumLoader()
        file_path = loader.verify_spectrum_path(datalabs_path, file_name)
        if file_path is None:
            raise FileNotFoundError(f"Spectra file not found: {os.path.join(datalabs_path, file_name)}")
        
        return loader.load_spectrum(file_path, hdu_index)
    
    def combine_spectra_to_fits(
        self,
        spectra_table: Table,
        output_file: Union[str, Path],
        source_id_col: str = 'source_id',
        datalabs_path_col: str = 'datalabs_path',
        file_name_col: str = 'file_name',
        hdu_index_col: str = 'hdu_index',
        max_spectra: Optional[int] = None
    ) -> str:
        """
        Combine individual spectra into a single FITS file.
        
        This function replicates the functionality from cell 23 of the 
        Spectra_visualization_catglobe.ipynb notebook, using the 
        existing SpectrumCompiler infrastructure.
        
        Parameters
        ----------
        spectra_table : Table
            Table containing spectra source information from query_spectra_sources
        output_file : str or Path
            Output FITS file path
        source_id_col : str, default 'source_id'
            Column name for source IDs
        datalabs_path_col : str, default 'datalabs_path'
            Column name for datalabs paths
        file_name_col : str, default 'file_name'
            Column name for file names
        hdu_index_col : str, default 'hdu_index'
            Column name for HDU indices
        max_spectra : int, optional
            Maximum number of spectra to include (for testing)
            
        Returns
        -------
        str
            Path to the created FITS file
        """
        from euclidqso.core.spectra import SpectrumCompiler
        
        # Limit number of spectra if requested
        if max_spectra is not None:
            spectra_subset = spectra_table[:max_spectra]
        else:
            spectra_subset = spectra_table
        
        logger.info(f"Combining {len(spectra_subset)} spectra into {output_file}")
        
        # Use SpectrumCompiler's convenience method for single-file output
        compiler = SpectrumCompiler()
        
        return compiler.compile_single_fits(
            spectra_subset,
            output_file=output_file,
            source_id_col=source_id_col,
            datalabs_path_col=datalabs_path_col,
            file_name_col=file_name_col,
            hdu_index_col=hdu_index_col,
            overwrite=True
        )
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.logout()
