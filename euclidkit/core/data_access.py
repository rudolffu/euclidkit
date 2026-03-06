"""
Data access module for euclidkit package.

Provides interfaces to Euclid science archive and data volumes.
"""

import os
import logging
import json
from datetime import datetime
from pathlib import Path
from typing import Union, Optional, List, Dict, Any, Tuple
import tempfile

import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

from euclidkit.config import config
from euclidkit.utils.io import load_table, save_table

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
        self.euclid = None
        self._logged_in = False
        
        logger.info(f"Initialized EuclidArchive for {environment} environment")

    def _ensure_client(self):
        """Lazily create the Euclid client to avoid network calls at import time."""
        if self.euclid is None:
            from astroquery.esa.euclid.core import EuclidClass
            self.euclid = EuclidClass(environment=self.environment)
        return self.euclid
    
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
        euclid = self._ensure_client()

        if credentials_file is None:
            credentials_file = os.environ.get('EUCLIDKIT_CREDENTIALS_FILE') or config.get('data.credentials_file')
        
        try:
            # 1) Credentials file (on-disk). Use only if explicitly provided or configured.
            if credentials_file and Path(credentials_file).exists():
                euclid.login(credentials_file=credentials_file)
                logger.info("Successfully logged in with credentials file")
            else:
                # 2) Environment variables
                env_user = os.environ.get('EUCLID_USER')
                env_pass = os.environ.get('EUCLID_PASSWORD')
                use_user = user or env_user

                if use_user and (password or env_pass):
                    euclid.login(user=use_user, password=password or env_pass)
                    logger.info(f"Successfully logged in as {use_user} via environment/user credentials")
                else:
                    # 3) OS keyring
                    if use_user and use_keyring and keyring is not None:
                        try:
                            kr_pass = keyring.get_password('euclidkit', use_user)
                        except Exception:
                            kr_pass = None
                        if kr_pass:
                            euclid.login(user=use_user, password=kr_pass)
                            logger.info(f"Successfully logged in as {use_user} via keyring")
                            self._logged_in = True
                            return
                    # 4) Prompt if allowed
                    if prompt and use_user and getpass is not None:
                        pw = getpass(f"Password for {use_user}: ")
                        if pw:
                            euclid.login(user=use_user, password=pw)
                            logger.info(f"Successfully logged in as {use_user} via prompt")
                            if store_in_keyring and keyring is not None:
                                try:
                                    keyring.set_password('euclidkit', use_user, pw)
                                except Exception:
                                    pass
                        else:
                            logger.warning("Empty password entered; skipping login")
                            return
                    else:
                        # 5) Fall back to interactive flow provided by astroquery (may open browser or prompt)
                        if use_user:
                            euclid.login(user=use_user)
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
            if self.euclid is None:
                return
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
        async_chunk_size: int = 500000,
    ) -> Union[Table, Dict[str, Any]]:
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
            When ``True``, ``source_id`` is also accepted as a fallback join key
            to MER ``object_id`` if ``object_id`` columns are absent.
        idr_field : {'WIDE', 'DEEP'}, optional
            IDR field selector. Only used when environment='IDR' and mer_table is not
            explicitly provided. Defaults to 'WIDE'.
        full_async : bool, optional
            If True, use asynchronous TAP submission mode.
            For tables larger than ``async_chunk_size``, the table is split into
            multiple async jobs and results are downloaded/merged automatically.
        async_chunk_size : int, optional
            Chunk size (rows per async job) used when ``full_async=True`` for
            large input tables. Default is 500000.
            
        Returns
        -------
        Table or dict
            Crossmatched results table, or async job metadata dictionary when
            ``full_async`` is True.
        """
        self._ensure_client()
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

        # Decide matching mode
        user_cols = list(user_data.colnames)
        has_oid = 'object_id' in user_cols
        has_oid_alt = 'object_id_euclid' in user_cols
        has_source_id = 'source_id' in user_cols
        want_oid = (
            (use_object_id is True and (has_oid or has_oid_alt or has_source_id))
            or (use_object_id is None and (has_oid or has_oid_alt))
        )
        if use_object_id is True and not (has_oid or has_oid_alt or has_source_id):
            logger.warning("use_object_id=True but no object_id/object_id_euclid/source_id found; using spatial match")
        elif use_object_id is True and has_source_id and not (has_oid or has_oid_alt):
            logger.info("use_object_id=True and only source_id present; joining source_id to MER object_id")

        logger.info(f"Using MER table: {mer_table}")
        if want_oid:
            logger.info(f"Crossmatching {len(user_data)} sources in object-id mode")
        else:
            logger.info(f"Crossmatching {len(user_data)} sources with radius {radius} arcsec")

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
        if full_async:
            if len(user_data) == 0:
                raise ValueError("User table is empty; nothing to submit in full_async mode.")
            if output_file is None:
                raise ValueError("output_file must be provided when full_async=True to store job metadata.")
            if async_chunk_size <= 0:
                raise ValueError("async_chunk_size must be a positive integer")

            # For very large tables, run async chunk jobs and merge downloaded results.
            if len(user_data) > async_chunk_size:
                logger.info(
                    "full_async enabled with %d rows; using chunked async mode with chunk_size=%d",
                    len(user_data), async_chunk_size
                )
                batches = [
                    (i, min(i + async_chunk_size, len(user_data)))
                    for i in range(0, len(user_data), async_chunk_size)
                ]
                chunk_results = []
                job_ids = []
                for idx, (start, end) in enumerate(batches, 1):
                    logger.info("Submitting async chunk %d/%d (%d rows)", idx, len(batches), end - start)
                    batch = user_data[start:end]
                    submission = self._crossmatch_batch(
                        batch,
                        ra_col,
                        dec_col,
                        radius,
                        mer_table,
                        use_object_id,
                        force_async=True,
                        fetch_results=False,
                    )
                    job = submission.get('job')
                    job_id = getattr(job, 'jobid', None)
                    if job_id:
                        job_ids.append(job_id)

                    try:
                        batch_result = job.get_results()
                        if batch_result is not None and len(batch_result) > 0:
                            chunk_results.append(batch_result)
                    except Exception as exc:
                        logger.warning("Failed to fetch async results for chunk %d: %s", idx, exc)

                if chunk_results:
                    final_result = Table(np.concatenate([r for r in chunk_results]))
                else:
                    final_result = Table()

                save_table(final_result, output_file)
                logger.info(
                    "Saved merged chunked async results (%d rows) to %s",
                    len(final_result), output_file
                )
                return {
                    'type': 'euclidkit_crossmatch_async_chunked',
                    'environment': self.environment,
                    'submitted_at_utc': datetime.utcnow().isoformat(timespec='seconds') + 'Z',
                    'row_count': len(user_data),
                    'chunk_size': async_chunk_size,
                    'chunk_count': len(batches),
                    'job_ids': job_ids,
                    'results_downloaded': True,
                    'result_row_count': len(final_result),
                    'output_file': str(output_file),
                    'idr_field': idr_field.upper() if idr_field else None,
                }

            submission = self._crossmatch_batch(
                user_data,
                ra_col,
                dec_col,
                radius,
                mer_table,
                use_object_id,
                force_async=True,
                fetch_results=False,
            )

            job_info = self._build_async_job_info(
                submission,
                mer_table=mer_table,
                radius=radius,
                ra_col=ra_col,
                dec_col=dec_col,
                idr_field=idr_field,
                row_count=len(user_data),
            )
            job = submission.get('job')
            results_downloaded = False
            result_row_count = None
            download_error = None

            # Try to fetch async results and write the requested output table.
            if job is not None:
                try:
                    logger.info(
                        "Attempting to download results for async job %s",
                        getattr(job, 'jobid', None),
                    )
                    async_result = job.get_results()
                    if async_result is not None:
                        save_table(async_result, output_file)
                        results_downloaded = True
                        result_row_count = len(async_result)
                        logger.info(
                            "Downloaded async results (%d rows) to %s",
                            result_row_count,
                            output_file,
                        )
                except Exception as exc:
                    download_error = str(exc)
                    logger.warning(
                        "Async job submitted but result download failed: %s",
                        exc,
                    )

            job_info['results_downloaded'] = results_downloaded
            job_info['result_row_count'] = result_row_count
            job_info['download_error'] = download_error

            if not results_downloaded:
                output_path = Path(output_file)
                if output_path.suffix.lower() == '.json':
                    job_info_path = output_path
                else:
                    job_info_path = output_path.with_name(output_path.name + '.job.json')
                self._write_job_info(job_info, job_info_path)
                job_info['job_info_file'] = str(job_info_path)
                logger.info(f"Saved async job info to {job_info_path}")

            return job_info

        crossmatch_results = []
        batch_size = 1000  # Process in batches
        use_async_batches = len(user_data) > 2000
        if use_async_batches:
            logger.info(
                "Input has %d rows (>2000); forcing asynchronous TAP jobs for all batches",
                len(user_data),
            )
        
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
                batch,
                ra_col,
                dec_col,
                radius,
                mer_table,
                use_object_id,
                force_async=use_async_batches,
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

    def upload_user_table(
        self,
        table: Union[str, Path, Table],
        table_name: str,
        description: Optional[str] = None,
        fmt: Optional[str] = None,
        overwrite: bool = False,
        verbose: bool = False,
    ) -> Dict[str, Any]:
        """
        Upload a table to the Euclid user workspace.

        Parameters
        ----------
        table : str, Path, or astropy Table
            Resource to upload. File paths are passed directly to TAP.
        table_name : str
            Destination table name (without schema qualifiers).
        description : str, optional
            Table description stored alongside the upload.
        fmt : str, optional
            Input format (e.g., 'votable', 'fits', 'csv'). If omitted and a file
            path is provided, the format is inferred from the extension.
        overwrite : bool, default False
            If True, delete any existing table with the same name before upload.
        verbose : bool, default False
            Emit TAP upload status messages.

        Returns
        -------
        dict
            Metadata describing the upload request, including any job ID.
        """
        self._ensure_client()
        if not table_name:
            raise ValueError("table_name is required for uploads")
        table_name = table_name.strip()
        if '.' in table_name:
            raise ValueError("table_name must not contain '.' (schema qualifiers are not allowed)")

        if not self._logged_in:
            self.login()

        resource = table
        inferred_format = fmt

        if isinstance(table, (str, Path)):
            path = Path(table)
            if not path.exists():
                raise FileNotFoundError(f"Table file not found: {path}")
            resource = str(path)
            if inferred_format is None:
                inferred_format = self._infer_upload_format(path.suffix.lower())
        else:
            # astropy Table upload
            if inferred_format is None:
                inferred_format = 'votable'

        if overwrite:
            try:
                self.euclid.delete_user_table(table_name=table_name, force_removal=True, verbose=verbose)
                logger.info(f"Removed existing table '{table_name}' prior to upload")
            except Exception as exc:  # pragma: no cover - network errors
                logger.warning(f"Failed to delete existing table '{table_name}': {exc}")

        job = self.euclid.upload_table(
            upload_resource=resource,
            table_name=table_name,
            table_description=description,
            format=inferred_format,
            verbose=verbose,
        )

        job_info = {
            'table_name': table_name,
            'format': inferred_format,
            'job_id': getattr(job, 'jobid', None),
            'job_phase': getattr(job, 'phase', None) or getattr(job, 'get_phase', lambda: None)(),
            'description': description,
            'overwrite': overwrite,
            'resource_type': 'table' if isinstance(table, Table) else 'file',
        }

        if job_info['job_id']:
            logger.info(f"Upload job '{job_info['job_id']}' created for table '{table_name}'")
        else:
            job_info['job_phase'] = job_info['job_phase'] or 'COMPLETED'
            logger.info(f"Uploaded table '{table_name}' synchronously")

        return job_info

    def _resolve_user_table_reference(self, user_table_name: str) -> str:
        """Resolve short user table names to fully qualified TAP schema names."""
        table_name = user_table_name.strip()
        if '.' in table_name:
            return table_name

        # astroquery TapPlus stores the logged-in username on this private attribute.
        username = getattr(self.euclid, "_TapPlus__user", None)
        if not username:
            username = os.environ.get("EUCLID_USER")
        if not username:
            raise ValueError(
                "Could not resolve username for user table. "
                "Provide fully qualified name like user_<username>.<table>."
            )
        return f"user_{username}.{table_name}"

    def _get_remote_table_columns(self, qualified_table_name: str) -> List[str]:
        """Fetch column names for a remote TAP table from TAP_SCHEMA."""
        if '.' not in qualified_table_name:
            raise ValueError("qualified_table_name must include schema and table name")
        schema_name, table_name = qualified_table_name.split('.', 1)
        query = (
            "SELECT column_name "
            "FROM TAP_SCHEMA.columns "
            f"WHERE lower(schema_name) = lower('{schema_name}') "
            f"AND lower(table_name) = lower('{table_name}')"
        )
        job = self.euclid.launch_job(query)
        result = job.get_results()
        return [str(name).lower() for name in result['column_name']]

    def crossmatch_user_table(
        self,
        user_table_name: str,
        radius: float = 1.0,
        mer_table: Optional[str] = None,
        output_file: Optional[Union[str, Path]] = None,
        ra_col: str = 'ra',
        dec_col: str = 'dec',
        max_sources: Optional[int] = None,
        use_object_id: Optional[bool] = None,
        idr_field: Optional[str] = None,
        full_async: bool = False,
    ) -> Union[Table, Dict[str, Any]]:
        """
        Crossmatch an existing archive user table against MER without re-upload.

        Parameters are similar to ``crossmatch_sources`` but ``user_table_name``
        refers to an already-uploaded TAP table.
        """
        self._ensure_client()
        if not self._logged_in:
            logger.warning("Not logged in - attempting login with default credentials")
            self.login()

        if mer_table is None:
            mer_table = self._get_mer_table_name(idr_field=idr_field)

        user_table_ref = self._resolve_user_table_reference(user_table_name)
        user_cols = self._get_remote_table_columns(user_table_ref)
        if not user_cols:
            raise ValueError(f"No columns found for remote table {user_table_ref}")
        if max_sources is not None and max_sources <= 0:
            raise ValueError("max_sources must be a positive integer")

        has_oid = 'object_id' in user_cols
        has_oid_alt = 'object_id_euclid' in user_cols
        has_source_id = 'source_id' in user_cols
        want_oid = (
            (use_object_id is True and (has_oid or has_oid_alt or has_source_id))
            or (use_object_id is None and (has_oid or has_oid_alt))
        )

        if use_object_id is True and not (has_oid or has_oid_alt or has_source_id):
            logger.warning("use_object_id=True but no object_id/object_id_euclid/source_id found; using spatial match")
        elif use_object_id is True and has_source_id and not (has_oid or has_oid_alt):
            logger.info("use_object_id=True and only source_id present; joining source_id to MER object_id")

        # Resolve RA/Dec columns for spatial mode only.
        if not want_oid:
            ra_candidates = [ra_col.lower(), 'right_ascension', 'ra', 'ra_deg', 'raj2000', 'right_ascension_euclid', 'ra_euclid']
            dec_candidates = [dec_col.lower(), 'declination', 'dec', 'dec_deg', 'dej2000', 'declination_euclid', 'dec_euclid']
            ra_use = next((c for c in ra_candidates if c in user_cols), None)
            dec_use = next((c for c in dec_candidates if c in user_cols), None)
            if not ra_use or not dec_use:
                raise ValueError(
                    f"Could not find RA/Dec columns in remote table {user_table_ref}. "
                    f"Available columns include: {user_cols[:20]}"
                )
            ra_col = ra_use
            dec_col = dec_use

        logger.info(f"Using MER table: {mer_table}")
        logger.info(f"Using remote user table: {user_table_ref}")
        if want_oid:
            logger.info("Crossmatching remote table in object-id mode")
        else:
            logger.info(f"Crossmatching remote table with radius {radius} arcsec")

        mer_oid_alias = "mer_object_id" if has_oid else "object_id"
        mer_select_parts = [
            f"m.object_id AS {mer_oid_alias}",
            "m.right_ascension AS mer_ra",
            "m.declination AS mer_dec",
            "m.mu_max AS mu_max",
            "m.mumax_minus_mag AS mumax_minus_mag",
            "m.point_like_prob AS point_like_prob",
            "m.extended_prob AS extended_prob",
            "m.kron_radius AS kron_radius",
            "m.kron_radius_err AS kron_radius_err",
            "m.gaia_id AS gaia_id",
            "m.gaia_match_quality AS gaia_match_quality",
            "m.flux_y_templfit AS flux_y_templfit",
            "m.flux_h_templfit AS flux_h_templfit",
            "m.flux_j_templfit AS flux_j_templfit",
            "m.flux_u_ext_decam_templfit AS flux_u_ext_decam_templfit",
            "m.flux_g_ext_decam_templfit AS flux_g_ext_decam_templfit",
            "m.flux_r_ext_decam_templfit AS flux_r_ext_decam_templfit",
            "m.flux_i_ext_decam_templfit AS flux_i_ext_decam_templfit",
            "m.flux_z_ext_decam_templfit AS flux_z_ext_decam_templfit",
            "m.flux_u_ext_megacam_templfit AS flux_u_ext_megacam_templfit",
            "m.flux_r_ext_megacam_templfit AS flux_r_ext_megacam_templfit",
            "m.flux_g_ext_jpcam_templfit AS flux_g_ext_jpcam_templfit",
            "m.flux_i_ext_panstarrs_templfit AS flux_i_ext_panstarrs_templfit",
            "m.flux_z_ext_panstarrs_templfit AS flux_z_ext_panstarrs_templfit",
            "m.flux_g_ext_hsc_templfit AS flux_g_ext_hsc_templfit",
            "m.flux_z_ext_hsc_templfit AS flux_z_ext_hsc_templfit",
            "m.fluxerr_y_templfit AS fluxerr_y_templfit",
            "m.fluxerr_j_templfit AS fluxerr_j_templfit",
            "m.fluxerr_h_templfit AS fluxerr_h_templfit",
            "m.fluxerr_r_ext_decam_templfit AS fluxerr_r_ext_decam_templfit",
            "m.fluxerr_i_ext_decam_templfit AS fluxerr_i_ext_decam_templfit",
            "m.fluxerr_z_ext_decam_templfit AS fluxerr_z_ext_decam_templfit",
            "m.fluxerr_u_ext_megacam_templfit AS fluxerr_u_ext_megacam_templfit",
            "m.fluxerr_r_ext_megacam_templfit AS fluxerr_r_ext_megacam_templfit",
            "m.fluxerr_g_ext_jpcam_templfit AS fluxerr_g_ext_jpcam_templfit",
            "m.fluxerr_i_ext_panstarrs_templfit AS fluxerr_i_ext_panstarrs_templfit",
            "m.fluxerr_z_ext_panstarrs_templfit AS fluxerr_z_ext_panstarrs_templfit",
            "m.fluxerr_g_ext_hsc_templfit AS fluxerr_g_ext_hsc_templfit",
            "m.fluxerr_z_ext_hsc_templfit AS fluxerr_z_ext_hsc_templfit",
            "m.fluxerr_u_ext_decam_templfit AS fluxerr_u_ext_decam_templfit",
            "m.fluxerr_g_ext_decam_templfit AS fluxerr_g_ext_decam_templfit",
            "m.flux_vis_psf AS flux_vis_psf",
            "m.fluxerr_vis_psf AS fluxerr_vis_psf",
            "m.segmentation_map_id AS segmentation_map_id",
            "m.segmentation_area AS segmentation_area",
        ]

        if want_oid:
            if has_oid:
                user_id_col = 'object_id'
            elif has_oid_alt:
                user_id_col = 'object_id_euclid'
            else:
                user_id_col = 'source_id'
            top_clause = f"TOP {int(max_sources)} " if max_sources is not None else ""
            query = (
                f"SELECT {top_clause}u.*, {', '.join(mer_select_parts)} "
                f"FROM {user_table_ref} AS u "
                f"JOIN {mer_table} AS m ON u.{user_id_col} = m.object_id "
                f"ORDER BY u.{user_id_col}"
            )
        else:
            radius_deg = radius / 3600.0
            top_clause = f"TOP {int(max_sources)} " if max_sources is not None else ""
            query = (
                f"SELECT {top_clause}u.*, {', '.join(mer_select_parts)} "
                f"FROM {user_table_ref} AS u "
                f"JOIN {mer_table} AS m "
                f"ON DISTANCE(u.{ra_col}, u.{dec_col}, m.right_ascension, m.declination) < {radius_deg} "
                f"ORDER BY u.{ra_col}"
            )

        if full_async:
            job = self.euclid.launch_job_async(query)
            result = job.get_results()
            if output_file:
                save_table(result, output_file)
                logger.info(f"Saved async crossmatch results to {output_file}")
            return {
                'type': 'euclidkit_crossmatch_remote_async',
                'environment': self.environment,
                'job_id': getattr(job, 'jobid', None),
                'results_downloaded': True,
                'result_row_count': len(result) if result is not None else 0,
                'output_file': str(output_file) if output_file else None,
                'query': query,
            }

        job = self.euclid.launch_job(query)
        result = job.get_results()
        if output_file:
            save_table(result, output_file)
            logger.info(f"Saved crossmatch results to {output_file}")
        return result
    
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

    def _infer_upload_format(self, suffix: str) -> str:
        """Infer TAP upload format from filename suffix."""
        mapping = {
            '.fits': 'fits',
            '.fit': 'fits',
            '.csv': 'csv',
            '.vot': 'votable',
            '.votable': 'votable',
            '.xml': 'votable',
            '.ecsv': 'csv',
            '.txt': 'csv',
        }
        return mapping.get(suffix.lower(), 'votable')

    def _get_spectra_source_table_name(self, idr_field: Optional[str] = None) -> str:
        """Get spectra source table name for current environment."""
        if self.environment == 'IDR':
            field = (idr_field or 'WIDE').upper()
            valid_fields = {'WIDE', 'DEEP'}
            if field not in valid_fields:
                raise ValueError(
                    f"Invalid IDR field '{idr_field}'. Expected one of {sorted(valid_fields)}"
                )
            return (
                'catalogue.spectra_source_wide'
                if field == 'WIDE'
                else 'catalogue.spectra_source_deep'
            )

        prefix_by_env = {
            'PDR': 'q1',
            'OTF': 'q1',
            'REG': 'dr1',
        }
        prefix = prefix_by_env.get(self.environment, 'q1')
        return f"{prefix}.spectra_source"
    
    def _crossmatch_batch(
        self, 
        batch: Table, 
        ra_col: str, 
        dec_col: str, 
        radius: float, 
        mer_table: str,
        use_object_id: Optional[bool] = None,
        force_async: bool = False,
        fetch_results: bool = True,
    ) -> Union[Table, Dict[str, Any]]:
        """Crossmatch a batch of sources."""
        
        self._ensure_client()

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
            has_source_id = 'source_id' in user_cols
            # Determine if we want to use object_id join
            want_oid = (
                (use_object_id is True and (has_oid or has_oid_alt or has_source_id))
                or (use_object_id is None and (has_oid or has_oid_alt))
            )

            # Columns we commonly want from MER
            mer_columns = [
                ('object_id', 'object_id'),
                ('right_ascension', 'mer_ra'),
                ('declination', 'mer_dec'),
                ('mu_max', 'mu_max'),
                ('mumax_minus_mag', 'mumax_minus_mag'),
                ('point_like_prob', 'point_like_prob'),
                ('extended_prob', 'extended_prob'),
                ('kron_radius', 'kron_radius'),
                ('kron_radius_err', 'kron_radius_err'),
                ('gaia_id', 'gaia_id'),
                ('gaia_match_quality', 'gaia_match_quality'),
                ('flux_y_templfit', 'flux_y_templfit'),
                ('flux_h_templfit', 'flux_h_templfit'),
                ('flux_j_templfit', 'flux_j_templfit'),
                ('flux_u_ext_decam_templfit', 'flux_u_ext_decam_templfit'), # south; inapplicapable
                ('flux_g_ext_decam_templfit', 'flux_g_ext_decam_templfit'), # south
                ('flux_r_ext_decam_templfit', 'flux_r_ext_decam_templfit'), # south
                ('flux_i_ext_decam_templfit', 'flux_i_ext_decam_templfit'), # south
                ('flux_z_ext_decam_templfit', 'flux_z_ext_decam_templfit'), # south
                ('flux_u_ext_megacam_templfit', 'flux_u_ext_megacam_templfit'), # north
                ('flux_r_ext_megacam_templfit', 'flux_r_ext_megacam_templfit'), # north
                ('flux_g_ext_jpcam_templfit', 'flux_g_ext_jpcam_templfit'), # north; JEDIS; not available yet
                ('flux_i_ext_panstarrs_templfit', 'flux_i_ext_panstarrs_templfit'), # north
                ('flux_z_ext_panstarrs_templfit', 'flux_z_ext_panstarrs_templfit'), # north; in DR2
                ('flux_g_ext_hsc_templfit', 'flux_g_ext_hsc_templfit'), # north
                ('flux_z_ext_hsc_templfit', 'flux_z_ext_hsc_templfit'), # north
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

            if want_oid and (has_oid or has_oid_alt or has_source_id):
                # Equality join on object_id
                if has_oid:
                    user_id_col = 'object_id'
                elif has_oid_alt:
                    user_id_col = 'object_id_euclid'
                else:
                    user_id_col = 'source_id'

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
                if use_object_id is True and not (has_oid or has_oid_alt or has_source_id):
                    logger.warning("use_object_id=True but no object_id/object_id_euclid/source_id found; falling back to spatial match")
                # Spatial crossmatch via distance predicate
                radius_deg = radius / 3600.0  # Convert to degrees
                query = f"""
                SELECT u.*, 
                       m.object_id, 
                       m.right_ascension AS mer_ra, 
                       m.declination AS mer_dec,
                       m.mu_max AS mu_max,
                       m.mumax_minus_mag AS mumax_minus_mag,
                       m.point_like_prob AS point_like_prob,
                       m.extended_prob AS extended_prob,
                       m.kron_radius AS kron_radius,
                       m.kron_radius_err AS kron_radius_err,
                       m.gaia_id AS gaia_id,
                       m.gaia_match_quality AS gaia_match_quality, 
                       flux_vis_psf, fluxerr_vis_psf,
                       flux_y_templfit, fluxerr_y_templfit,
                       flux_j_templfit, fluxerr_j_templfit,
                       flux_h_templfit, fluxerr_h_templfit, 
                       flux_g_ext_decam_templfit, fluxerr_g_ext_decam_templfit, 
                       flux_r_ext_decam_templfit, fluxerr_r_ext_decam_templfit,
                       flux_i_ext_decam_templfit, fluxerr_i_ext_decam_templfit,
                       flux_z_ext_decam_templfit, fluxerr_z_ext_decam_templfit,
                       flux_u_ext_megacam_templfit, fluxerr_u_ext_megacam_templfit,
                       flux_r_ext_megacam_templfit, fluxerr_r_ext_megacam_templfit,
                       flux_i_ext_panstarrs_templfit, fluxerr_i_ext_panstarrs_templfit,
                       flux_g_ext_hsc_templfit, fluxerr_g_ext_hsc_templfit,
                       flux_z_ext_hsc_templfit, fluxerr_z_ext_hsc_templfit,
                       m.segmentation_map_id, m.segmentation_area
                FROM TAP_UPLOAD.{upload_name} AS u
                JOIN {mer_table} AS m
                    ON DISTANCE(u.{ra_col}, u.{dec_col}, m.right_ascension, m.declination) < {radius_deg}
                ORDER BY u.{ra_col}
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
            
            if not fetch_results:
                return {
                    'job': job,
                    'query': query.strip(),
                    'upload_table_name': upload_name,
                    'row_count': len(batch),
                }

            result = job.get_results()
            
            # Convert separation or compute locally
            if len(result) > 0 and 'separation_deg' in result.colnames:
                result['separation_arcsec'] = result['separation_deg'] * 3600.0
                result.remove_column('separation_deg')
            elif len(result) > 0 and not want_oid:
                try:
                    user_coords = SkyCoord(result[ra_col] * u.deg, result[dec_col] * u.deg, frame='icrs')
                    mer_coords = SkyCoord(result['mer_ra'] * u.deg, result['mer_dec'] * u.deg, frame='icrs')
                    result['separation_arcsec'] = user_coords.separation(mer_coords).arcsec
                except Exception as exc:
                    logger.warning(f"Failed to compute separations locally: {exc}")

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

    def _build_async_job_info(
        self,
        submission: Dict[str, Any],
        mer_table: str,
        radius: float,
        ra_col: str,
        dec_col: str,
        idr_field: Optional[str],
        row_count: int,
    ) -> Dict[str, Any]:
        """Create metadata dictionary for an asynchronous TAP job."""
        job = submission.get('job')

        def _get_attr(obj, attr):
            return getattr(obj, attr, None) if obj is not None else None

        job_id = _get_attr(job, 'jobid')
        job_phase = _get_attr(job, 'phase')
        job_url = _get_attr(job, 'url')
        results_url = _get_attr(job, 'remote_results_location')

        job_info = {
            'type': 'euclidkit_crossmatch_async_job',
            'environment': self.environment,
            'submitted_at_utc': datetime.utcnow().isoformat(timespec='seconds') + 'Z',
            'job_id': job_id,
            'job_phase': job_phase,
            'job_url': job_url,
            'results_url': results_url,
            'upload_table': submission.get('upload_table_name'),
            'row_count': row_count,
            'mer_table': mer_table,
            'radius_arcsec': radius,
            'ra_column': ra_col,
            'dec_column': dec_col,
            'idr_field': idr_field.upper() if idr_field else None,
            'query': submission.get('query'),
        }
        return job_info

    def _write_job_info(self, job_info: Dict[str, Any], output_path: Union[str, Path]):
        """Write async job metadata to disk."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as fh:
            json.dump(job_info, fh, indent=2)
    
    def _query_spectra_batch(self, batch: Table, idr_field: Optional[str] = None) -> Table:
        """Query spectra for a batch of object IDs."""
        
        # Upload user table to archive for querying
        with tempfile.NamedTemporaryFile(suffix='.vot', delete=False) as tmp_file:
            batch.write(tmp_file.name, format='votable', overwrite=True)
            tmp_name = tmp_file.name
        
        try:
            # Use temporary upload with launch_job
            upload_name = f"user_spectra_batch_{np.random.randint(10000, 99999)}"
            spectra_source_table = self._get_spectra_source_table_name(idr_field=idr_field)
            
            # Construct spectra query
            query = f"""
            SELECT s.source_id, s.ra_obj, s.dec_obj, s.datalabs_path, 
                   s.file_name, s.hdu_index,
                   u.object_id
            FROM TAP_UPLOAD.{upload_name} AS u
            JOIN {spectra_source_table} AS s
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
        output_file: Optional[Union[str, Path]] = None,
        idr_field: Optional[str] = None,
    ) -> Table:
        """
        Query spectral sources from Euclid archive.
        
        Parameters
        ----------
        crossmatch_table : Table, optional
            Table with object_id column from crossmatch results
        output_file : str or Path, optional
            Output file path to save results
        idr_field : {'WIDE', 'DEEP'}, optional
            IDR field selector. Used only when environment='IDR'. Defaults to
            'WIDE' if not provided.
            
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
            
            batch_result = self._query_spectra_batch(batch, idr_field=idr_field)
            
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
        from euclidkit.core.spectra import SpectrumLoader
        
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
        from euclidkit.core.spectra import SpectrumCompiler
        
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
