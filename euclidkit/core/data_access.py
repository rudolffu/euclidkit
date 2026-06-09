"""
Data access module for euclidkit package.

Provides interfaces to Euclid science archive and data volumes.
"""

import os
import logging
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Union, Optional, List, Dict, Any, Tuple
import tempfile

import numpy as np
import pandas as pd
from astropy.table import Table, join, vstack
from astropy.coordinates import SkyCoord
from astropy import units as u

from euclidkit.config import config
from euclidkit.utils.io import load_table, save_table

logger = logging.getLogger(__name__)


def _utc_timestamp() -> str:
    """Return an ISO UTC timestamp with the existing trailing-Z format."""
    return datetime.now(timezone.utc).isoformat(timespec='seconds').replace('+00:00', 'Z')

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

    def _sanitize_column_name(self, name: str) -> str:
        """Convert a column name to an ADQL-safe simple identifier."""
        cleaned = re.sub(r"[^0-9A-Za-z_]", "_", str(name)).strip("_").lower()
        if not cleaned:
            cleaned = "col"
        if not re.match(r"^[A-Za-z_]", cleaned):
            cleaned = f"c_{cleaned}"
        return cleaned

    @staticmethod
    def _stringify_upload_value(value: Any) -> str:
        """Convert object-dtype upload values to VOTable-safe scalar strings."""
        if np.ma.is_masked(value) or value is None:
            return ""
        if isinstance(value, bytes):
            return value.decode("utf-8", errors="replace")
        if isinstance(value, (np.ndarray, list, tuple, dict)):
            try:
                return json.dumps(value.tolist() if isinstance(value, np.ndarray) else value)
            except TypeError:
                return str(value)
        try:
            if np.isscalar(value) and pd.isna(value):
                return ""
        except (TypeError, ValueError):
            pass
        return str(value)

    def _sanitize_upload_table_columns(self, table: Table) -> Tuple[Table, Dict[str, str]]:
        """
        Return a copy of a table with ADQL-safe upload column names.

        Names are normalized to lower-case and restricted to [A-Za-z_][A-Za-z0-9_]*.
        Object-dtype columns are converted to strings so Astropy can serialize
        masked/null scalar values to VOTable for TAP_UPLOAD.
        """
        renamed = table.copy()
        mapping: Dict[str, str] = {}
        used: set[str] = set()

        for original in table.colnames:
            base = self._sanitize_column_name(original)
            candidate = base
            suffix = 1
            while candidate in used:
                suffix += 1
                candidate = f"{base}_{suffix}"
            used.add(candidate)
            mapping[str(original)] = candidate

        changes = {old: new for old, new in mapping.items() if old != new}
        for old, new in changes.items():
            renamed.rename_column(old, new)

        object_columns = []
        for colname in renamed.colnames:
            column = renamed[colname]
            dtype = getattr(column, "dtype", None)
            if dtype is not None and dtype.kind == "O":
                renamed[colname] = [self._stringify_upload_value(value) for value in column]
                object_columns.append(colname)

        if changes:
            logger.info("Sanitized %d upload column name(s) for ADQL compatibility", len(changes))
            for old, new in changes.items():
                logger.info("  %s -> %s", old, new)

        if object_columns:
            logger.info(
                "Converted %d object-dtype upload column(s) to strings for VOTable compatibility: %s",
                len(object_columns),
                ", ".join(object_columns),
            )

        return renamed, mapping

    @staticmethod
    def _is_missing_table_value(value: Any) -> bool:
        """Return True for null-like scalar values, preserving falsey real values."""
        if np.ma.is_masked(value) or value is None:
            return True

        try:
            missing = pd.isna(value)
        except (TypeError, ValueError):
            return False

        if isinstance(missing, (np.ndarray, list, tuple)):
            return bool(np.all(missing))
        try:
            if hasattr(missing, 'item'):
                missing = missing.item()
            return bool(missing)
        except (TypeError, ValueError):
            return False

    def _drop_empty_table_columns(self, table: Table) -> Tuple[Table, List[str]]:
        """
        Remove columns where every row is null/missing.

        Empty strings, zero, and False are retained as real values. Zero-row
        tables are returned unchanged to avoid dropping schema by vacuous truth.
        """
        if table is None or len(table) == 0:
            return table.copy() if table is not None else Table(), []

        dropped = [
            colname
            for colname in table.colnames
            if all(self._is_missing_table_value(value) for value in table[colname])
        ]
        if not dropped:
            return table.copy(), []

        cleaned = table.copy()
        cleaned.remove_columns(dropped)
        logger.info(
            "Dropped %d entirely empty crossmatch result column(s): %s",
            len(dropped),
            ", ".join(dropped),
        )
        return cleaned, dropped

    def _finalize_crossmatch_result(
        self,
        table: Table,
        drop_empty_columns: bool,
    ) -> Tuple[Table, List[str]]:
        """Apply optional final crossmatch result cleanup before save/return."""
        if not drop_empty_columns:
            return table, []
        return self._drop_empty_table_columns(table)
    
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
        idr_deep_partition: str = 'survey',
        full_async: bool = False,
        async_chunk_size: int = 500000,
        drop_empty_columns: bool = False,
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
        idr_deep_partition : {'survey', 'mode', 'both'}, default 'survey'
            IDR DEEP MER partition to query when ``idr_field='DEEP'``.
        full_async : bool, optional
            If True, use asynchronous TAP submission mode.
            For tables larger than ``async_chunk_size``, the table is split into
            multiple async jobs and results are downloaded/merged automatically.
        async_chunk_size : int, optional
            Chunk size (rows per async job) used when ``full_async=True`` for
            large input tables. Default is 500000.
        drop_empty_columns : bool, default False
            If True, drop columns that are entirely null/missing from the final
            crossmatch result before saving and returning it.
            
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

        # Normalize column names so TAP_UPLOAD queries are robust against non-simple names.
        user_data, colmap = self._sanitize_upload_table_columns(user_data)
        ra_col = colmap.get(ra_col, ra_col)
        dec_col = colmap.get(dec_col, dec_col)
        
        if max_sources:
            user_data = user_data[:max_sources]
            logger.info(f"Limited to {max_sources} sources")
        
        # Determine MER table name(s) based on environment. IDR WIDE spans two
        # MER tables; explicit mer_table overrides keep single-table behavior.
        mer_tables = [mer_table] if mer_table is not None else self._get_mer_table_names(
            idr_field=idr_field,
            idr_deep_partition=idr_deep_partition,
        )
        use_wide_fallback = self._uses_idr_wide_fallback(mer_tables)

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

        logger.info("Using MER table(s): %s", ", ".join(mer_tables))
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
                output_path = Path(output_file)
                output_path.parent.mkdir(parents=True, exist_ok=True)
                chunk_suffix = output_path.suffix or ".fits"
                manifest_path = output_path.with_suffix(output_path.suffix + ".manifest.json")
                batches = [
                    (i, min(i + async_chunk_size, len(user_data)))
                    for i in range(0, len(user_data), async_chunk_size)
                ]
                manifest: Dict[str, Any] = {
                    'type': 'euclidkit_crossmatch_async_chunked',
                    'environment': self.environment,
                    'submitted_at_utc': _utc_timestamp(),
                    'row_count': len(user_data),
                    'chunk_size': async_chunk_size,
                    'chunk_count': len(batches),
                    'query_count': 0,
                    'mer_tables': mer_tables,
                    'wide_fallback_mode': use_wide_fallback,
                    'idr_field': idr_field.upper() if idr_field else None,
                    'idr_deep_partition': idr_deep_partition.lower(),
                    'output_file': str(output_path),
                    'chunks': [],
                }
                part_files: List[Path] = []
                query_idx = 0
                for idx, (start, end) in enumerate(batches, 1):
                    logger.info("Submitting async chunk %d/%d (%d rows)", idx, len(batches), end - start)
                    batch = user_data[start:end]
                    tables_for_chunk = [mer_tables[0]]
                    if use_wide_fallback:
                        logger.info(
                            "IDR WIDE fallback chunk %d: querying wide_survey first",
                            idx,
                        )
                    else:
                        tables_for_chunk = mer_tables

                    for table_name in tables_for_chunk:
                        query_idx += 1
                        submission = self._crossmatch_batch(
                            batch,
                            ra_col,
                            dec_col,
                            radius,
                            table_name,
                            use_object_id,
                            force_async=True,
                            fetch_results=False,
                        )
                        job = submission.get('job')
                        job_id = getattr(job, 'jobid', None)
                        part_path = output_path.with_name(
                            f"{output_path.stem}_part_{query_idx:04d}{chunk_suffix}"
                        )
                        chunk_rows = 0
                        status = 'FAILED'
                        error_msg = None
                        try:
                            batch_result = job.get_results()
                            chunk_rows = len(batch_result) if batch_result is not None else 0
                            save_table(batch_result, part_path)
                            self._delete_async_job(job_id)
                            part_files.append(part_path)
                            status = 'COMPLETED'
                            logger.info(
                                "Saved async chunk %d for %s (%d rows) to %s",
                                query_idx, table_name, chunk_rows, part_path
                            )
                        except Exception as exc:
                            error_msg = str(exc)
                            logger.warning(
                                "Failed to fetch async results for chunk %d table %s: %s",
                                idx, table_name, exc,
                            )
                            manifest['chunks'].append({
                                'index': query_idx,
                                'batch_index': idx,
                                'row_start': start,
                                'row_end': end,
                                'mer_table': table_name,
                                'job_id': job_id,
                                'status': status,
                                'rows': chunk_rows,
                                'file': str(part_path),
                                'query': submission.get('query'),
                                'error': error_msg,
                            })
                            with open(manifest_path, 'w', encoding='utf-8') as fh:
                                json.dump(manifest, fh, indent=2)
                            raise

                        manifest['chunks'].append({
                            'index': query_idx,
                            'batch_index': idx,
                            'row_start': start,
                            'row_end': end,
                            'mer_table': table_name,
                            'stage': 'wide_survey' if use_wide_fallback else 'single_table',
                            'input_rows': len(batch),
                            'unmatched_rows': None,
                            'job_id': job_id,
                            'status': status,
                            'rows': chunk_rows,
                            'file': str(part_path),
                            'query': submission.get('query'),
                            'error': error_msg,
                        })
                        with open(manifest_path, 'w', encoding='utf-8') as fh:
                            json.dump(manifest, fh, indent=2)

                    if use_wide_fallback:
                        survey_result = load_table(part_files[-1]) if part_files else Table()
                        unmatched_batch = self._filter_unmatched_input_rows(
                            batch,
                            survey_result,
                            use_object_id=want_oid,
                            ra_col=ra_col,
                            dec_col=dec_col,
                        )
                        logger.info(
                            "IDR WIDE fallback chunk %d: %d/%d rows unmatched in wide_survey",
                            idx,
                            len(unmatched_batch),
                            len(batch),
                        )
                        if len(unmatched_batch) == 0:
                            manifest['chunks'][-1]['unmatched_rows'] = 0
                            with open(manifest_path, 'w', encoding='utf-8') as fh:
                                json.dump(manifest, fh, indent=2)
                            continue

                        query_idx += 1
                        mode_table = mer_tables[1]
                        submission = self._crossmatch_batch(
                            unmatched_batch,
                            ra_col,
                            dec_col,
                            radius,
                            mode_table,
                            use_object_id,
                            force_async=True,
                            fetch_results=False,
                        )
                        job = submission.get('job')
                        job_id = getattr(job, 'jobid', None)
                        part_path = output_path.with_name(
                            f"{output_path.stem}_part_{query_idx:04d}{chunk_suffix}"
                        )
                        chunk_rows = 0
                        status = 'FAILED'
                        error_msg = None
                        try:
                            batch_result = job.get_results()
                            chunk_rows = len(batch_result) if batch_result is not None else 0
                            save_table(batch_result, part_path)
                            self._delete_async_job(job_id)
                            part_files.append(part_path)
                            status = 'COMPLETED'
                            logger.info(
                                "Saved IDR WIDE fallback chunk %d for %s (%d rows) to %s",
                                idx, mode_table, chunk_rows, part_path
                            )
                        except Exception as exc:
                            error_msg = str(exc)
                            logger.warning(
                                "Failed to fetch async fallback results for chunk %d table %s: %s",
                                idx, mode_table, exc,
                            )
                            manifest['chunks'].append({
                                'index': query_idx,
                                'batch_index': idx,
                                'row_start': start,
                                'row_end': end,
                                'mer_table': mode_table,
                                'stage': 'wide_mode_unmatched',
                                'input_rows': len(unmatched_batch),
                                'unmatched_rows': len(unmatched_batch),
                                'job_id': job_id,
                                'status': status,
                                'rows': chunk_rows,
                                'file': str(part_path),
                                'query': submission.get('query'),
                                'error': error_msg,
                            })
                            with open(manifest_path, 'w', encoding='utf-8') as fh:
                                json.dump(manifest, fh, indent=2)
                            raise

                        manifest['chunks'].append({
                            'index': query_idx,
                            'batch_index': idx,
                            'row_start': start,
                            'row_end': end,
                            'mer_table': mode_table,
                            'stage': 'wide_mode_unmatched',
                            'input_rows': len(unmatched_batch),
                            'unmatched_rows': len(unmatched_batch),
                            'job_id': job_id,
                            'status': status,
                            'rows': chunk_rows,
                            'file': str(part_path),
                            'query': submission.get('query'),
                            'error': error_msg,
                        })
                        with open(manifest_path, 'w', encoding='utf-8') as fh:
                            json.dump(manifest, fh, indent=2)

                if part_files:
                    merge_tables = [load_table(p) for p in part_files]
                    final_result = self._combine_mer_results(merge_tables)
                else:
                    final_result = Table()
                final_result, dropped_empty_columns = self._finalize_crossmatch_result(
                    final_result,
                    drop_empty_columns,
                )
                manifest['query_count'] = len(manifest['chunks'])
                manifest['dropped_empty_columns'] = dropped_empty_columns
                manifest['dropped_empty_column_count'] = len(dropped_empty_columns)
                with open(manifest_path, 'w', encoding='utf-8') as fh:
                    json.dump(manifest, fh, indent=2)

                save_table(final_result, output_path)
                logger.info(
                    "Saved merged chunked async results (%d rows) to %s",
                    len(final_result), output_path
                )
                logger.info("Chunk manifest saved to %s", manifest_path)
                return {
                    'type': 'euclidkit_crossmatch_async_chunked',
                    'environment': self.environment,
                    'submitted_at_utc': manifest['submitted_at_utc'],
                    'row_count': len(user_data),
                    'chunk_size': async_chunk_size,
                    'chunk_count': len(batches),
                    'query_count': len(manifest['chunks']),
                    'mer_tables': mer_tables,
                    'job_ids': [c.get('job_id') for c in manifest['chunks'] if c.get('job_id')],
                    'results_downloaded': True,
                    'result_row_count': len(final_result),
                    'output_file': str(output_path),
                    'manifest_file': str(manifest_path),
                    'idr_field': idr_field.upper() if idr_field else None,
                    'idr_deep_partition': idr_deep_partition.lower(),
                    'dropped_empty_columns': dropped_empty_columns,
                    'dropped_empty_column_count': len(dropped_empty_columns),
                }

            job_infos = []
            async_results = []
            download_errors = []
            successful_job_ids = []
            tables_to_query = [mer_tables[0]] if use_wide_fallback else mer_tables
            for table_name in tables_to_query:
                submission = self._crossmatch_batch(
                    user_data,
                    ra_col,
                    dec_col,
                    radius,
                    table_name,
                    use_object_id,
                    force_async=True,
                    fetch_results=False,
                )

                job_info = self._build_async_job_info(
                    submission,
                    mer_table=table_name,
                    radius=radius,
                    ra_col=ra_col,
                    dec_col=dec_col,
                    idr_field=idr_field,
                    idr_deep_partition=idr_deep_partition,
                    row_count=len(user_data),
                )
                job_infos.append(job_info)
                job = submission.get('job')
                try:
                    logger.info(
                        "Attempting to download results for async job %s against %s",
                        getattr(job, 'jobid', None),
                        table_name,
                    )
                    async_result = job.get_results()
                    if async_result is not None:
                        async_results.append(async_result)
                        successful_job_ids.append(getattr(job, 'jobid', None))
                        logger.info(
                            "Downloaded async results (%d rows) from %s",
                            len(async_result),
                            table_name,
                        )
                except Exception as exc:
                    download_errors.append(str(exc))
                    logger.warning(
                        "Async job submitted for %s but result download failed: %s",
                        table_name, exc,
                    )

            if use_wide_fallback and async_results:
                unmatched_data = self._filter_unmatched_input_rows(
                    user_data,
                    async_results[0],
                    use_object_id=want_oid,
                    ra_col=ra_col,
                    dec_col=dec_col,
                )
                logger.info(
                    "IDR WIDE fallback async mode: %d/%d rows unmatched in wide_survey",
                    len(unmatched_data),
                    len(user_data),
                )
                if len(unmatched_data) > 0:
                    table_name = mer_tables[1]
                    submission = self._crossmatch_batch(
                        unmatched_data,
                        ra_col,
                        dec_col,
                        radius,
                        table_name,
                        use_object_id,
                        force_async=True,
                        fetch_results=False,
                    )
                    job_info = self._build_async_job_info(
                        submission,
                        mer_table=table_name,
                        radius=radius,
                        ra_col=ra_col,
                        dec_col=dec_col,
                        idr_field=idr_field,
                        idr_deep_partition=idr_deep_partition,
                        row_count=len(unmatched_data),
                    )
                    job_infos.append(job_info)
                    job = submission.get('job')
                    try:
                        logger.info(
                            "Attempting to download IDR WIDE fallback async job %s against %s",
                            getattr(job, 'jobid', None),
                            table_name,
                        )
                        async_result = job.get_results()
                        if async_result is not None:
                            async_results.append(async_result)
                            successful_job_ids.append(getattr(job, 'jobid', None))
                            logger.info(
                                "Downloaded fallback async results (%d rows) from %s",
                                len(async_result),
                                table_name,
                            )
                    except Exception as exc:
                        download_errors.append(str(exc))
                        logger.warning(
                            "Async fallback job submitted for %s but result download failed: %s",
                            table_name, exc,
                        )

            final_result = self._combine_mer_results(async_results)
            final_result, dropped_empty_columns = self._finalize_crossmatch_result(
                final_result,
                drop_empty_columns and len(download_errors) == 0,
            )
            results_downloaded = len(download_errors) == 0
            result_row_count = len(final_result)
            combined_job_info = {
                'type': 'euclidkit_crossmatch_async_job',
                'environment': self.environment,
                'submitted_at_utc': _utc_timestamp(),
                'job_id': job_infos[0].get('job_id') if len(job_infos) == 1 else None,
                'job_ids': [info.get('job_id') for info in job_infos if info.get('job_id')],
                'row_count': len(user_data),
                'mer_table': mer_tables[0] if len(mer_tables) == 1 else None,
                'mer_tables': mer_tables,
                'radius_arcsec': radius,
                'ra_column': ra_col,
                'dec_column': dec_col,
                'idr_field': idr_field.upper() if idr_field else None,
                'idr_deep_partition': idr_deep_partition.lower(),
                'query': job_infos[0].get('query') if len(job_infos) == 1 else None,
                'queries': [info.get('query') for info in job_infos],
                'results_downloaded': results_downloaded,
                'result_row_count': result_row_count,
                'download_error': '; '.join(download_errors) if download_errors else None,
                'dropped_empty_columns': dropped_empty_columns,
                'dropped_empty_column_count': len(dropped_empty_columns),
            }

            if results_downloaded:
                save_table(final_result, output_file)
                for job_id in successful_job_ids:
                    self._delete_async_job(job_id)
                logger.info(
                    "Downloaded and merged async results (%d rows) to %s",
                    result_row_count,
                    output_file,
                )
            else:
                output_path = Path(output_file)
                if output_path.suffix.lower() == '.json':
                    job_info_path = output_path
                else:
                    job_info_path = output_path.with_name(output_path.name + '.job.json')
                self._write_job_info(combined_job_info, job_info_path)
                combined_job_info['job_info_file'] = str(job_info_path)
                logger.info(f"Saved async job info to {job_info_path}")

            return combined_job_info

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
            
            tables_for_batch = [mer_tables[0]] if use_wide_fallback else mer_tables
            survey_result = Table()
            for table_name in tables_for_batch:
                batch_result = self._crossmatch_batch(
                    batch,
                    ra_col,
                    dec_col,
                    radius,
                    table_name,
                    use_object_id,
                    force_async=use_async_batches,
                )
                if use_wide_fallback:
                    survey_result = batch_result
                
                if len(batch_result) > 0:
                    crossmatch_results.append(batch_result)

            if use_wide_fallback:
                unmatched_batch = self._filter_unmatched_input_rows(
                    batch,
                    survey_result,
                    use_object_id=want_oid,
                    ra_col=ra_col,
                    dec_col=dec_col,
                )
                logger.info(
                    "IDR WIDE fallback batch %d/%d: %d/%d rows unmatched in wide_survey",
                    idx,
                    total_batches,
                    len(unmatched_batch),
                    len(batch),
                )
                if len(unmatched_batch) > 0:
                    mode_result = self._crossmatch_batch(
                        unmatched_batch,
                        ra_col,
                        dec_col,
                        radius,
                        mer_tables[1],
                        use_object_id,
                        force_async=use_async_batches,
                    )
                    if len(mode_result) > 0:
                        crossmatch_results.append(mode_result)
        
        if crossmatch_results:
            final_result = self._combine_mer_results(crossmatch_results)
        else:
            logger.warning("No crossmatches found")
            return Table()
        final_result, _ = self._finalize_crossmatch_result(
            final_result,
            drop_empty_columns,
        )
        
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
            Resource to upload. Input is normalized to ADQL-safe column names
            before upload.
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
        temp_upload_path: Optional[Path] = None

        if isinstance(table, (str, Path)):
            path = Path(table)
            if not path.exists():
                raise FileNotFoundError(f"Table file not found: {path}")
            loaded = load_table(path)
            sanitized, _ = self._sanitize_upload_table_columns(loaded)
            with tempfile.NamedTemporaryFile(suffix='.vot', delete=False) as tmp_file:
                sanitized.write(tmp_file.name, format='votable', overwrite=True)
                temp_upload_path = Path(tmp_file.name)
            resource = str(temp_upload_path)
            inferred_format = 'votable'
        else:
            # astropy Table upload
            resource, _ = self._sanitize_upload_table_columns(table)
            if inferred_format is None:
                inferred_format = 'votable'

        if overwrite:
            try:
                self.euclid.delete_user_table(table_name=table_name, force_removal=True, verbose=verbose)
                logger.info(f"Removed existing table '{table_name}' prior to upload")
            except Exception as exc:  # pragma: no cover - network errors
                logger.warning(f"Failed to delete existing table '{table_name}': {exc}")

        try:
            job = self.euclid.upload_table(
                upload_resource=resource,
                table_name=table_name,
                table_description=description,
                format=inferred_format,
                verbose=verbose,
            )
        finally:
            if temp_upload_path is not None and temp_upload_path.exists():
                temp_upload_path.unlink()

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
        tap_queries = [
            (
                "SELECT column_name "
                "FROM TAP_SCHEMA.columns "
                f"WHERE lower(schema_name) = lower('{schema_name}') "
                f"AND lower(table_name) = lower('{table_name}')"
            ),
            (
                "SELECT column_name "
                "FROM TAP_SCHEMA.columns "
                f"WHERE lower(table_name) = lower('{qualified_table_name}')"
            ),
        ]

        for query in tap_queries:
            try:
                job = self.euclid.launch_job(query)
                result = job.get_results()
                if result is not None and len(result) > 0 and 'column_name' in result.colnames:
                    cols = [str(name).lower() for name in result['column_name']]
                    if cols:
                        return cols
            except Exception as exc:
                logger.debug("TAP_SCHEMA column lookup failed for %s: %s", qualified_table_name, exc)

        # Final fallback: query the table directly and infer columns from result schema.
        # This handles archive setups where TAP_SCHEMA metadata for user tables is incomplete.
        fallback_query = f"SELECT TOP 1 * FROM {qualified_table_name}"
        try:
            job = self.euclid.launch_job(fallback_query)
            result = job.get_results()
            if result is not None:
                cols = [str(name).lower() for name in result.colnames]
                if cols:
                    logger.info(
                        "Resolved columns for %s using direct table introspection fallback",
                        qualified_table_name,
                    )
                    return cols
        except Exception as exc:
            logger.debug(
                "Direct table introspection failed for %s with query `%s`: %s",
                qualified_table_name,
                fallback_query,
                exc,
            )

        return []

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
        idr_deep_partition: str = 'survey',
        full_async: bool = False,
        drop_empty_columns: bool = False,
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

        mer_tables = [mer_table] if mer_table is not None else self._get_mer_table_names(
            idr_field=idr_field,
            idr_deep_partition=idr_deep_partition,
        )
        use_wide_fallback = self._uses_idr_wide_fallback(mer_tables)

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

        logger.info("Using MER table(s): %s", ", ".join(mer_tables))
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
            "m.det_quality_flag AS det_quality_flag",
            "m.parent_id AS parent_id",
            "m.spurious_flag AS spurious_flag",
            "m.vis_det AS vis_det",
            "m.flag_vis AS flag_vis",
            "m.flag_y AS flag_y",
            "m.flag_j AS flag_j",
            "m.flag_h AS flag_h",
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
            order_col = user_id_col
        else:
            radius_deg = radius / 3600.0
            order_col = ra_col

        # User-table crossmatch always runs via async TAP execution.
        # Preflight row-count check is mandatory to decide whether to chunk.
        short_table_name = user_table_ref.split('.', 1)[1]
        oid_col = f"{short_table_name}_oid"
        if oid_col.lower() not in user_cols:
            raise ValueError(
                f"Expected indexed OID column '{oid_col}' in {user_table_ref}, "
                f"but it was not found."
            )

        preflight_query = (
            f"SELECT MIN(u.{oid_col}) AS min_oid, "
            f"MAX(u.{oid_col}) AS max_oid, "
            "COUNT(*) AS n_rows "
            f"FROM {user_table_ref} AS u"
        )
        preflight_job = self.euclid.launch_job(preflight_query)
        preflight_result = preflight_job.get_results()
        if preflight_result is None or len(preflight_result) == 0:
            raise RuntimeError(f"Could not determine row count for remote table {user_table_ref}")

        def _to_int(value: Any) -> Optional[int]:
            if value is None:
                return None
            if hasattr(value, 'item'):
                value = value.item()
            try:
                return int(value)
            except Exception:
                return None

        table_rows = _to_int(preflight_result['n_rows'][0]) or 0
        min_oid = _to_int(preflight_result['min_oid'][0])
        max_oid = _to_int(preflight_result['max_oid'][0])
        effective_rows = min(table_rows, int(max_sources)) if max_sources is not None else table_rows

        logger.info(
            "User table preflight: n_rows=%d, min_oid=%s, max_oid=%s",
            table_rows,
            str(min_oid),
            str(max_oid),
        )

        chunk_threshold = 2_000_000
        chunk_size = 500_000
        use_chunked = effective_rows >= chunk_threshold
        logger.info(
            "User-table async mode decision: effective_rows=%d, threshold=%d -> %s",
            effective_rows,
            chunk_threshold,
            "chunked async" if use_chunked else "single async query",
        )

        def _build_remote_query(
            table_name: str,
            row_filter: Optional[str] = None,
            top_limit: Optional[int] = None,
            exclude_survey_matches: bool = False,
        ) -> str:
            top_clause = f"TOP {int(top_limit)} " if top_limit is not None else ""
            filters = []
            if row_filter:
                filters.append(f"({row_filter})")
            if exclude_survey_matches:
                filters.append(
                    self._wide_fallback_exclusion_filter(
                        mer_tables[0],
                        want_oid=want_oid,
                        user_id_col=user_id_col if want_oid else None,
                        ra_col=ra_col if not want_oid else None,
                        dec_col=dec_col if not want_oid else None,
                        radius_deg=radius_deg if not want_oid else None,
                    )
                )
            where_clause = f"WHERE {' AND '.join(filters)} " if filters else ""
            if want_oid:
                return (
                    f"SELECT {top_clause}u.*, {', '.join(mer_select_parts)} "
                    f"FROM {user_table_ref} AS u "
                    f"JOIN {table_name} AS m ON u.{user_id_col} = m.object_id "
                    f"{where_clause}"
                    f"ORDER BY u.{order_col}"
                )
            return (
                f"SELECT {top_clause}u.*, {', '.join(mer_select_parts)} "
                f"FROM {user_table_ref} AS u "
                f"JOIN {table_name} AS m "
                f"ON DISTANCE(u.{ra_col}, u.{dec_col}, m.right_ascension, m.declination) < {radius_deg} "
                f"{where_clause}"
                f"ORDER BY u.{order_col}"
            )

        if use_chunked:
            if output_file is None:
                raise ValueError(
                    "output_file must be provided for automatic chunked user-table crossmatch."
                )
            if min_oid is None or max_oid is None:
                raise RuntimeError(
                    f"Could not determine OID range ({oid_col}) for chunked query on {user_table_ref}."
                )

            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            chunk_suffix = output_path.suffix or ".fits"
            manifest_path = output_path.with_suffix(output_path.suffix + ".manifest.json")

            oid_ranges = [
                (start_oid, min(start_oid + chunk_size, max_oid + 1))
                for start_oid in range(min_oid, max_oid + 1, chunk_size)
            ]
            logger.info(
                "Chunked user-table async mode enabled: %d chunks, chunk_size=%d (OID range windows)",
                len(oid_ranges),
                chunk_size,
            )

            manifest: Dict[str, Any] = {
                'type': 'euclidkit_crossmatch_remote_async_chunked',
                'environment': self.environment,
                'submitted_at_utc': _utc_timestamp(),
                'user_table': user_table_ref,
                'mer_table': mer_tables[0] if len(mer_tables) == 1 else None,
                'mer_tables': mer_tables,
                'wide_fallback_mode': use_wide_fallback,
                'idr_field': idr_field.upper() if idr_field else None,
                'idr_deep_partition': idr_deep_partition.lower(),
                'oid_column': oid_col,
                'n_rows_detected': table_rows,
                'n_rows_effective': effective_rows,
                'chunk_threshold': chunk_threshold,
                'chunk_size': chunk_size,
                'max_sources': int(max_sources) if max_sources is not None else None,
                'min_oid': min_oid,
                'max_oid': max_oid,
                'chunks': [],
            }

            part_files: List[Path] = []
            query_idx = 0

            for idx, (start_oid, end_oid) in enumerate(oid_ranges, 1):
                row_filter = f"u.{oid_col} >= {start_oid} AND u.{oid_col} < {end_oid}"
                top_limit = int(max_sources) if max_sources is not None else None
                for table_index, table_name in enumerate(mer_tables):
                    exclude_survey_matches = use_wide_fallback and table_index == 1
                    query_idx += 1
                    part_path = output_path.with_name(
                        f"{output_path.stem}_part_{query_idx:04d}{chunk_suffix}"
                    )
                    query = _build_remote_query(
                        table_name,
                        row_filter=row_filter,
                        top_limit=top_limit,
                        exclude_survey_matches=exclude_survey_matches,
                    )

                    logger.info(
                        "Running query %d/%d: oid_range=[%d,%d), table=%s, top_limit=%s",
                        query_idx,
                        len(oid_ranges) * len(mer_tables),
                        start_oid,
                        end_oid,
                        table_name,
                        str(top_limit),
                    )

                    job_id = None
                    chunk_rows = 0
                    status = 'FAILED'
                    error_msg = None
                    try:
                        job = self.euclid.launch_job_async(query)
                        job_id = getattr(job, 'jobid', None)
                        chunk_result = job.get_results()
                        chunk_rows = len(chunk_result) if chunk_result is not None else 0
                        save_table(chunk_result, part_path)
                        self._delete_async_job(job_id)
                        part_files.append(part_path)
                        status = 'COMPLETED'
                        logger.info(
                            "Chunk query %d completed: rows=%d, output=%s",
                            query_idx,
                            chunk_rows,
                            part_path,
                        )
                    except Exception as exc:
                        error_msg = str(exc)
                        logger.error("Chunk query %d failed: %s", query_idx, error_msg)
                        manifest['chunks'].append({
                            'index': query_idx,
                            'oid_chunk_index': idx,
                            'oid_start': start_oid,
                            'oid_end': end_oid,
                            'mer_table': table_name,
                            'stage': 'wide_mode_unmatched' if exclude_survey_matches else (
                                'wide_survey' if use_wide_fallback else 'single_table'
                            ),
                            'query': query,
                            'job_id': job_id,
                            'status': status,
                            'rows': chunk_rows,
                            'file': str(part_path),
                            'error': error_msg,
                        })
                        with open(manifest_path, 'w', encoding='utf-8') as fh:
                            json.dump(manifest, fh, indent=2)
                        raise

                    manifest['chunks'].append({
                        'index': query_idx,
                        'oid_chunk_index': idx,
                        'oid_start': start_oid,
                        'oid_end': end_oid,
                        'mer_table': table_name,
                        'stage': 'wide_mode_unmatched' if exclude_survey_matches else (
                            'wide_survey' if use_wide_fallback else 'single_table'
                        ),
                        'query': query,
                        'job_id': job_id,
                        'status': status,
                        'rows': chunk_rows,
                        'file': str(part_path),
                        'error': error_msg,
                    })
                    with open(manifest_path, 'w', encoding='utf-8') as fh:
                        json.dump(manifest, fh, indent=2)

            if part_files:
                merge_tables = [load_table(p) for p in part_files]
                final_result = self._combine_mer_results(merge_tables)
            else:
                final_result = Table()
            final_result, dropped_empty_columns = self._finalize_crossmatch_result(
                final_result,
                drop_empty_columns,
            )
            manifest['dropped_empty_columns'] = dropped_empty_columns
            manifest['dropped_empty_column_count'] = len(dropped_empty_columns)
            with open(manifest_path, 'w', encoding='utf-8') as fh:
                json.dump(manifest, fh, indent=2)

            save_table(final_result, output_path)
            logger.info(
                "Merged %d chunk files into %s (final rows=%d)",
                len(manifest['chunks']),
                output_path,
                len(final_result),
            )
            logger.info("Chunk manifest saved to %s", manifest_path)
            logger.info("Chunked user-table crossmatch final merged rows: %d", len(final_result))

            if full_async:
                return {
                    'type': 'euclidkit_crossmatch_remote_async_chunked',
                    'environment': self.environment,
                    'results_downloaded': True,
                    'result_row_count': len(final_result),
                    'chunk_count': len(manifest['chunks']),
                    'chunk_size': chunk_size,
                    'mer_tables': mer_tables,
                    'idr_field': idr_field.upper() if idr_field else None,
                    'idr_deep_partition': idr_deep_partition.lower(),
                    'output_file': str(output_path),
                    'manifest_file': str(manifest_path),
                    'job_ids': [c.get('job_id') for c in manifest['chunks'] if c.get('job_id')],
                    'dropped_empty_columns': dropped_empty_columns,
                    'dropped_empty_column_count': len(dropped_empty_columns),
                }
            return final_result

        results = []
        job_ids = []
        queries = []
        successful_job_ids = []
        for table_index, table_name in enumerate(mer_tables):
            exclude_survey_matches = use_wide_fallback and table_index == 1
            query = _build_remote_query(
                table_name,
                top_limit=max_sources,
                exclude_survey_matches=exclude_survey_matches,
            )
            queries.append(query)
            logger.info("Running single async query for user-table crossmatch against %s", table_name)
            job = self.euclid.launch_job_async(query)
            result = job.get_results()
            job_id = getattr(job, 'jobid', None)
            if job_id:
                job_ids.append(job_id)
                successful_job_ids.append(job_id)
            if result is not None and len(result) > 0:
                results.append(result)

        result = self._combine_mer_results(results)
        result, dropped_empty_columns = self._finalize_crossmatch_result(
            result,
            drop_empty_columns,
        )

        if output_file:
            save_table(result, output_file)
            logger.info("Saved async crossmatch results to %s", output_file)
        for job_id in successful_job_ids:
            self._delete_async_job(job_id)

        if full_async:
            return {
                'type': 'euclidkit_crossmatch_remote_async',
                'environment': self.environment,
                'job_id': job_ids[0] if len(job_ids) == 1 else None,
                'job_ids': job_ids,
                'mer_tables': mer_tables,
                'idr_field': idr_field.upper() if idr_field else None,
                'idr_deep_partition': idr_deep_partition.lower(),
                'results_downloaded': True,
                'result_row_count': len(result) if result is not None else 0,
                'output_file': str(output_file) if output_file else None,
                'query': queries[0] if len(queries) == 1 else None,
                'queries': queries,
                'dropped_empty_columns': dropped_empty_columns,
                'dropped_empty_column_count': len(dropped_empty_columns),
            }

        return result
    
    def _get_mer_table_name(
        self,
        idr_field: Optional[str] = None,
        idr_deep_partition: str = 'survey',
    ) -> str:
        """Get primary MER catalogue table name for current environment."""
        return self._get_mer_table_names(
            idr_field=idr_field,
            idr_deep_partition=idr_deep_partition,
        )[0]

    def _get_mer_table_names(
        self,
        idr_field: Optional[str] = None,
        idr_deep_partition: str = 'survey',
    ) -> List[str]:
        """Get MER catalogue table names for current environment."""
        if self.environment == 'IDR':
            field = (idr_field or 'WIDE').upper()
            valid_fields = {'WIDE', 'DEEP'}
            if field not in valid_fields:
                raise ValueError(
                    f"Invalid IDR field '{idr_field}'. Expected one of {sorted(valid_fields)}"
                )
            if field == 'WIDE':
                return [
                    'catalogue.mer_catalogue_wide_survey',
                    'catalogue.mer_catalogue_wide_mode',
                ]
            partition = (idr_deep_partition or 'survey').lower()
            valid_partitions = {'survey', 'mode', 'both'}
            if partition not in valid_partitions:
                raise ValueError(
                    f"Invalid IDR DEEP partition '{idr_deep_partition}'. "
                    f"Expected one of {sorted(valid_partitions)}"
                )
            if partition == 'survey':
                return ['catalogue.mer_catalogue_deep_survey']
            if partition == 'mode':
                return ['catalogue.mer_catalogue_deep_mode']
            return [
                'catalogue.mer_catalogue_deep_survey',
                'catalogue.mer_catalogue_deep_mode',
            ]
        
        table_names = {
            'PDR': 'catalogue.mer_catalogue',
            'OTF': 'catalogue.mer_catalogue',
            'REG': 'catalogue.mer_final_catalog_fits_file_regreproc1_r2'
        }
        return [table_names.get(self.environment, 'catalogue.mer_catalogue')]

    @staticmethod
    def _uses_idr_wide_fallback(mer_tables: List[str]) -> bool:
        """Return True when MER tables represent the IDR WIDE survey/mode fallback pair."""
        return mer_tables == [
            'catalogue.mer_catalogue_wide_survey',
            'catalogue.mer_catalogue_wide_mode',
        ]

    @staticmethod
    def _value_key(value: Any) -> str:
        """Normalize scalar table values for robust membership comparisons."""
        if hasattr(value, 'item'):
            try:
                value = value.item()
            except Exception:
                pass
        return str(value)

    def _input_result_match_columns(
        self,
        input_cols: List[str],
        use_object_id: bool,
        ra_col: str,
        dec_col: str,
    ) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
        """Resolve input/result columns used to identify already matched rows."""
        if use_object_id:
            if 'object_id' in input_cols:
                return 'object_id', 'object_id_user', None, None
            if 'object_id_euclid' in input_cols:
                return 'object_id_euclid', 'object_id_euclid', None, None
            if 'source_id' in input_cols:
                return 'source_id', 'source_id', None, None
            return None, None, None, None
        return None, None, ra_col, dec_col

    def _filter_unmatched_input_rows(
        self,
        input_table: Table,
        result_table: Table,
        use_object_id: bool,
        ra_col: str,
        dec_col: str,
    ) -> Table:
        """Return input rows not represented in a MER result table."""
        if result_table is None or len(result_table) == 0:
            return input_table.copy()

        input_id_col, result_id_col, input_ra_col, input_dec_col = self._input_result_match_columns(
            list(input_table.colnames),
            use_object_id,
            ra_col,
            dec_col,
        )

        if input_id_col and result_id_col and result_id_col in result_table.colnames:
            matched = {self._value_key(value) for value in result_table[result_id_col]}
            mask = [self._value_key(value) not in matched for value in input_table[input_id_col]]
            return input_table[mask]

        if (
            input_ra_col
            and input_dec_col
            and input_ra_col in input_table.colnames
            and input_dec_col in input_table.colnames
            and input_ra_col in result_table.colnames
            and input_dec_col in result_table.colnames
        ):
            matched = {
                (self._value_key(row[input_ra_col]), self._value_key(row[input_dec_col]))
                for row in result_table
            }
            mask = [
                (
                    self._value_key(row[input_ra_col]),
                    self._value_key(row[input_dec_col]),
                ) not in matched
                for row in input_table
            ]
            return input_table[mask]

        logger.warning(
            "Could not identify matched input rows from wide_survey result; "
            "falling back to querying wide_mode with the full batch."
        )
        return input_table.copy()

    def _wide_fallback_exclusion_filter(
        self,
        survey_table: str,
        want_oid: bool,
        user_id_col: Optional[str] = None,
        ra_col: Optional[str] = None,
        dec_col: Optional[str] = None,
        radius_deg: Optional[float] = None,
    ) -> str:
        """Build a server-side anti-match predicate against the primary WIDE survey table."""
        if want_oid:
            if user_id_col is None:
                raise ValueError("user_id_col is required for object-id fallback exclusion")
            return (
                "NOT EXISTS ("
                f"SELECT 1 FROM {survey_table} AS s "
                f"WHERE s.object_id = u.{user_id_col}"
                ")"
            )

        if ra_col is None or dec_col is None or radius_deg is None:
            raise ValueError("ra_col, dec_col, and radius_deg are required for spatial fallback exclusion")
        return (
            "NOT EXISTS ("
            f"SELECT 1 FROM {survey_table} AS s "
            f"WHERE DISTANCE(u.{ra_col}, u.{dec_col}, s.right_ascension, s.declination) < {radius_deg}"
            ")"
        )

    @staticmethod
    def _combine_mer_results(tables: List[Table]) -> Table:
        """Combine MER query results, keeping the first row for duplicate object IDs."""
        non_empty = [table for table in tables if table is not None and len(table) > 0]
        if not non_empty:
            return Table()

        combined = non_empty[0] if len(non_empty) == 1 else vstack(
            non_empty,
            join_type='outer',
            metadata_conflicts='silent',
        )
        id_col = 'object_id' if 'object_id' in combined.colnames else None
        if id_col is None and 'mer_object_id' in combined.colnames:
            id_col = 'mer_object_id'
        if id_col is None:
            return combined

        seen = set()
        keep_indices = []
        for idx, value in enumerate(combined[id_col]):
            key = str(value)
            if key in seen:
                continue
            seen.add(key)
            keep_indices.append(idx)
        return combined[keep_indices] if keep_indices else combined[:0]

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
        batch, colmap = self._sanitize_upload_table_columns(batch)
        ra_col = colmap.get(ra_col, ra_col)
        dec_col = colmap.get(dec_col, dec_col)

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
                ('det_quality_flag', 'det_quality_flag'),
                ('parent_id', 'parent_id'),
                ('spurious_flag', 'spurious_flag'),
                ('vis_det', 'vis_det'),
                ('flag_vis', 'flag_vis'),
                ('flag_y', 'flag_y'),
                ('flag_j', 'flag_j'),
                ('flag_h', 'flag_h'),
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
                       m.det_quality_flag AS det_quality_flag,
                       m.parent_id AS parent_id,
                       m.spurious_flag AS spurious_flag,
                       m.vis_det AS vis_det,
                       m.flag_vis AS flag_vis,
                       m.flag_y AS flag_y,
                       m.flag_j AS flag_j,
                       m.flag_h AS flag_h,
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
        idr_deep_partition: str,
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
            'submitted_at_utc': _utc_timestamp(),
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
            'idr_deep_partition': idr_deep_partition.lower(),
            'query': submission.get('query'),
        }
        return job_info

    def _write_job_info(self, job_info: Dict[str, Any], output_path: Union[str, Path]):
        """Write async job metadata to disk."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as fh:
            json.dump(job_info, fh, indent=2)

    def _delete_async_job(self, job_id: Optional[str]) -> bool:
        """Best-effort deletion of an async TAP job by ID."""
        if not job_id:
            return False
        try:
            self.euclid.remove_jobs([job_id])
            logger.info("Deleted async TAP job: %s", job_id)
            return True
        except Exception as exc:
            logger.warning("Failed to delete async TAP job %s: %s", job_id, exc)
            return False
    
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

    def _get_zspe_candidate_table_names(
        self,
        object_type: str = "qso",
        idr_field: Optional[str] = None,
    ) -> List[Tuple[str, str]]:
        """Resolve IDR SPE redshift candidate tables and provenance labels."""
        if self.environment != "IDR":
            raise ValueError("query-zspe is only supported for the IDR environment")

        normalized_type = object_type.lower()
        if normalized_type not in {"qso", "galaxy"}:
            raise ValueError("object_type must be one of: qso, galaxy")

        field = (idr_field or "WIDE").upper()
        if field not in {"WIDE", "DEEP"}:
            raise ValueError("idr_field must be one of: WIDE, DEEP")

        prefix = f"catalogue.spectro_z_spe_{normalized_type}_candidates"
        if field == "WIDE":
            return [
                (f"{prefix}_wide_survey", "wide_survey"),
                (f"{prefix}_wide", "wide"),
            ]
        return [(f"{prefix}_deep", "deep")]

    def _query_zspe_batch(self, batch: Table, table_name: str, source_label: str) -> Table:
        """Query one SPE redshift candidate table for a batch of object IDs."""
        batch, _ = self._sanitize_upload_table_columns(batch)

        with tempfile.NamedTemporaryFile(suffix=".vot", delete=False) as tmp_file:
            batch.write(tmp_file.name, format="votable", overwrite=True)
            tmp_name = tmp_file.name

        try:
            upload_name = f"user_zspe_batch_{np.random.randint(10000, 99999)}"
            query = f"""
            SELECT
                u.object_id,
                z.spe_rank,
                z.spe_z,
                z.spe_z_err,
                '{source_label}' AS source_table
            FROM TAP_UPLOAD.{upload_name} AS u
            JOIN {table_name} AS z
                ON u.object_id = z.object_id
            """

            if len(batch) < 2000:
                job = self.euclid.launch_job(
                    query,
                    upload_resource=tmp_name,
                    upload_table_name=upload_name,
                )
            else:
                job = self.euclid.launch_job_async(
                    query,
                    upload_resource=tmp_name,
                    upload_table_name=upload_name,
                )
            result = job.get_results()
            return result if result is not None else Table()
        except Exception as e:
            logger.error(f"Error querying SPE redshift batch: {e}")
            return Table()
        finally:
            if os.path.exists(tmp_name):
                os.unlink(tmp_name)

    def _build_zspe_async_query(
        self,
        upload_name: str,
        table_names: List[Tuple[str, str]],
    ) -> str:
        """Build a server-side SPE redshift query, including WIDE fallback when needed."""
        first_table, first_label = table_names[0]
        first_select = f"""
            SELECT
                u.object_id,
                z.spe_rank,
                z.spe_z,
                z.spe_z_err,
                '{first_label}' AS source_table
            FROM TAP_UPLOAD.{upload_name} AS u
            JOIN {first_table} AS z
                ON u.object_id = z.object_id
        """

        if len(table_names) == 1:
            return first_select

        fallback_table, fallback_label = table_names[1]
        fallback_select = f"""
            SELECT
                u.object_id,
                z.spe_rank,
                z.spe_z,
                z.spe_z_err,
                '{fallback_label}' AS source_table
            FROM TAP_UPLOAD.{upload_name} AS u
            JOIN {fallback_table} AS z
                ON u.object_id = z.object_id
            WHERE NOT EXISTS (
                SELECT 1
                FROM {first_table} AS s
                WHERE s.object_id = u.object_id
            )
        """
        return f"{first_select}\nUNION ALL\n{fallback_select}"

    def _submit_zspe_async_batch(
        self,
        batch: Table,
        table_names: List[Tuple[str, str]],
    ) -> Dict[str, Any]:
        """Submit one SPE redshift async query for a batch of object IDs."""
        batch, _ = self._sanitize_upload_table_columns(batch)

        with tempfile.NamedTemporaryFile(suffix=".vot", delete=False) as tmp_file:
            batch.write(tmp_file.name, format="votable", overwrite=True)
            tmp_name = tmp_file.name

        try:
            upload_name = f"user_zspe_async_{np.random.randint(10000, 99999)}"
            query = self._build_zspe_async_query(upload_name, table_names)
            job = self.euclid.launch_job_async(
                query,
                upload_resource=tmp_name,
                upload_table_name=upload_name,
            )
            return {
                "job": job,
                "query": query.strip(),
                "upload_table_name": upload_name,
                "row_count": len(batch),
            }
        finally:
            if os.path.exists(tmp_name):
                os.unlink(tmp_name)

    def _join_zspe_with_input(self, input_table: Table, zspe_table: Table) -> Table:
        """Inner join input rows with SPE redshift columns on object_id."""
        if zspe_table is None or len(zspe_table) == 0:
            return Table()
        if "object_id" not in input_table.colnames:
            raise ValueError("input_table must contain 'object_id'")
        if "object_id" not in zspe_table.colnames:
            return Table()

        spe_cols = ["object_id", "spe_rank", "spe_z", "spe_z_err"]
        missing = [col for col in spe_cols if col not in zspe_table.colnames]
        if missing:
            raise ValueError(f"SPE redshift result missing required columns: {missing}")

        left = input_table.copy()
        for col in ["spe_rank", "spe_z", "spe_z_err"]:
            if col in left.colnames:
                left.remove_column(col)
        left["object_id"] = [int(value) for value in left["object_id"]]

        right = zspe_table[spe_cols].copy()
        right["object_id"] = [int(value) for value in right["object_id"]]
        return join(left, right, keys="object_id", join_type="inner", metadata_conflicts="silent")

    def query_zspe_candidates(
        self,
        crossmatch_table: Optional[Table] = None,
        output_file: Optional[Union[str, Path]] = None,
        object_type: str = "qso",
        idr_field: Optional[str] = None,
        full_async: bool = False,
        async_chunk_size: int = 500000,
    ) -> Union[Table, Dict[str, Any]]:
        """
        Query IDR SPE redshift candidate catalogues for object IDs.

        WIDE queries use wide_survey first, then submit only unmatched object IDs
        to the wide table.
        """
        if not self._logged_in:
            logger.warning("Not logged in - attempting login with default credentials")
            self.login()

        if crossmatch_table is None:
            raise ValueError("Must provide crossmatch_table with 'object_id'")
        if "object_id" not in crossmatch_table.colnames:
            raise ValueError("crossmatch_table must contain 'object_id'")

        object_ids = []
        seen = set()
        for value in crossmatch_table["object_id"]:
            key = self._value_key(value)
            if key not in seen:
                seen.add(key)
                object_ids.append(int(value))

        logger.info(
            "Querying SPE redshift candidates for %d unique objects",
            len(object_ids),
        )

        user_table = Table()
        user_table["object_id"] = object_ids
        table_names = self._get_zspe_candidate_table_names(
            object_type=object_type,
            idr_field=idr_field,
        )

        normalized_object_type = object_type.lower()
        normalized_idr_field = (idr_field or "WIDE").upper()

        if full_async:
            if len(user_table) == 0:
                raise ValueError("Object ID table is empty; nothing to submit in full_async mode.")
            if output_file is None:
                raise ValueError("output_file must be provided when full_async=True.")
            if async_chunk_size <= 0:
                raise ValueError("async_chunk_size must be a positive integer")

            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            if len(user_table) > async_chunk_size:
                chunk_suffix = output_path.suffix or ".fits"
                manifest_path = output_path.with_suffix(output_path.suffix + ".manifest.json")
                batches = [
                    (i, min(i + async_chunk_size, len(user_table)))
                    for i in range(0, len(user_table), async_chunk_size)
                ]
                manifest: Dict[str, Any] = {
                    "type": "euclidkit_zspe_async_chunked",
                    "environment": self.environment,
                    "submitted_at_utc": _utc_timestamp(),
                    "row_count": len(user_table),
                    "chunk_size": async_chunk_size,
                    "chunk_count": len(batches),
                    "query_count": 0,
                    "object_type": normalized_object_type,
                    "idr_field": normalized_idr_field,
                    "tables": [table for table, _ in table_names],
                    "output_file": str(output_path),
                    "chunks": [],
                }
                part_files: List[Path] = []

                for idx, (start, end) in enumerate(batches, 1):
                    batch = user_table[start:end]
                    submission = self._submit_zspe_async_batch(batch, table_names)
                    job = submission.get("job")
                    job_id = getattr(job, "jobid", None)
                    part_path = output_path.with_name(
                        f"{output_path.stem}_part_{idx:04d}{chunk_suffix}"
                    )
                    rows = 0
                    status = "FAILED"
                    error_msg = None
                    try:
                        chunk_result = job.get_results()
                        rows = len(chunk_result) if chunk_result is not None else 0
                        save_table(chunk_result if chunk_result is not None else Table(), part_path)
                        self._delete_async_job(job_id)
                        part_files.append(part_path)
                        status = "COMPLETED"
                    except Exception as exc:
                        error_msg = str(exc)
                        manifest["chunks"].append({
                            "index": idx,
                            "row_start": start,
                            "row_end": end,
                            "job_id": job_id,
                            "status": status,
                            "rows": rows,
                            "file": str(part_path),
                            "query": submission.get("query"),
                            "error": error_msg,
                        })
                        with open(manifest_path, "w", encoding="utf-8") as fh:
                            json.dump(manifest, fh, indent=2)
                        raise

                    manifest["chunks"].append({
                        "index": idx,
                        "row_start": start,
                        "row_end": end,
                        "job_id": job_id,
                        "status": status,
                        "rows": rows,
                        "file": str(part_path),
                        "query": submission.get("query"),
                        "error": error_msg,
                    })
                    with open(manifest_path, "w", encoding="utf-8") as fh:
                        json.dump(manifest, fh, indent=2)

                merge_tables = [load_table(path) for path in part_files] if part_files else []
                zspe_result = vstack(merge_tables, metadata_conflicts="silent") if merge_tables else Table()
                final_result = self._join_zspe_with_input(crossmatch_table, zspe_result)
                manifest["query_count"] = len(manifest["chunks"])
                with open(manifest_path, "w", encoding="utf-8") as fh:
                    json.dump(manifest, fh, indent=2)
                save_table(final_result, output_path)

                return {
                    "type": "euclidkit_zspe_async_chunked",
                    "environment": self.environment,
                    "submitted_at_utc": manifest["submitted_at_utc"],
                    "row_count": len(user_table),
                    "chunk_size": async_chunk_size,
                    "chunk_count": len(batches),
                    "query_count": len(manifest["chunks"]),
                    "object_type": normalized_object_type,
                    "idr_field": normalized_idr_field,
                    "tables": [table for table, _ in table_names],
                    "job_ids": [c.get("job_id") for c in manifest["chunks"] if c.get("job_id")],
                    "results_downloaded": True,
                    "result_row_count": len(final_result),
                    "output_file": str(output_path),
                    "manifest_file": str(manifest_path),
                }

            submission = self._submit_zspe_async_batch(user_table, table_names)
            job = submission.get("job")
            job_id = getattr(job, "jobid", None)

            def _job_attr(attr: str) -> Any:
                value = getattr(job, attr, None)
                if isinstance(value, (str, int, float, bool)) or value is None:
                    return value
                return None

            job_info: Dict[str, Any] = {
                "type": "euclidkit_zspe_async_job",
                "environment": self.environment,
                "submitted_at_utc": _utc_timestamp(),
                "job_id": job_id,
                "job_phase": _job_attr("phase"),
                "job_url": _job_attr("url"),
                "results_url": _job_attr("remote_results_location"),
                "upload_table": submission.get("upload_table_name"),
                "row_count": len(user_table),
                "object_type": normalized_object_type,
                "idr_field": normalized_idr_field,
                "tables": [table for table, _ in table_names],
                "query": submission.get("query"),
                "results_downloaded": False,
                "result_row_count": 0,
                "download_error": None,
            }

            try:
                async_result = job.get_results()
                zspe_result = async_result if async_result is not None else Table()
                final_result = self._join_zspe_with_input(crossmatch_table, zspe_result)
                save_table(final_result, output_path)
                self._delete_async_job(job_id)
                job_info["results_downloaded"] = True
                job_info["result_row_count"] = len(final_result)
                job_info["output_file"] = str(output_path)
            except Exception as exc:
                job_info["download_error"] = str(exc)
                job_info_path = output_path.with_name(output_path.name + ".job.json")
                self._write_job_info(job_info, job_info_path)
                job_info["job_info_file"] = str(job_info_path)

            return job_info

        batch_size = 1000
        all_results: List[Table] = []

        for i in range(0, len(user_table), batch_size):
            batch = user_table[i:i + batch_size]
            logger.info(
                "Querying SPE redshift batch %d/%d",
                i // batch_size + 1,
                (len(user_table) - 1) // batch_size + 1,
            )

            first_table, first_label = table_names[0]
            first_result = self._query_zspe_batch(batch, first_table, first_label)
            if len(first_result) > 0:
                all_results.append(first_result)

            if len(table_names) == 2:
                matched = (
                    {self._value_key(value) for value in first_result["object_id"]}
                    if len(first_result) > 0 and "object_id" in first_result.colnames
                    else set()
                )
                mask = [self._value_key(value) not in matched for value in batch["object_id"]]
                unmatched = batch[mask]
                logger.info(
                    "IDR WIDE SPE fallback batch: %d/%d objects unmatched in wide_survey",
                    len(unmatched),
                    len(batch),
                )
                if len(unmatched) > 0:
                    fallback_table, fallback_label = table_names[1]
                    fallback_result = self._query_zspe_batch(
                        unmatched,
                        fallback_table,
                        fallback_label,
                    )
                    if len(fallback_result) > 0:
                        all_results.append(fallback_result)

        if all_results:
            zspe_result = vstack(all_results, metadata_conflicts="silent")
            final_result = self._join_zspe_with_input(crossmatch_table, zspe_result)
        else:
            logger.warning("No SPE redshift candidates found")
            final_result = Table()

        if output_file:
            save_table(final_result, output_file)
            logger.info(f"Saved SPE redshift candidate results to {output_file}")

        return final_result
    
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
