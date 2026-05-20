"""
Image cutout generation module for euclidkit package.

Provides functionality to generate cutouts from Euclid data products.
"""

import os
import re
import logging
import tempfile
from pathlib import Path
from typing import Union, Optional, List, Dict, Any, Tuple

import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy import wcs
from astropy.nddata.utils import Cutout2D
from skimage.transform import resize
from tqdm import tqdm

from euclidkit.core.data_access import EuclidArchive
from euclidkit.utils.io import load_table, save_table

logger = logging.getLogger(__name__)


class CutoutGenerator:
    """Generator for Euclid image cutouts."""
    
    def __init__(self, archive: Optional[EuclidArchive] = None, environment: str = 'PDR'):
        """
        Initialize cutout generator.
        
        Parameters
        ----------
        archive : EuclidArchive, optional
            Pre-initialized archive connection. If None, creates new one.
        environment : str
            Archive environment: 'PDR', 'OTF', 'REG', or 'IDR'
        """
        if archive is None:
            self.archive = EuclidArchive(environment=environment)
        else:
            self.archive = archive
        self._last_cutana_mosaic_query: Optional[str] = None
        
        logger.info(f"Initialized CutoutGenerator for {environment} environment")
    
    def make_cutouts(
        self,
        sources: Union[str, pd.DataFrame, Table, List[Dict], Dict],
        output_location: str,
        cutout_size: Union[float, int, u.Quantity] = 10 * u.arcsec,
        instrument_name: str = 'VIS',
        product_type: str = 'mosaic',
        nisp_filters: Optional[List[str]] = None,
        **kwargs
    ) -> int:
        """
        Generate cutouts for sources.
        
        Parameters
        ----------
        sources : str, DataFrame, Table, List[Dict], or Dict
            Source information. Can be:
            - object_id (str): Single Euclid object ID
            - Dict with 'object_id' key
            - Dict with 'ra' and 'dec' keys (in degrees)
            - DataFrame/Table with either object_id or ra/dec columns
            - List of above
        output_location : str
            Directory path where cutouts will be saved
        cutout_size : float, int, or u.Quantity, default 10*arcsec
            Size of cutout. If no units, assumed to be arcseconds.
        instrument_name : str, default 'VIS'
            Instrument name: 'VIS' or 'NISP'
        product_type : str, default 'mosaic'
            Product type: 'mosaic', 'calibrated_frame', or 'stacked_frame'
        nisp_filters : List[str], optional
            NISP filter names: ['NIR_H', 'NIR_J', 'NIR_Y']
        
        Returns
        -------
        int
            Number of cutouts successfully generated
        
        Examples
        --------
        >>> cutgen = CutoutGenerator()
        >>> # Single object by ID
        >>> cutgen.make_cutouts('123456789', '/path/to/output/')
        >>> 
        >>> # Single object by coordinates
        >>> cutgen.make_cutouts({'ra': 150.0, 'dec': 2.5}, '/path/to/output/')
        >>> 
        >>> # Multiple sources
        >>> sources = [
        ...     {'object_id': '123456789'},
        ...     {'ra': 150.0, 'dec': 2.5}
        ... ]
        >>> cutgen.make_cutouts(sources, '/path/to/output/')
        """
        # Ensure output directory exists
        os.makedirs(output_location, exist_ok=True)
        
        # Normalize input to DataFrame format
        df = self._normalize_sources(sources)
        
        if len(df) == 0:
            logger.warning("No valid sources provided")
            return 0
        
        # Ensure cutout_size has units
        if not isinstance(cutout_size, u.Quantity):
            cutout_size = cutout_size * u.arcsec
            
        logger.info(f"Generating cutouts for {len(df)} sources with {cutout_size} size")
        
        # Get file matches for sources
        df_with_files = self._get_file_matches(
            df, instrument_name, product_type, nisp_filters
        )
        
        if len(df_with_files) == 0:
            logger.warning("No file matches found for sources")
            return 0
        
        # Generate cutouts
        return self._generate_cutouts(
            df_with_files, output_location, cutout_size, product_type, **kwargs
        )

    def generate_cutana_input(
        self,
        sources: Union[str, pd.DataFrame, Table, List[Dict], Dict],
        output_file: Union[str, Path],
        instrument_name: str = 'VIS',
        nisp_filters: Optional[List[str]] = None,
        cutout_size: str = 'arcsec',
        cutout_size_value: float = 15.0,
        drop_noncutana_cols: bool = True,
        idr_field: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Generate a Cutana-compatible source catalogue from source coordinates.

        Notes
        -----
        This workflow is adapted from the Bulk_mosaic_cutouts notebook by
        Kristin Anett Remmelgas.

        Parameters
        ----------
        sources : str, DataFrame, Table, List[Dict], or Dict
            Source information with either ``object_id`` or RA/Dec columns.
        output_file : str or Path
            Output CSV path for the Cutana input file.
        instrument_name : str, default 'VIS'
            Instrument selection: 'VIS' or 'NISP'.
        nisp_filters : List[str], optional
            Optional subset of NISP filters when instrument is 'NISP'.
        cutout_size : str, default 'arcsec'
            Cutana diameter unit, either 'pixel' or 'arcsec'.
        cutout_size_value : float, default 15.0
            Constant diameter value applied to each source.
        drop_noncutana_cols : bool, default True
            When False, keep additional user columns in output.
        idr_field : {'WIDE', 'DEEP'}, optional
            IDR field selector used when resolving ``object_id`` values to
            coordinates in ``environment='IDR'``. Defaults to 'WIDE' when
            omitted, matching ``EuclidArchive.crossmatch_sources`` behavior.

        Returns
        -------
        pd.DataFrame
            DataFrame written to ``output_file``.
        """
        unit = cutout_size.lower()
        if unit not in {'pixel', 'arcsec'}:
            raise ValueError("cutout_size must be either 'pixel' or 'arcsec'")
        if instrument_name not in {'VIS', 'NISP'}:
            raise ValueError("instrument_name must be 'VIS' or 'NISP'")
        if instrument_name != 'NISP' and nisp_filters:
            raise ValueError("nisp_filters can only be used with instrument_name='NISP'")

        df_input = self._normalize_sources(sources)
        if len(df_input) == 0:
            raise ValueError("No sources found in input")

        df_coords = self._resolve_coordinates(df_input, idr_field=idr_field)
        if len(df_coords) == 0:
            raise ValueError("Could not resolve source coordinates from input")

        radec_colnames = self._get_radec_colnames(df_coords.columns.tolist())
        ra_col = radec_colnames['ra_colname']
        dec_col = radec_colnames['dec_colname']

        if 'segmentation_map_id' not in df_coords.columns:
            raise ValueError("segmentation_map_id is required to query mosaic_product")

        df_files = self._get_files_mosaic_from_segmentation(
            df_coords,
            instrument_name=instrument_name,
            nisp_filters=nisp_filters,
            radec_colnames=radec_colnames,
        )
        if len(df_files) == 0:
            query_msg = self._last_cutana_mosaic_query or "<unavailable>"
            raise ValueError(
                "No mosaic files matched the provided sources.\n"
                f"Second query was:\n{query_msg}"
            )

        df_files = df_files.copy()
        df_files['fits_file_paths'] = (
            df_files['datalabs_path'].astype(str).str.rstrip('/')
            + '/'
            + df_files['file_name'].astype(str)
        )

        # Build a stable source identifier for Cutana rows.
        if 'object_id' in df_files.columns:
            source_ids = df_files['object_id'].astype(str)
            df_files['SourceID'] = 'OBJ_' + source_ids.str.replace(r'\.0$', '', regex=True)
        else:
            df_files['SourceID'] = (
                'SRC_'
                + df_files[ra_col].astype(str).str.replace('.', '', regex=False)
                + '_'
                + df_files[dec_col].astype(str).str.replace('.', '', regex=False)
            )

        group_cols = ['SourceID', ra_col, dec_col]
        if 'object_id' in df_files.columns:
            group_cols.append('object_id')

        df_cutana = (
            df_files[group_cols + ['fits_file_paths']]
            .groupby(group_cols, as_index=False)['fits_file_paths']
            .agg(lambda values: str(sorted(set(values))))
        )

        df_cutana.insert(
            loc=3,
            column=f'diameter_{unit}',
            value=float(cutout_size_value),
        )

        if not drop_noncutana_cols:
            keep_cols = [c for c in df_coords.columns if c not in df_cutana.columns]
            if keep_cols:
                df_cutana = df_cutana.merge(
                    df_coords[[ra_col, dec_col] + keep_cols].drop_duplicates(),
                    on=[ra_col, dec_col],
                    how='left',
                )

        df_cutana = df_cutana.rename(columns={ra_col: 'RA', dec_col: 'Dec'})

        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df_cutana.to_csv(output_path, index=False)

        logger.info(
            "Created Cutana input with %d sources at %s",
            len(df_cutana),
            output_path,
        )
        return df_cutana
    
    def _normalize_sources(self, sources) -> pd.DataFrame:
        """Convert various source input formats to standardized DataFrame."""
        
        # Handle string input (object_id)
        if isinstance(sources, str):
            return pd.DataFrame({'object_id': [sources]})
        
        # Handle single dict input
        if isinstance(sources, dict):
            return pd.DataFrame([sources])
        
        # Handle list input
        if isinstance(sources, list):
            # Check if list of dicts
            if all(isinstance(item, dict) for item in sources):
                return pd.DataFrame(sources)
            else:
                # Assume list of object_ids
                return pd.DataFrame({'object_id': sources})
        
        # Handle DataFrame/Table input
        if isinstance(sources, (pd.DataFrame, Table)):
            if isinstance(sources, Table):
                sources = sources.to_pandas()
            return sources
        
        raise ValueError(f"Unsupported sources format: {type(sources)}")
    
    def _get_file_matches(
        self, 
        df: pd.DataFrame, 
        instrument_name: str,
        product_type: str,
        nisp_filters: Optional[List[str]]
    ) -> pd.DataFrame:
        """Find files matching the sources."""
        
        # Ensure we have coordinates
        df_coords = self._resolve_coordinates(df)
        
        if len(df_coords) == 0:
            return pd.DataFrame()
        
        # Check coordinate column names
        radec_colnames = self._get_radec_colnames(df_coords.columns.tolist())
        
        # Get file matches using appropriate method
        if product_type == 'mosaic':
            df_files = self._get_files_mosaic(
                df_coords, instrument_name, nisp_filters, radec_colnames
            )
        elif product_type == 'calibrated_frame':
            df_files = self._get_files_calibrated(
                df_coords, instrument_name, nisp_filters, radec_colnames
            )
        elif product_type == 'stacked_frame':
            df_files = self._get_files_stacked(
                df_coords, instrument_name, nisp_filters, radec_colnames
            )
        else:
            raise ValueError(f"Invalid product_type: {product_type}")
        
        return df_files
    
    def _resolve_coordinates(self, df: pd.DataFrame, idr_field: Optional[str] = None) -> pd.DataFrame:
        """Resolve source rows to include coordinates and segmentation_map_id."""
        
        # If object IDs are present, always resolve via MER to obtain segmentation_map_id.
        if 'object_id' in df.columns:
            return self._lookup_mer_source_info(df, idr_field=idr_field)

        # If we only have coordinates, spatially crossmatch once against MER to obtain
        # object_id and segmentation_map_id, then continue with segmentation-based joins.
        if self._has_coordinates(df):
            return self._lookup_mer_source_info_by_position(df, idr_field=idr_field)
        
        logger.error("No coordinates or object_ids found in sources")
        return pd.DataFrame()
    
    def _has_coordinates(self, df: pd.DataFrame) -> bool:
        """Check if DataFrame has coordinate columns."""
        cols = df.columns.tolist()
        return (('ra' in cols and 'dec' in cols) or 
                ('right_ascension' in cols and 'declination' in cols))
    
    def _lookup_mer_source_info(self, df: pd.DataFrame, idr_field: Optional[str] = None) -> pd.DataFrame:
        """Look up MER source metadata from object_ids via archive query."""
        
        if not self.archive._logged_in:
            self.archive.login()
        
        df = df.copy()

        # Reuse environment-specific MER selection from EuclidArchive.
        table_names = self.archive._get_mer_table_names(idr_field=idr_field)

        # Upload object IDs and perform a server-side join, matching data_access.py pattern.
        upload_df = df[['object_id']].drop_duplicates().copy()
        try:
            def _query_ids(ids_df: pd.DataFrame, table_name: str) -> Table:
                upload_table = Table.from_pandas(ids_df)
                with tempfile.NamedTemporaryFile(suffix='.vot', delete=False) as tmp_file:
                    upload_table.write(tmp_file.name, format='votable', overwrite=True)
                    tmp_name = tmp_file.name

                upload_name = f"user_cutana_oid_{np.random.randint(10000, 99999)}"
                try:
                    query = f"""
                    SELECT u.object_id, m.right_ascension, m.declination, m.segmentation_map_id
                    FROM TAP_UPLOAD.{upload_name} AS u
                    JOIN {table_name} AS m ON u.object_id = m.object_id
                    ORDER BY u.object_id
                    """
                    job = self.archive.euclid.launch_job_async(
                        query,
                        upload_resource=tmp_name,
                        upload_table_name=upload_name,
                    )
                    if job is None:
                        logger.error("Failed to resolve object_ids to MER metadata from %s", table_name)
                        return Table()
                    result = job.get_results()
                    return result if result is not None else Table()
                finally:
                    if os.path.exists(tmp_name):
                        os.unlink(tmp_name)

            result_tables = []
            if (
                len(table_names) == 2
                and table_names[0] == 'catalogue.mer_catalogue_wide_survey'
                and table_names[1] == 'catalogue.mer_catalogue_wide_mode'
            ):
                survey_result = _query_ids(upload_df, table_names[0])
                if survey_result is not None and len(survey_result) > 0:
                    result_tables.append(survey_result)
                matched_ids = set()
                if survey_result is not None and 'object_id' in survey_result.colnames:
                    matched_ids = {str(value) for value in survey_result['object_id']}
                unmatched_df = upload_df[
                    ~upload_df['object_id'].astype(str).isin(matched_ids)
                ].copy()
                logger.info(
                    "IDR WIDE Cutana fallback: %d/%d object_ids unmatched in wide_survey",
                    len(unmatched_df),
                    len(upload_df),
                )
                if len(unmatched_df) > 0:
                    mode_result = _query_ids(unmatched_df, table_names[1])
                    if mode_result is not None and len(mode_result) > 0:
                        result_tables.append(mode_result)
            else:
                for table_name in table_names:
                    result = _query_ids(upload_df, table_name)
                    if result is not None and len(result) > 0:
                        result_tables.append(result)

            coords_table = self.archive._combine_mer_results(result_tables)
            if len(coords_table) == 0:
                logger.error("Failed to resolve object_ids to MER metadata")
                return pd.DataFrame()

            coords_df = coords_table.to_pandas()
            # Ensure matching dtype to avoid pandas merge type issues.
            if 'object_id' in coords_df.columns and 'object_id' in df.columns:
                coords_df['object_id'] = coords_df['object_id'].astype(str)
                df['object_id'] = df['object_id'].astype(str)
            logger.info(
                "Resolved MER metadata for %d of %d object_ids",
                len(coords_df),
                len(upload_df),
            )
            
            # Merge back with original data. Drop potentially stale MER-style
            # columns first to avoid pandas suffixing to *_x/*_y.
            drop_if_present = [
                col for col in ['right_ascension', 'declination', 'segmentation_map_id']
                if col in df.columns
            ]
            if drop_if_present:
                df = df.drop(columns=drop_if_present)
            return df.merge(coords_df, on='object_id', how='inner')
            
        except Exception as e:
            logger.error(f"Error resolving MER metadata: {e}")
            return pd.DataFrame()

    def _lookup_mer_source_info_by_position(
        self,
        df: pd.DataFrame,
        idr_field: Optional[str] = None,
        radius_arcsec: float = 1.0,
    ) -> pd.DataFrame:
        """
        Resolve MER source metadata from input coordinates via archive crossmatch.

        Uses EuclidArchive.crossmatch_sources to reuse environment-specific MER table
        selection and ensure segmentation_map_id is retrieved.
        """
        if not self.archive._logged_in:
            self.archive.login()

        radec_colnames = self._get_radec_colnames(df.columns.tolist())
        ra_col = radec_colnames['ra_colname']
        dec_col = radec_colnames['dec_colname']

        try:
            match_table = self.archive.crossmatch_sources(
                user_table=Table.from_pandas(df.copy()),
                radius=radius_arcsec,
                ra_col=ra_col,
                dec_col=dec_col,
                use_object_id=False,
                idr_field=idr_field,
            )
            if len(match_table) == 0:
                return pd.DataFrame()

            match_df = match_table.to_pandas()
            # Keep nearest match per input coordinate.
            if 'separation_arcsec' in match_df.columns:
                match_df = match_df.sort_values('separation_arcsec').drop_duplicates(
                    subset=[ra_col, dec_col],
                    keep='first',
                )

            # Use MER coordinates when available for consistency downstream.
            if 'mer_ra' in match_df.columns:
                match_df['right_ascension'] = match_df['mer_ra']
            elif ra_col in match_df.columns:
                match_df['right_ascension'] = match_df[ra_col]

            if 'mer_dec' in match_df.columns:
                match_df['declination'] = match_df['mer_dec']
            elif dec_col in match_df.columns:
                match_df['declination'] = match_df[dec_col]

            if 'segmentation_map_id' not in match_df.columns:
                logger.error("Crossmatch did not return segmentation_map_id")
                return pd.DataFrame()

            # Preserve user columns by merging back on input coordinate columns.
            keep_cols = [c for c in ['object_id', 'right_ascension', 'declination', 'segmentation_map_id'] if c in match_df.columns]
            select_cols = []
            for col in [ra_col, dec_col] + keep_cols:
                if col not in select_cols:
                    select_cols.append(col)
            merged = df.merge(
                match_df[select_cols].drop_duplicates(),
                on=[ra_col, dec_col],
                how='inner',
            )
            return merged
        except Exception as e:
            logger.error(f"Error resolving MER metadata by position: {e}")
            return pd.DataFrame()

    def _get_files_mosaic_from_segmentation(
        self,
        df: pd.DataFrame,
        instrument_name: str,
        nisp_filters: Optional[List[str]],
        radec_colnames: Dict[str, str],
    ) -> pd.DataFrame:
        """
        Get mosaic files using segmentation_map_id -> tile_index shortcut.

        This avoids per-source positional CONTAINS crossmatch by joining source rows
        to ``sedm.mosaic_product`` on the first 9 digits of segmentation_map_id.
        """
        if not self.archive._logged_in:
            self.archive.login()
        self.archive._ensure_client()

        required_cols = [radec_colnames['ra_colname'], radec_colnames['dec_colname'], 'segmentation_map_id']
        missing_cols = [c for c in required_cols if c not in df.columns]
        if missing_cols:
            self._last_cutana_mosaic_query = (
                f"<skipped: missing columns for segmentation shortcut: {missing_cols}; "
                f"available columns={list(df.columns)}>"
            )
            logger.warning("Missing columns for segmentation shortcut: %s", missing_cols)
            return pd.DataFrame()

        df_work = df[required_cols + ([c for c in ['object_id'] if c in df.columns])].copy()
        df_work = df_work[df_work['segmentation_map_id'].notna()]
        if len(df_work) == 0:
            self._last_cutana_mosaic_query = (
                "<skipped: segmentation_map_id present but all values are null/NaN after preprocessing>"
            )
            logger.warning("Segmentation shortcut skipped: all segmentation_map_id values are null/NaN.")
            return pd.DataFrame()

        upload_table = Table.from_pandas(df_work)
        with tempfile.NamedTemporaryFile(suffix='.vot', delete=False) as tmp_file:
            upload_table.write(tmp_file.name, format='votable', overwrite=True)
            tmp_name = tmp_file.name

        ra_col = radec_colnames['ra_colname']
        dec_col = radec_colnames['dec_colname']
        upload_name = f"user_cutana_seg_{np.random.randint(10000, 99999)}"
        object_id_select = "u.object_id," if 'object_id' in df_work.columns else ""

        query = f"""
        SELECT
            {object_id_select}
            u.{ra_col} AS {ra_col},
            u.{dec_col} AS {dec_col},
            u.segmentation_map_id,
            mosaics.file_name, mosaics.file_path, mosaics.datalabs_path,
            mosaics.mosaic_product_oid, mosaics.tile_index, mosaics.instrument_name,
            mosaics.filter_name, mosaics.ra AS image_ra, mosaics.dec AS image_dec
        FROM TAP_UPLOAD.{upload_name} AS u
        JOIN sedm.mosaic_product AS mosaics
          ON CAST(mosaics.tile_index AS CHAR(25)) = SUBSTRING(CAST(u.segmentation_map_id AS CHAR(25)), 1, 9)
        WHERE mosaics.datalabs_path IS NOT NULL
          AND mosaics.instrument_name='{instrument_name}'
        """

        if instrument_name == 'NISP' and nisp_filters:
            filter_conditions = " OR ".join([f"mosaics.filter_name='{f}'" for f in nisp_filters])
            query += f" AND ({filter_conditions})"
        query += f" ORDER BY u.{ra_col} ASC"
        self._last_cutana_mosaic_query = query.strip()
        logger.warning("Cutana mosaic query:\n%s", self._last_cutana_mosaic_query)

        try:
            job = self.archive.euclid.launch_job_async(
                query,
                upload_resource=tmp_name,
                upload_table_name=upload_name,
            )

            if job is None:
                logger.error("Query for segmentation-based mosaic files failed")
                return pd.DataFrame()

            df_result = job.get_results().to_pandas()
            if len(df_result) == 0:
                logger.warning("Cutana mosaic query returned 0 rows.")
                return pd.DataFrame()

            # Keep one record per source/filter when duplicates are returned.
            df_result = df_result.drop_duplicates(subset=[ra_col, 'filter_name'], keep='first')
            logger.info("Found %d mosaic file matches via segmentation shortcut", len(df_result))
            return df_result
        except Exception as e:
            logger.error(f"Error querying segmentation-based mosaic files: {e}")
            return pd.DataFrame()
        finally:
            if os.path.exists(tmp_name):
                os.unlink(tmp_name)
    
    def _get_radec_colnames(self, columns: List[str]) -> Dict[str, str]:
        """Determine RA/Dec column names."""
        if 'ra' in columns and 'dec' in columns:
            return {"ra_colname": 'ra', "dec_colname": 'dec'}
        else:
            return {"ra_colname": 'right_ascension', "dec_colname": 'declination'}
    
    def _build_source_query(self, df: pd.DataFrame, radec_colnames: Dict[str, str]) -> str:
        """Legacy helper retained for compatibility; no longer used.

        Some TAP services used by Euclid do not accept ADQL with VALUES,
        so we now query per-source instead of joining a values table.
        """
        return ""
    
    def _get_files_mosaic(
        self, 
        df: pd.DataFrame,
        instrument_name: str, 
        nisp_filters: Optional[List[str]],
        radec_colnames: Dict[str, str],
        searchdist: float = 0.5
    ) -> pd.DataFrame:
        """Get mosaic files for sources."""
        
        if not self.archive._logged_in:
            self.archive.login()
        
        all_results: List[pd.DataFrame] = []
        radius_deg = float(searchdist) / 60.0
        
        for _, row in df.iterrows():
            ra = row[radec_colnames['ra_colname']]
            dec = row[radec_colnames['dec_colname']]
            query = f"""
            SELECT {ra} AS {radec_colnames['ra_colname']},
                   {dec} AS {radec_colnames['dec_colname']},
                   mosaics.file_name, mosaics.file_path, mosaics.datalabs_path,
                   mosaics.mosaic_product_oid, mosaics.tile_index, mosaics.instrument_name,
                   mosaics.filter_name, mosaics.ra AS image_ra, mosaics.dec AS image_dec
            FROM sedm.mosaic_product AS mosaics
            WHERE mosaics.datalabs_path IS NOT NULL
              AND mosaics.instrument_name='{instrument_name}'
              AND (mosaics.fov IS NOT NULL AND CONTAINS(CIRCLE('ICRS', {ra}, {dec}, {radius_deg}), mosaics.fov)=1)
            """
            if instrument_name == 'NISP' and nisp_filters:
                filter_conditions = " OR ".join([f"mosaics.filter_name='{f}'" for f in nisp_filters])
                query += f" AND ({filter_conditions})"
            try:
                job = self.archive.euclid.launch_job_async(query, verbose=False)
                if job is None:
                    logger.error("Query for mosaic files failed for one source")
                    continue
                df_result = job.get_results().to_pandas()
                if len(df_result) > 0:
                    all_results.append(df_result)
            except Exception as e:
                logger.error(f"Error querying mosaic files for RA={ra}, DEC={dec}: {e}")
                continue
        
        if not all_results:
            return pd.DataFrame()
        
        df_result = pd.concat(all_results, ignore_index=True)
        logger.info(f"Found {len(df_result)} mosaic file matches")
        
        if len(df_result) > 0:
            df_result = df_result.drop_duplicates(
                subset=[radec_colnames["ra_colname"], 'filter_name'],
                keep='first',
            )
        
        return df_result
    
    def _get_files_calibrated(
        self, 
        df: pd.DataFrame,
        instrument_name: str, 
        nisp_filters: Optional[List[str]],
        radec_colnames: Dict[str, str],
        searchdist: float = 0.5
    ) -> pd.DataFrame:
        """Get calibrated frame files for sources."""
        
        if not self.archive._logged_in:
            self.archive.login()
        
        all_results: List[pd.DataFrame] = []
        radius_deg = float(searchdist) / 60.0
        
        for _, row in df.iterrows():
            ra = row[radec_colnames['ra_colname']]
            dec = row[radec_colnames['dec_colname']]
            query = f"""
            SELECT {ra} AS {radec_colnames['ra_colname']},
                   {dec} AS {radec_colnames['dec_colname']},
                   frames.file_name, frames.file_path, frames.datalabs_path,
                   frames.instrument_name, frames.filter_name, frames.observation_id,
                   frames.ra AS image_ra, frames.dec AS image_dec
            FROM sedm.calibrated_frame AS frames
            WHERE frames.datalabs_path IS NOT NULL
              AND frames.instrument_name='{instrument_name}'
              AND (frames.fov IS NOT NULL AND CONTAINS(CIRCLE('ICRS', {ra}, {dec}, {radius_deg}), frames.fov)=1)
            """
            if instrument_name == 'NISP' and nisp_filters:
                filter_conditions = " OR ".join([f"frames.filter_name='{f}'" for f in nisp_filters])
                query += f" AND ({filter_conditions})"
            try:
                job = self.archive.euclid.launch_job_async(query, verbose=False)
                if job is None:
                    logger.error("Query for calibrated frame files failed for one source")
                    continue
                df_result = job.get_results().to_pandas()
                if len(df_result) > 0:
                    all_results.append(df_result)
            except Exception as e:
                logger.error(f"Error querying calibrated frame files for RA={ra}, DEC={dec}: {e}")
                continue
        
        if not all_results:
            return pd.DataFrame()
        
        df_result = pd.concat(all_results, ignore_index=True)
        logger.info(f"Found {len(df_result)} calibrated frame file matches")
        return df_result
    
    def _get_files_stacked(
        self, 
        df: pd.DataFrame,
        instrument_name: str, 
        nisp_filters: Optional[List[str]],
        radec_colnames: Dict[str, str],
        searchdist: float = 0.5
    ) -> pd.DataFrame:
        """Get stacked frame files for sources."""
        
        if not self.archive._logged_in:
            self.archive.login()
        
        all_results: List[pd.DataFrame] = []
        radius_deg = float(searchdist) / 60.0
        
        for _, row in df.iterrows():
            ra = row[radec_colnames['ra_colname']]
            dec = row[radec_colnames['dec_colname']]
            query = f"""
            SELECT {ra} AS {radec_colnames['ra_colname']},
                   {dec} AS {radec_colnames['dec_colname']},
                   frames.file_name, frames.file_path, frames.datalabs_path,
                   frames.instrument_name, frames.filter_name, frames.observation_id,
                   frames.ra AS image_ra, frames.dec AS image_dec
            FROM observation_stack AS frames
            WHERE frames.datalabs_path IS NOT NULL
              AND frames.instrument_name='{instrument_name}'
              AND (frames.fov IS NOT NULL AND CONTAINS(CIRCLE('ICRS', {ra}, {dec}, {radius_deg}), frames.fov)=1)
            """
            try:
                job = self.archive.euclid.launch_job_async(query, verbose=False)
                if job is None:
                    logger.error("Query for stacked frame files failed for one source")
                    continue
                df_result = job.get_results().to_pandas()
                if len(df_result) > 0:
                    all_results.append(df_result)
            except Exception as e:
                logger.error(f"Error querying stacked frame files for RA={ra}, DEC={dec}: {e}")
                continue
        
        if not all_results:
            return pd.DataFrame()
        
        df_result = pd.concat(all_results, ignore_index=True)
        logger.info(f"Found {len(df_result)} stacked frame file matches")
        
        if len(df_result) > 0:
            df_result = df_result.drop_duplicates(
                subset=[radec_colnames["ra_colname"], 'filter_name'],
                keep='first',
            )
        
        return df_result
    
    def _generate_cutouts(
        self,
        df: pd.DataFrame,
        output_location: str,
        cutout_size: u.Quantity,
        product_type: str,
        **kwargs
    ) -> int:
        """Generate the actual cutout files."""
        
        radec_colnames = self._get_radec_colnames(df.columns.tolist())
        total_cutouts = 0
        
        # Group by file to minimize I/O
        for file_name, sub_df in tqdm(df.groupby('file_name'), desc="Processing files"):
            cutouts_made = self._make_cutouts_from_file(
                sub_df, output_location, cutout_size, product_type, radec_colnames, **kwargs
            )
            total_cutouts += cutouts_made
        
        logger.info(f"Generated {total_cutouts} cutouts total")
        return total_cutouts
    
    def _make_cutouts_from_file(
        self,
        sub_df: pd.DataFrame,
        output_location: str,
        cutout_size: u.Quantity,
        product_type: str,
        radec_colnames: Dict[str, str],
        **kwargs
    ) -> int:
        """Make cutouts from a single file."""
        
        file_path = sub_df['datalabs_path'].iloc[0] + '/' + sub_df['file_name'].iloc[0]
        file_name = sub_df['file_name'].iloc[0]
        
        try:
            with fits.open(file_path) as hdul:
                cutouts_made = 0
                
                for _, row in sub_df.iterrows():
                    ra = row[radec_colnames['ra_colname']]
                    dec = row[radec_colnames['dec_colname']]
                    
                    success = self._make_single_cutout(
                        ra, dec, file_name, hdul, product_type, cutout_size, output_location, **kwargs
                    )
                    if success:
                        cutouts_made += 1
                
                return cutouts_made
                
        except Exception as e:
            logger.error(f"Error processing file {file_name}: {e}")
            return 0
    
    def _make_single_cutout(
        self,
        ra: float,
        dec: float,
        file_name: str,
        hdul: fits.HDUList,
        product_type: str,
        cutout_size: u.Quantity,
        output_location: str,
        **kwargs
    ) -> bool:
        """Make a single cutout from an opened FITS file."""
        
        # Determine which extension to use
        if product_type == 'mosaic':
            ext_nr = 0  # Use primary extension for mosaics
        elif product_type in ['calibrated_frame', 'frame']:
            # Find appropriate detector extension
            if "_NIR_" in file_name:
                ext_nr = self._get_nisp_extension(hdul, ra, dec)
            elif "_VIS_" in file_name:
                ext_nr = self._get_vis_extension(hdul, ra, dec)
            else:
                ext_nr = -1
        elif product_type in ['stacked_frame', 'stack']:
            try:
                ext_nr = hdul.index_of("SCI")
            except KeyError:
                ext_nr = 1  # Fallback to extension 1
        else:
            ext_nr = -1
        
        # Check if source is in a valid region
        if ext_nr == -1:
            logger.warning(f"Source at ({ra:.6f}, {dec:.6f}) not in valid detector region")
            return False
        
        try:
            # Get data and header
            data = hdul[ext_nr].data
            header = hdul[ext_nr].header
            
            # Create cutout
            position = SkyCoord(ra, dec, frame='icrs', unit="deg")
            cutout = Cutout2D(
                data=data, 
                position=position, 
                size=cutout_size, 
                wcs=wcs.WCS(header),
                mode='trim'  # Default to trim mode
            )
            
            # Create output HDU
            hdu_cutout = fits.PrimaryHDU(cutout.data)
            hdu_cutout.header.update(cutout.wcs.to_header())
            
            # Generate output filename
            if product_type == 'mosaic':
                match = re.search(r"\-(.*\-.*?)_", file_name)
                if match:
                    base_name = match.group(1)
                else:
                    base_name = file_name.split('.')[0]
            else:
                match = re.search(r"\_(.*)\_", file_name)
                if match:
                    base_name = match.group(1)
                else:
                    base_name = file_name.split('.')[0]
            
            cutout_name = f"{base_name}-CUTOUT_{ra:.7f}_{dec:.7f}.fits"
            output_path = os.path.join(output_location, cutout_name)
            
            # Save cutout
            hdu_cutout.writeto(output_path, overwrite=True)
            logger.debug(f"Created cutout: {cutout_name}")
            
            return True
            
        except Exception as e:
            logger.error(f"Error creating cutout for source at ({ra:.6f}, {dec:.6f}): {e}")
            return False
    
    def _get_vis_extension(self, hdul: fits.HDUList, ra: float, dec: float) -> int:
        """Get appropriate VIS detector extension for coordinates."""
        
        # VIS detector parameters (from reference code)
        mean_gap_x, mean_gap_y = 134.7317209303794, 612.8441331180677
        det_x, det_y = 2048, 2066
        
        # Get detector index
        det_idx_x, det_idx_y = self._get_detector_index(
            hdul, ra, dec, mean_gap_x, mean_gap_y, 2*det_x, 2*det_y
        )
        
        if det_idx_x > 6 or det_idx_x < 1 or det_idx_y > 6 or det_idx_y < 1:
            return -1
        
        try:
            # Find E quadrant extension
            ext_idx_E = hdul.index_of(f"{det_idx_y}-{det_idx_x}.E.SCI")
            
            # Check all 4 quadrants to find the one containing the source
            for ext in range(ext_idx_E, ext_idx_E+12, 3):
                if ext >= len(hdul):
                    continue
                    
                w = wcs.WCS(hdul[ext].header)
                px, py = w.wcs_world2pix(ra, dec, 1)
                if 0 < py < det_y and 0 < px < det_x:
                    return ext
                    
        except (KeyError, IndexError):
            pass
        
        return -1
    
    def _get_nisp_extension(self, hdul: fits.HDUList, ra: float, dec: float) -> int:
        """Get appropriate NISP detector extension for coordinates."""
        
        # NISP detector parameters (from reference code)
        mean_gap_x, mean_gap_y = 164.90065775050047, 328.63495367441845
        det_dim = 2040
        
        # Get detector index
        det_idx_x, det_idx_y = self._get_detector_index(
            hdul, ra, dec, mean_gap_x, mean_gap_y, det_dim, det_dim
        )
        
        if det_idx_x > 4 or det_idx_x < 1 or det_idx_y > 4 or det_idx_y < 1:
            return -1
        
        try:
            # Find detector extension
            ext_idx = hdul.index_of(f"DET{det_idx_y}{det_idx_x}.SCI")
            
            # Check if source is within detector bounds
            w = wcs.WCS(hdul[ext_idx].header)
            px, py = w.wcs_world2pix(ra, dec, 1)
            if 0 < py < det_dim and 0 < px < det_dim:
                return ext_idx
                
        except (KeyError, IndexError):
            pass
        
        return -1
    
    def _get_detector_index(
        self,
        hdul: fits.HDUList,
        ra: float,
        dec: float,
        mean_gap_x: float,
        mean_gap_y: float,
        det_width: float,
        det_height: float
    ) -> Tuple[int, int]:
        """Calculate detector index for given coordinates."""
        
        # Get pixel coordinates in first detector frame
        w1 = wcs.WCS(hdul[1].header)
        px, py = w1.wcs_world2pix(ra, dec, 1)
        
        # Check if coordinates are valid
        if px < 0 or py < 0:
            return -1, -1
        
        # Adjust for gap offset
        px = px + mean_gap_x / 2
        py = py + mean_gap_y / 2
        
        # Calculate detector indices
        det_idx_x = int(np.ceil(px / (det_width + mean_gap_x)))
        det_idx_y = int(np.ceil(py / (det_height + mean_gap_y)))
        
        return det_idx_x, det_idx_y


# Convenience functions for backward compatibility and ease of use
def make_cutouts(
    sources: Union[str, pd.DataFrame, Table, List[Dict], Dict],
    output_location: str,
    cutout_size: Union[float, int, u.Quantity] = 10 * u.arcsec,
    instrument_name: str = 'VIS',
    environment: str = 'PDR',
    **kwargs
) -> int:
    """
    Convenience function to generate cutouts.
    
    Parameters
    ----------
    sources : various
        Source information (see CutoutGenerator.make_cutouts for details)
    output_location : str
        Output directory path
    cutout_size : float, int, or u.Quantity, default 10*arcsec
        Cutout size
    instrument_name : str, default 'VIS'
        Instrument name
    environment : str, default 'PDR'
        Archive environment
    **kwargs
        Additional arguments passed to CutoutGenerator.make_cutouts
    
    Returns
    -------
    int
        Number of cutouts generated
    
    Examples
    --------
    >>> from euclidkit.core.cutouts import make_cutouts
    >>> make_cutouts('123456789', '/tmp/cutouts/')
    >>> make_cutouts({'ra': 150.0, 'dec': 2.5}, '/tmp/cutouts/', cutout_size=20)
    """
    generator = CutoutGenerator(environment=environment)
    return generator.make_cutouts(
        sources, output_location, cutout_size, instrument_name, **kwargs
    )
