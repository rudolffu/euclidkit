"""
Spectral data handling module for euclidkit package.

Provides functionality for loading, processing, and combining spectral data.
"""

import os
import io
import re
import zipfile
import logging
import tempfile
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from collections import defaultdict
from typing import Union, Optional, List, Dict, Any, TYPE_CHECKING, Iterable

if TYPE_CHECKING:
    from astropy.io.fits.hdu.base import ExtensionHDU
import warnings

import numpy as np
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm

from euclidkit.utils.io import ensure_dir

logger = logging.getLogger(__name__)

try:
    import fitsio  # type: ignore
except Exception:  # pragma: no cover
    fitsio = None


class SpectrumLoader:
    """Loader for individual Euclid spectra."""

    @staticmethod
    def _as_hdu(data: np.ndarray) -> 'ExtensionHDU':
        """Convert a numpy payload to an Astropy extension HDU."""
        if getattr(data, 'dtype', None) is not None and data.dtype.names:
            return fits.BinTableHDU(data=data)
        return fits.ImageHDU(data=data)

    def load_spectra_batch(self, file_path: str, hdu_indices: Iterable[int]) -> Dict[int, 'ExtensionHDU']:
        """
        Load multiple HDUs from a single FITS file with one open call.

        Uses fitsio when available for faster random extension reads, with
        astropy fallback for compatibility.
        """
        unique_indices = sorted({int(idx) for idx in hdu_indices})
        loaded: Dict[int, 'ExtensionHDU'] = {}

        if not unique_indices:
            return loaded

        if fitsio is not None:
            try:
                with fitsio.FITS(file_path, 'r') as ff:  # type: ignore[attr-defined]
                    n_hdus = len(ff)
                    for idx in unique_indices:
                        if idx < 0 or idx >= n_hdus:
                            raise ValueError(f"HDU index {idx} out of range for {file_path}")
                        data = ff[idx].read()
                        loaded[idx] = self._as_hdu(data)
                return loaded
            except Exception as e:
                logger.warning("fitsio batch read failed for %s: %s; falling back to astropy", file_path, e)

        with fits.open(file_path, memmap=True, lazy_load_hdus=True, do_not_scale_image_data=True) as hdul:
            n_hdus = len(hdul)
            for idx in unique_indices:
                if idx < 0 or idx >= n_hdus:
                    raise ValueError(f"HDU index {idx} out of range for {file_path}")
                loaded[idx] = hdul[idx].copy()
        return loaded

    def load_spectrum_from_datalink(
        self,
        euclid_client: Any,
        source_id: Union[str, int],
        schema: str = 'sedm',
        retrieval_type: str = 'SPECTRA_RGS',
    ) -> 'ExtensionHDU':
        """
        Retrieve a spectrum via Euclid datalink without extracting per-source files.

        This bypasses ``get_spectrum`` directory checks by using ``load_data`` and
        reading the returned zip in-memory.
        """
        params_dict = {
            'ID': f"{schema} {source_id}",
            'SCHEMA': schema,
            'RETRIEVAL_TYPE': str(retrieval_type),
            'USE_ZIP_ALWAYS': 'true',
            'TAPCLIENT': 'ASTROQUERY',
        }

        def _download_with_available_clients(output_file: str) -> None:
            """
            Use astroquery's private data client path first.

            This mirrors get_spectrum internals more closely than calling the
            public EuclidClass.load_data wrapper.
            """
            private_client = getattr(euclid_client, '_EuclidClass__eucliddata', None)
            if private_client is not None:
                private_client.load_data(params_dict=params_dict, output_file=output_file, verbose=False)
                return

            logger.warning(
                "Private __eucliddata client unavailable for source_id=%s; falling back to public load_data",
                source_id
            )
            euclid_client.load_data(params_dict=params_dict, output_file=output_file, verbose=False)

        def _read_zip_spectrum(zip_path: str) -> 'ExtensionHDU':
            with zipfile.ZipFile(zip_path, 'r') as zf:
                fits_members = [n for n in zf.namelist() if n.lower().endswith('.fits')]
                if not fits_members:
                    raise ValueError(f"No FITS member found in datalink payload for source_id={source_id}")

                member_name = fits_members[0]
                payload = zf.read(member_name)
                with fits.open(io.BytesIO(payload), memmap=False, lazy_load_hdus=True, do_not_scale_image_data=True) as hdul:
                    if len(hdul) > 1:
                        return hdul[1].copy()
                    return hdul[0].copy()

        with tempfile.NamedTemporaryFile(suffix='.fits.zip', delete=False) as tmp_file:
            tmp_name = tmp_file.name

        try:
            _download_with_available_clients(tmp_name)
            # Some server-side 500s can still leave empty files; validate and retry once.
            if (not os.path.exists(tmp_name)) or os.path.getsize(tmp_name) == 0:
                private_client = getattr(euclid_client, '_EuclidClass__eucliddata', None)
                if private_client is not None:
                    logger.warning(
                        "Empty datalink payload for source_id=%s; retrying with private __eucliddata",
                        source_id
                    )
                    private_client.load_data(params_dict=params_dict, output_file=tmp_name, verbose=False)
                else:
                    euclid_client.load_data(params_dict=params_dict, output_file=tmp_name, verbose=False)
            return _read_zip_spectrum(tmp_name)
        finally:
            if os.path.exists(tmp_name):
                os.unlink(tmp_name)
    
    def load_spectrum(self, file_path: str, hdu_index: int = 1) -> 'ExtensionHDU':
        """
        Load a single spectrum from FITS file.
        
        Parameters
        ----------
        file_path : str
            Path to FITS file
        hdu_index : int
            HDU index containing the spectrum
            
        Returns
        -------
        ExtensionHDU
            Spectrum HDU
        """
        try:
            return self.load_spectra_batch(file_path, [hdu_index])[int(hdu_index)]
        except Exception as e:
            logger.error(f"Failed to load spectrum from {file_path}[{hdu_index}]: {e}")
            raise
    
    def verify_spectrum_path(self, datalabs_path: str, file_name: str) -> str:
        """
        Verify and construct full spectrum file path.
        
        Parameters
        ----------
        datalabs_path : str
            Datalabs path from archive query
        file_name : str
            FITS file name
            
        Returns
        -------
        str
            Full file path
        """
        full_path = os.path.join(datalabs_path, file_name)
        
        if not os.path.exists(full_path):
            logger.warning(f"Spectrum file not found: {full_path}")
            return None
        
        return full_path


class SpectrumCompiler:
    """Compiler for creating multi-extension FITS files from individual spectra."""
    
    def __init__(self, max_extensions: int = 1000):
        """
        Initialize spectrum compiler.
        
        Parameters
        ----------
        max_extensions : int, default 1000
            Maximum number of extensions per output file
        """
        self.max_extensions = max_extensions
        self.loader = SpectrumLoader()

    @staticmethod
    def _normalize_spectrum_filename(file_name: str) -> str:
        """Normalize FITS names so table '.fits' and XML '.fits.gz' match."""
        if file_name is None:
            return ""
        name = Path(str(file_name)).name.strip().lower()
        if name.endswith(".gz"):
            name = name[:-3]
        return name

    @staticmethod
    def _resolve_source_id_column(
        table: Table,
        preferred: str = "source_id",
    ) -> str:
        """Resolve the identifier column used for spectrum compilation."""
        if preferred in table.colnames:
            return preferred
        if preferred == "source_id" and "object_id" in table.colnames:
            logger.info(
                "Column 'source_id' not found; using 'object_id' as spectrum identifier"
            )
            return "object_id"
        raise ValueError(f"Missing required column: {preferred}")

    @staticmethod
    def deduplicate_by_source_id(
        spectra_table: Table,
        source_id_col: str = "source_id",
    ) -> tuple[Table, Dict[str, int]]:
        """
        Keep the first row for each source_id, preserving input order.

        Returns
        -------
        tuple
            (deduplicated_table, stats)
        """
        source_id_col = SpectrumCompiler._resolve_source_id_column(
            spectra_table,
            preferred=source_id_col,
        )

        seen = set()
        keep_indices: List[int] = []
        duplicate_rows = 0
        for idx, row in enumerate(spectra_table):
            sid = str(row[source_id_col])
            if sid in seen:
                duplicate_rows += 1
                continue
            seen.add(sid)
            keep_indices.append(idx)

        dedup = spectra_table[keep_indices] if keep_indices else spectra_table[:0]
        stats = {
            "input_rows": len(spectra_table),
            "unique_sources": len(dedup),
            "duplicate_rows_removed": duplicate_rows,
        }
        return dedup, stats

    @staticmethod
    def _resolve_radec_columns(
        table: Table,
        ra_candidates: Optional[List[str]] = None,
        dec_candidates: Optional[List[str]] = None,
    ) -> tuple[Optional[str], Optional[str]]:
        """Resolve RA/Dec columns from common aliases."""
        if ra_candidates is None:
            ra_candidates = [
                'ra_obj', 'right_ascension', 'RA', 'ra', 'ra_deg', 'RAJ2000',
                'right_ascension_euclid', 'ra_euclid',
            ]
        if dec_candidates is None:
            dec_candidates = [
                'dec_obj', 'declination', 'DEC', 'dec', 'dec_deg', 'DEJ2000',
                'declination_euclid', 'dec_euclid',
            ]

        ra_col = next((c for c in ra_candidates if c in table.colnames), None)
        dec_col = next((c for c in dec_candidates if c in table.colnames), None)
        return ra_col, dec_col

    def _parse_lambda_xml(self, xml_path: Path) -> Dict[str, Any]:
        """Parse one XML file and return LambdaRange + normalized FITS file names."""
        result = {"lambda_range": None, "file_names": []}
        try:
            root = ET.parse(xml_path).getroot()
        except Exception as exc:
            logger.warning("Failed to parse XML '%s': %s", xml_path, exc)
            return result

        lambda_range = None
        file_names: List[str] = []
        for elem in root.iter():
            tag = elem.tag.split("}")[-1]
            text = (elem.text or "").strip()
            if not text:
                continue
            if tag == "LambdaRange":
                lambda_range = text.upper()
            elif tag == "FileName" and ".fits" in text.lower():
                file_names.append(self._normalize_spectrum_filename(text))

        if lambda_range not in {"RGS", "BGS"}:
            lambda_range = None

        result["lambda_range"] = lambda_range
        result["file_names"] = file_names
        return result

    def annotate_lambda_range_from_xml(
        self,
        spectra_table: Table,
        datalabs_path_col: str = "datalabs_path",
        file_name_col: str = "file_name",
        lambda_col: str = "lambda_range",
    ) -> tuple[Table, Dict[str, int]]:
        """
        Annotate local spectra rows with LambdaRange by parsing XML in each datalabs path.

        Returns
        -------
        tuple
            (annotated_table, stats)
        """
        if datalabs_path_col not in spectra_table.colnames or file_name_col not in spectra_table.colnames:
            raise ValueError(
                f"Missing required columns for XML lambda resolution: {datalabs_path_col}, {file_name_col}"
            )

        table = spectra_table.copy()
        path_to_index: Dict[str, Dict[str, set[str]]] = {}
        unresolved = 0
        ambiguous = 0
        resolved = 0
        rgs_count = 0
        bgs_count = 0

        lambda_values: List[str] = []
        for row in table:
            dpath = str(row[datalabs_path_col])
            fname_norm = self._normalize_spectrum_filename(str(row[file_name_col]))

            if dpath not in path_to_index:
                mapping: Dict[str, set[str]] = defaultdict(set)
                xml_files = sorted(Path(dpath).glob("*.xml"))
                if not xml_files:
                    logger.warning("No XML files found under datalabs path: %s", dpath)
                for xml_file in xml_files:
                    parsed = self._parse_lambda_xml(xml_file)
                    lam = parsed.get("lambda_range")
                    if lam is None:
                        continue
                    for xml_name in parsed.get("file_names", []):
                        mapping[xml_name].add(lam)
                path_to_index[dpath] = mapping

            lambdas = path_to_index[dpath].get(fname_norm, set())
            if len(lambdas) == 1:
                lam = next(iter(lambdas))
                lambda_values.append(lam)
                resolved += 1
                if lam == "RGS":
                    rgs_count += 1
                elif lam == "BGS":
                    bgs_count += 1
            elif len(lambdas) > 1:
                lambda_values.append("AMBIGUOUS")
                ambiguous += 1
            else:
                lambda_values.append("UNKNOWN")
                unresolved += 1

        table[lambda_col] = np.asarray(lambda_values, dtype="U16")
        stats = {
            "total": len(table),
            "resolved": resolved,
            "rgs": rgs_count,
            "bgs": bgs_count,
            "unresolved": unresolved,
            "ambiguous": ambiguous,
        }
        return table, stats

    def filter_table_by_lambda_range(
        self,
        spectra_table: Table,
        lambda_range: str = "RGS",
        lambda_col: str = "lambda_range",
    ) -> tuple[Table, Dict[str, int]]:
        """Filter an annotated table by LambdaRange selection."""
        if lambda_col not in spectra_table.colnames:
            raise ValueError(f"Missing required lambda column: {lambda_col}")

        target = str(lambda_range).upper()
        if target not in {"RGS", "BGS", "BOTH"}:
            raise ValueError("lambda_range must be one of: RGS, BGS, BOTH")

        values = np.asarray(spectra_table[lambda_col]).astype(str)
        total = len(spectra_table)
        is_resolved = np.isin(values, ["RGS", "BGS"])
        selected_mask = is_resolved if target == "BOTH" else (values == target)

        selected = spectra_table[selected_mask]
        selected_values = np.asarray(selected[lambda_col]).astype(str) if len(selected) > 0 else np.asarray([])
        stats = {
            "total": total,
            "resolved": int(np.sum(is_resolved)),
            "selected": int(np.sum(selected_mask)),
            "selected_rgs": int(np.sum(selected_values == "RGS")),
            "selected_bgs": int(np.sum(selected_values == "BGS")),
            "unresolved": int(np.sum(values == "UNKNOWN")),
            "ambiguous": int(np.sum(values == "AMBIGUOUS")),
            "skipped_other_band": int(np.sum(is_resolved & ~selected_mask)),
        }
        return selected, stats

    def _compile_single_chunk_local(
        self,
        chunk_number: int,
        chunk: Table,
        output_file: Path,
        source_id_col: str,
        datalabs_path_col: str,
        file_name_col: str,
        hdu_index_col: str,
        overwrite: bool,
        show_progress: bool,
    ) -> Dict[str, Any]:
        """Compile one local-files chunk and return chunk stats."""
        if output_file.exists() and not overwrite:
            logger.info(f"Skipping existing chunk (resume mode): {output_file.name}")
            return {'ok': True, 'failed': 0, 'file': str(output_file), 'skipped': True}

        logger.info(f"Creating chunk {chunk_number}: {len(chunk)} spectra -> {output_file.name}")

        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header['COMMENT'] = f'Compiled Euclid spectra - Chunk {chunk_number}'
        primary_hdu.header['NSPECTRA'] = len(chunk)
        primary_hdu.header['CREATOR'] = 'euclidkit.core.spectra.SpectrumCompiler'
        hdul_new = fits.HDUList([primary_hdu])

        grouped_rows: Dict[str, List[Any]] = defaultdict(list)
        ra_col, dec_col = self._resolve_radec_columns(chunk)
        chunk_failed = 0
        grouped_count = 0
        for source in chunk:
            file_path = self.loader.verify_spectrum_path(
                source[datalabs_path_col],
                source[file_name_col],
            )
            if file_path is None:
                chunk_failed += 1
                continue
            grouped_rows[file_path].append(source)
            grouped_count += 1

        logger.info(
            "Chunk %d grouping active: %d spectra mapped to %d unique FITS source files",
            chunk_number, grouped_count, len(grouped_rows)
        )

        progress = tqdm(total=grouped_count, desc=f"Chunk {chunk_number}") if show_progress else None
        for file_path, sources_in_file in grouped_rows.items():
            try:
                hdu_indices = [int(src[hdu_index_col]) for src in sources_in_file]
                spectra_by_index = self.loader.load_spectra_batch(file_path, hdu_indices)
            except Exception as e:
                logger.warning("Failed to load grouped spectra from %s: %s", file_path, e)
                chunk_failed += len(sources_in_file)
                if progress is not None:
                    progress.update(len(sources_in_file))
                continue

            for source in sources_in_file:
                try:
                    spectrum_hdu = spectra_by_index[int(source[hdu_index_col])]
                    spectrum_hdu.name = str(source[source_id_col])
                    spectrum_hdu.header['SOURC_ID'] = source[source_id_col]

                    if ra_col is not None:
                        spectrum_hdu.header['RA'] = source[ra_col]
                    if dec_col is not None:
                        spectrum_hdu.header['DEC'] = source[dec_col]

                    hdul_new.append(spectrum_hdu)
                except Exception as e:
                    logger.warning(f"Failed to add spectrum {source[source_id_col]}: {e}")
                    chunk_failed += 1
                finally:
                    if progress is not None:
                        progress.update(1)

        if progress is not None:
            progress.close()

        try:
            hdul_new.writeto(output_file, overwrite=overwrite)
            logger.info(f"Wrote {len(hdul_new) - 1} spectra to {output_file.name}")
            if chunk_failed > 0:
                logger.warning(f"Failed to include {chunk_failed} spectra in chunk {chunk_number}")
            return {'ok': True, 'failed': chunk_failed, 'file': str(output_file), 'skipped': False}
        except Exception as e:
            logger.error(f"Failed to write chunk {chunk_number}: {e}")
            if output_file.exists():
                output_file.unlink()
            return {'ok': False, 'failed': chunk_failed, 'file': str(output_file), 'skipped': False}
        finally:
            hdul_new.close()

    @staticmethod
    def _parse_chunk_number(file_name: str, output_prefix: str) -> Optional[int]:
        """Extract chunk number from '<prefix>_chunk_NNN.fits' file names."""
        pattern = rf"^{re.escape(output_prefix)}_chunk_(\d+)\.fits$"
        match = re.match(pattern, file_name)
        if not match:
            return None
        return int(match.group(1))

    @staticmethod
    def _count_chunk_extensions(file_path: Union[str, Path]) -> int:
        """Count spectra extensions in a compiled chunk FITS file."""
        with fits.open(file_path, memmap=True, lazy_load_hdus=True) as hdul:
            return max(0, len(hdul) - 1)

    def _discover_existing_chunks(
        self,
        output_dir: Path,
        output_prefix: str,
    ) -> List[Dict[str, Any]]:
        """
        Discover contiguous existing chunk files (starting from chunk 1).

        Returns a list of dictionaries with keys:
        - chunk_number
        - file
        - nspectra
        """
        discovered: List[Dict[str, Any]] = []
        all_chunk_files: List[tuple[int, Path]] = []

        for file_path in output_dir.glob(f"{output_prefix}_chunk_*.fits"):
            chunk_num = self._parse_chunk_number(file_path.name, output_prefix)
            if chunk_num is None:
                continue
            all_chunk_files.append((chunk_num, file_path))

        all_chunk_files.sort(key=lambda x: x[0])
        expected = 1
        for chunk_num, file_path in all_chunk_files:
            if chunk_num != expected:
                break
            try:
                nspectra = self._count_chunk_extensions(file_path)
            except Exception as e:
                logger.warning("Could not inspect existing chunk %s: %s", file_path.name, e)
                break

            discovered.append({
                'chunk_number': chunk_num,
                'file': str(file_path),
                'nspectra': int(nspectra),
            })
            expected += 1

        return discovered
    
    def compile_spectra(
        self,
        spectra_table: Table,
        output_dir: Union[str, Path],
        output_prefix: str = "compiled_spectra",
        source_id_col: str = 'source_id',
        datalabs_path_col: str = 'datalabs_path',
        file_name_col: str = 'file_name', 
        hdu_index_col: str = 'hdu_index',
        overwrite: bool = True,
        workers: int = 1,
    ) -> List[str]:
        """
        Compile spectra into chunked FITS files.
        
        Parameters
        ----------
        spectra_table : Table
            Table with spectral source information
        output_dir : str or Path
            Output directory for compiled FITS files
        output_prefix : str, default "compiled_spectra"
            Prefix for output file names
        source_id_col : str, default 'source_id'
            Column name for source IDs
        datalabs_path_col : str, default 'datalabs_path'
            Column name for datalabs paths
        file_name_col : str, default 'file_name'
            Column name for file names
        hdu_index_col : str, default 'hdu_index'
            Column name for HDU indices
        overwrite : bool, default True
            Whether to overwrite existing files. If False, existing chunk
            files are kept and skipped (resume mode).
        workers : int, default 1
            Number of chunk workers. Each worker handles different chunks.
            
        Returns
        -------
        List[str]
            List of output file paths
        """
        output_dir = Path(output_dir)
        ensure_dir(output_dir)

        source_id_col = self._resolve_source_id_column(
            spectra_table,
            preferred=source_id_col,
        )
        
        # Validate input table
        required_cols = [source_id_col, datalabs_path_col, file_name_col, hdu_index_col]
        missing_cols = [col for col in required_cols if col not in spectra_table.colnames]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        logger.info(f"Compiling {len(spectra_table)} spectra into FITS files")
        logger.info(f"Max extensions per file: {self.max_extensions}")
        if workers < 1:
            raise ValueError("workers must be >= 1")
        
        output_files: List[str] = []
        failed_count = 0
        resume_rows = 0
        next_chunk_number = 1

        if not overwrite:
            existing_chunks = self._discover_existing_chunks(output_dir, output_prefix)
            if existing_chunks:
                output_files.extend([c['file'] for c in existing_chunks])
                resume_rows = int(sum(c['nspectra'] for c in existing_chunks))
                next_chunk_number = int(existing_chunks[-1]['chunk_number']) + 1

                existing_chunk_sizes = [c['nspectra'] for c in existing_chunks if c['nspectra'] > 0]
                if len(existing_chunk_sizes) > 1:
                    reference_size = existing_chunk_sizes[0]
                    if reference_size != self.max_extensions:
                        logger.info(
                            "Existing chunk size (%d) differs from current max_extensions (%d); "
                            "resume will skip %d already compiled rows and continue.",
                            reference_size, self.max_extensions, resume_rows
                        )

                logger.info(
                    "Resume mode: found %d existing contiguous chunk files with %d compiled spectra",
                    len(existing_chunks), resume_rows
                )

        chunk_specs = []
        if resume_rows < len(spectra_table):
            remaining_table = spectra_table[resume_rows:]
            chunk_number = next_chunk_number
            for i in range(0, len(remaining_table), self.max_extensions):
                chunk = remaining_table[i:i + self.max_extensions]
                output_file = output_dir / f"{output_prefix}_chunk_{chunk_number:03d}.fits"
                chunk_specs.append((chunk_number, chunk, output_file))
                chunk_number += 1

        if workers == 1:
            for cnum, chunk, output_file in chunk_specs:
                result = self._compile_single_chunk_local(
                    chunk_number=cnum,
                    chunk=chunk,
                    output_file=output_file,
                    source_id_col=source_id_col,
                    datalabs_path_col=datalabs_path_col,
                    file_name_col=file_name_col,
                    hdu_index_col=hdu_index_col,
                    overwrite=overwrite,
                    show_progress=True,
                )
                failed_count += int(result['failed'])
                if result['ok']:
                    output_files.append(result['file'])
        elif chunk_specs:
            logger.info("Processing %d chunks with %d workers", len(chunk_specs), workers)
            with ThreadPoolExecutor(max_workers=workers) as executor:
                future_to_chunk = {
                    executor.submit(
                        self._compile_single_chunk_local,
                        chunk_number=cnum,
                        chunk=chunk,
                        output_file=output_file,
                        source_id_col=source_id_col,
                        datalabs_path_col=datalabs_path_col,
                        file_name_col=file_name_col,
                        hdu_index_col=hdu_index_col,
                        overwrite=overwrite,
                        show_progress=False,
                    ): cnum
                    for cnum, chunk, output_file in chunk_specs
                }

                progress = tqdm(total=len(chunk_specs), desc="Chunks")
                results_by_chunk: Dict[int, Dict[str, Any]] = {}
                for future in as_completed(future_to_chunk):
                    cnum = future_to_chunk[future]
                    try:
                        result = future.result()
                    except Exception as e:
                        logger.error("Chunk %d failed with exception: %s", cnum, e)
                        result = {'ok': False, 'failed': 0, 'file': '', 'skipped': False}
                    results_by_chunk[cnum] = result
                    failed_count += int(result['failed'])
                    progress.update(1)
                progress.close()

            for cnum, _, _ in chunk_specs:
                result = results_by_chunk.get(cnum, {'ok': False, 'file': ''})
                if result.get('ok'):
                    output_files.append(result['file'])
        
        logger.info(f"Compilation complete: {len(output_files)} files created")
        if failed_count > 0:
            logger.warning(f"Total failed spectra: {failed_count}")
        
        return output_files
    
    def create_metadata_table(
        self,
        spectra_table: Table,
        output_files: List[str],
        output_dir: Union[str, Path],
        output_name: str = "spectra_metadata.fits"
    ) -> str:
        """
        Create metadata table for compiled spectra.
        
        Parameters
        ----------
        spectra_table : Table
            Original spectral source table
        output_files : List[str]
            List of output FITS files
        output_dir : str or Path
            Output directory
        output_name : str, default "spectra_metadata.fits"
            Metadata file name
            
        Returns
        -------
        str
            Path to metadata file
        """
        output_path = Path(output_dir) / output_name
        
        # Create extended metadata table
        metadata = spectra_table.copy()
        
        # Add chunk information by scanning output files in sequence. This
        # supports resume scenarios where older chunk sizes differ.
        chunk_numbers: List[int] = [0] * len(metadata)
        chunk_files: List[str] = [""] * len(metadata)
        extension_numbers: List[int] = [0] * len(metadata)

        row_cursor = 0
        for order_idx, chunk_file in enumerate(output_files, start=1):
            file_name = Path(chunk_file).name
            parsed_chunk_number = self._parse_chunk_number(file_name, file_name.split("_chunk_")[0] if "_chunk_" in file_name else "")
            chunk_number = int(parsed_chunk_number) if parsed_chunk_number is not None else int(order_idx)
            try:
                nspectra = self._count_chunk_extensions(chunk_file)
            except Exception as e:
                logger.warning("Could not inspect output chunk %s for metadata: %s", file_name, e)
                nspectra = 0

            for ext_idx in range(nspectra):
                if row_cursor >= len(metadata):
                    break
                chunk_numbers[row_cursor] = chunk_number
                chunk_files[row_cursor] = file_name
                extension_numbers[row_cursor] = ext_idx + 1
                row_cursor += 1

        # Add new columns with explicit FITS-friendly dtypes.
        metadata['chunk_number'] = np.asarray(chunk_numbers, dtype=np.int32)
        metadata['chunk_file'] = np.asarray(chunk_files, dtype='U256')
        metadata['extension_number'] = np.asarray(extension_numbers, dtype=np.int32)
        
        # Save metadata
        metadata.write(output_path, format='fits', overwrite=True)
        logger.info(f"Saved metadata table to {output_path}")
        
        return str(output_path)
    
    def compile_single_fits(
        self,
        spectra_table: Table,
        output_file: Union[str, Path],
        source_id_col: str = 'source_id',
        datalabs_path_col: str = 'datalabs_path',
        file_name_col: str = 'file_name',
        hdu_index_col: str = 'hdu_index',
        overwrite: bool = True
    ) -> str:
        """
        Compile spectra into a single FITS file.
        
        This is a convenience method for creating a single FITS file
        containing all spectra, similar to the notebook cell 23 pattern.
        
        Parameters
        ----------
        spectra_table : Table
            Table with spectral source information
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
        overwrite : bool, default True
            Whether to overwrite existing files
            
        Returns
        -------
        str
            Path to output file
        """
        output_path = Path(output_file)
        output_dir = output_path.parent
        output_stem = output_path.stem
        
        # Temporarily set max_extensions to accommodate all spectra
        original_max = self.max_extensions
        self.max_extensions = len(spectra_table) + 1000
        
        try:
            output_files = self.compile_spectra(
                spectra_table,
                output_dir=output_dir,
                output_prefix=output_stem,
                source_id_col=source_id_col,
                datalabs_path_col=datalabs_path_col,
                file_name_col=file_name_col,
                hdu_index_col=hdu_index_col,
                overwrite=overwrite
            )
            
            if len(output_files) == 1:
                # Rename the chunked file to the requested name
                chunked_file = Path(output_files[0])
                if chunked_file != output_path:
                    chunked_file.rename(output_path)
                return str(output_path)
            else:
                logger.warning(f"Expected single file but got {len(output_files)} files")
                return output_files[0] if output_files else None
                
        finally:
            # Restore original max_extensions
            self.max_extensions = original_max

    def compile_spectra_datalink(
        self,
        spectra_table: Table,
        euclid_client: Any,
        output_dir: Union[str, Path],
        output_prefix: str = "compiled_spectra",
        source_id_col: str = 'source_id',
        retrieval_type: str = 'SPECTRA_RGS',
        schema: str = 'sedm',
        overwrite: bool = True,
    ) -> List[str]:
        """
        Compile spectra by retrieving each source via Euclid datalink.

        If ``overwrite`` is False, existing chunk files are kept and skipped
        (resume mode).
        """
        output_dir = Path(output_dir)
        ensure_dir(output_dir)

        source_id_col = self._resolve_source_id_column(
            spectra_table,
            preferred=source_id_col,
        )

        dedup_table, dedup_stats = self.deduplicate_by_source_id(
            spectra_table=spectra_table,
            source_id_col=source_id_col,
        )
        if dedup_stats["duplicate_rows_removed"] > 0:
            logger.info(
                "Datalink deduplication: %d input rows -> %d unique sources "
                "(removed %d duplicate rows)",
                dedup_stats["input_rows"],
                dedup_stats["unique_sources"],
                dedup_stats["duplicate_rows_removed"],
            )
        spectra_table = dedup_table

        output_files = []
        chunk_number = 1
        failed_count = 0
        ra_col, dec_col = self._resolve_radec_columns(spectra_table)

        for i in range(0, len(spectra_table), self.max_extensions):
            chunk = spectra_table[i:i + self.max_extensions]
            output_file = output_dir / f"{output_prefix}_chunk_{chunk_number:03d}.fits"
            output_files.append(str(output_file))

            if output_file.exists() and not overwrite:
                logger.info(f"Skipping existing chunk (resume mode): {output_file.name}")
                chunk_number += 1
                continue

            primary_hdu = fits.PrimaryHDU()
            primary_hdu.header['COMMENT'] = f'Compiled Euclid spectra (datalink) - Chunk {chunk_number}'
            primary_hdu.header['NSPECTRA'] = len(chunk)
            primary_hdu.header['CREATOR'] = 'euclidkit.core.spectra.SpectrumCompiler'
            hdul_new = fits.HDUList([primary_hdu])

            cache: Dict[tuple[str, str], 'ExtensionHDU'] = {}
            chunk_failed = 0
            added_count = 0

            retrieval_upper = str(retrieval_type).upper()
            retrieval_sequence = (
                ['SPECTRA_RGS', 'SPECTRA_BGS']
                if retrieval_upper == 'ALL'
                else [retrieval_upper]
            )

            for source in tqdm(chunk, desc=f"Chunk {chunk_number}"):
                sid = str(source[source_id_col])
                for retrieval_arm in retrieval_sequence:
                    try:
                        cache_key = (sid, retrieval_arm)
                        if cache_key not in cache:
                            cache[cache_key] = self.loader.load_spectrum_from_datalink(
                                euclid_client=euclid_client,
                                source_id=sid,
                                schema=schema,
                                retrieval_type=retrieval_arm,
                            )
                        spectrum_hdu = cache[cache_key].copy()
                        if len(retrieval_sequence) > 1:
                            arm_tag = 'RGS' if retrieval_arm == 'SPECTRA_RGS' else 'BGS'
                            spectrum_hdu.name = f"{sid}_{arm_tag}"
                        else:
                            spectrum_hdu.name = sid
                        spectrum_hdu.header['SOURC_ID'] = sid
                        spectrum_hdu.header['LRANGE'] = 'RGS' if retrieval_arm == 'SPECTRA_RGS' else 'BGS'
                        spectrum_hdu.header['RETRTYPE'] = retrieval_arm

                        if ra_col is not None:
                            spectrum_hdu.header['RA'] = source[ra_col]
                        if dec_col is not None:
                            spectrum_hdu.header['DEC'] = source[dec_col]

                        hdul_new.append(spectrum_hdu)
                        added_count += 1
                    except Exception as e:
                        logger.warning(
                            "Failed to add datalink spectrum %s (%s): %s",
                            sid,
                            retrieval_arm,
                            e,
                        )
                        chunk_failed += 1

            primary_hdu.header['NSPECTRA'] = added_count

            try:
                hdul_new.writeto(output_file, overwrite=overwrite)
                if chunk_failed > 0:
                    failed_count += chunk_failed
            except Exception as e:
                logger.error(f"Failed to write chunk {chunk_number}: {e}")
                if output_file.exists():
                    output_file.unlink()
                output_files.pop()
            finally:
                hdul_new.close()

            chunk_number += 1

        if failed_count > 0:
            logger.warning(f"Total failed spectra (datalink): {failed_count}")
        return output_files


class SpectrumProcessor:
    """Processor for spectral analysis and manipulation."""
    
    def __init__(self):
        pass
    
    def validate_spectrum(self, spectrum_hdu: 'ExtensionHDU') -> Dict[str, Any]:
        """
        Validate a spectrum HDU.
        
        Parameters
        ----------
        spectrum_hdu : ExtensionHDU
            Spectrum HDU to validate
            
        Returns
        -------
        Dict[str, Any]
            Validation results
        """
        validation = {
            'valid': True,
            'issues': [],
            'n_pixels': 0,
            'wavelength_range': None,
            'flux_range': None
        }
        
        try:
            data = spectrum_hdu.data
            if data is None:
                validation['valid'] = False
                validation['issues'].append('No data in HDU')
                return validation
            
            # Check for required columns
            required_cols = ['WAVELENGTH', 'SIGNAL']  # Standard Euclid spectrum columns
            missing_cols = []
            
            if hasattr(data, 'columns'):  # Binary table
                existing_cols = [col.name for col in data.columns]
                missing_cols = [col for col in required_cols if col not in existing_cols]
                if missing_cols:
                    validation['issues'].append(f'Missing columns: {missing_cols}')
                else:
                    validation['n_pixels'] = len(data)
                    validation['wavelength_range'] = (
                        float(np.min(data['WAVELENGTH'])), 
                        float(np.max(data['WAVELENGTH']))
                    )
                    validation['flux_range'] = (
                        float(np.min(data['SIGNAL'])),
                        float(np.max(data['SIGNAL']))
                    )
            else:
                validation['issues'].append('Unexpected data format')
            
            if validation['issues']:
                validation['valid'] = False
                
        except Exception as e:
            validation['valid'] = False
            validation['issues'].append(f'Validation error: {str(e)}')
        
        return validation
