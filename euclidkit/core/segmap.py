"""Local MER segmentation-map cutout compilation utilities."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata.utils import Cutout2D
from astropy.table import Table
from astropy.wcs import WCS
from tqdm import tqdm

from euclidkit.utils.io import load_table

logger = logging.getLogger(__name__)


@dataclass
class SegmapCutoutStats:
    """Summary returned by ``compile_segmap_cutouts``."""

    requested_rows: int = 0
    written_rows: int = 0
    failed_rows: int = 0
    skipped_existing_rows: int = 0
    output_files: List[str] = field(default_factory=list)
    failures: List[Dict[str, Any]] = field(default_factory=list)

    @property
    def output_file_count(self) -> int:
        """Number of cutout FITS files written in this run."""
        return len(self.output_files)


@dataclass(frozen=True)
class _SegmapRow:
    row_index: int
    object_id: Any
    segmentation_map_id: Any
    ra: float
    dec: float
    tile_index: Optional[Any]
    file_path: Path
    file_name: str
    output_path: Path


def _column_lookup(table: Table) -> Dict[str, str]:
    return {str(name).lower(): str(name) for name in table.colnames}


def _resolve_column(table: Table, *candidates: str) -> Optional[str]:
    lookup = _column_lookup(table)
    for candidate in candidates:
        actual = lookup.get(candidate.lower())
        if actual is not None:
            return actual
    return None


def _is_missing(value: Any) -> bool:
    if np.ma.is_masked(value) or value is None:
        return True
    try:
        missing = np.isnan(value)
    except (TypeError, ValueError):
        return False
    try:
        if hasattr(missing, "item"):
            missing = missing.item()
        return bool(missing)
    except (TypeError, ValueError):
        return False


def _as_clean_text(value: Any) -> str:
    if hasattr(value, "item"):
        try:
            value = value.item()
        except Exception:
            pass
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    return str(value).strip()


def _safe_filename_token(value: Any) -> str:
    text = _as_clean_text(value)
    text = re.sub(r"[^0-9A-Za-z_.-]+", "-", text).strip("-._")
    return text or "unknown"


def _header_value(value: Any) -> Any:
    if hasattr(value, "item"):
        try:
            value = value.item()
        except Exception:
            pass
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    if isinstance(value, (np.integer, int)):
        return int(value) if abs(int(value)) <= 2_147_483_647 else str(int(value))
    if isinstance(value, (np.floating, float)):
        return float(value)
    return str(value)


def _required_columns(table: Table) -> Tuple[str, str, str, str, str, str, Optional[str]]:
    datalabs_path_col = _resolve_column(table, "datalabs_path")
    file_name_col = _resolve_column(table, "file_name")
    seg_col = _resolve_column(table, "segmentation_map_id", "SEGMENTATION_MAP_ID")
    object_id_col = _resolve_column(table, "object_id")
    ra_col = _resolve_column(table, "ra", "mer_ra")
    dec_col = _resolve_column(table, "dec", "mer_dec")
    tile_col = _resolve_column(table, "tile_index")

    missing = []
    if datalabs_path_col is None:
        missing.append("datalabs_path")
    if file_name_col is None:
        missing.append("file_name")
    if seg_col is None:
        missing.append("segmentation_map_id")
    if object_id_col is None:
        missing.append("object_id")
    if ra_col is None:
        missing.append("ra or mer_ra")
    if dec_col is None:
        missing.append("dec or mer_dec")
    if missing:
        raise ValueError("Missing required column(s) for compile-segmap: " + ", ".join(missing))

    return (
        datalabs_path_col,
        file_name_col,
        seg_col,
        object_id_col,
        ra_col,
        dec_col,
        tile_col,
    )


def _prepare_rows(table: Table, output_dir: Path) -> List[_SegmapRow]:
    (
        datalabs_path_col,
        file_name_col,
        seg_col,
        object_id_col,
        ra_col,
        dec_col,
        tile_col,
    ) = _required_columns(table)

    rows: List[_SegmapRow] = []
    for row_index, row in enumerate(table):
        datalabs_path = _as_clean_text(row[datalabs_path_col])
        file_name = _as_clean_text(row[file_name_col])
        if not datalabs_path or not file_name:
            raise ValueError(f"Row {row_index} has empty datalabs_path or file_name")
        ra = float(row[ra_col])
        dec = float(row[dec_col])
        object_id = row[object_id_col]
        seg_id = row[seg_col]
        tile_index = row[tile_col] if tile_col is not None and not _is_missing(row[tile_col]) else None
        output_name = (
            f"segmap_cutout_row{row_index:06d}_"
            f"obj{_safe_filename_token(object_id)}_"
            f"seg{_safe_filename_token(seg_id)}.fits"
        )
        rows.append(_SegmapRow(
            row_index=row_index,
            object_id=object_id,
            segmentation_map_id=seg_id,
            ra=ra,
            dec=dec,
            tile_index=tile_index,
            file_path=Path(datalabs_path) / file_name,
            file_name=file_name,
            output_path=output_dir / output_name,
        ))
    return rows


def _write_cutout(row: _SegmapRow, hdul: fits.HDUList, size_arcsec: float, overwrite: bool) -> None:
    hdu = hdul[0]
    if hdu.data is None:
        raise ValueError(f"Segmentation map has no primary image data: {row.file_path}")

    image_wcs = WCS(hdu.header)
    position = SkyCoord(row.ra * u.deg, row.dec * u.deg, frame="icrs")
    cutout = Cutout2D(
        hdu.data,
        position=position,
        size=(size_arcsec * u.arcsec, size_arcsec * u.arcsec),
        wcs=image_wcs,
        mode="partial",
        fill_value=0,
        copy=True,
    )

    out_hdu = fits.PrimaryHDU(np.asarray(cutout.data))
    out_hdu.header.update(cutout.wcs.to_header())
    out_hdu.header["OBJECTID"] = (_header_value(row.object_id), "Input object_id")
    out_hdu.header["SEGMAPID"] = (_header_value(row.segmentation_map_id), "Input segmentation_map_id")
    out_hdu.header["SRC_RA"] = (float(row.ra), "Source RA, deg")
    out_hdu.header["SRC_DEC"] = (float(row.dec), "Source Dec, deg")
    if row.tile_index is not None:
        out_hdu.header["TILEIND"] = (_header_value(row.tile_index), "Segmentation tile index")
    out_hdu.header["SRCFILE"] = _as_clean_text(row.file_name)[:68]
    out_hdu.header["CUTSIZE"] = (float(size_arcsec), "Cutout size, arcsec")
    row.output_path.parent.mkdir(parents=True, exist_ok=True)
    out_hdu.writeto(row.output_path, overwrite=overwrite)


def compile_segmap_cutouts(
    input_table: Union[str, Path, Table],
    output_dir: Union[str, Path],
    size_arcsec: float = 10.0,
    overwrite: bool = False,
    on_error: str = "fail",
    show_progress: bool = True,
) -> SegmapCutoutStats:
    """Create FITS cutouts from local MER segmentation-map files."""
    if size_arcsec <= 0:
        raise ValueError("size_arcsec must be positive")
    if on_error not in {"fail", "skip"}:
        raise ValueError("on_error must be one of: fail, skip")

    table = load_table(input_table) if not isinstance(input_table, Table) else input_table
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    rows = _prepare_rows(table, output_path)
    stats = SegmapCutoutStats(requested_rows=len(rows))

    rows_by_file: Dict[Path, List[_SegmapRow]] = {}
    for row in rows:
        rows_by_file.setdefault(row.file_path, []).append(row)

    file_iter: Iterable[Tuple[Path, List[_SegmapRow]]] = rows_by_file.items()
    if show_progress:
        file_iter = tqdm(file_iter, total=len(rows_by_file), desc="Segmentation maps")

    for file_path, file_rows in file_iter:
        pending_rows: List[_SegmapRow] = []
        for row in file_rows:
            if row.output_path.exists() and not overwrite:
                stats.skipped_existing_rows += 1
            else:
                pending_rows.append(row)
        if not pending_rows:
            continue

        try:
            if not file_path.exists():
                raise FileNotFoundError(f"Segmentation map file not found: {file_path}")
            with fits.open(file_path, memmap=True) as hdul:
                for row in pending_rows:
                    try:
                        _write_cutout(row, hdul, size_arcsec=size_arcsec, overwrite=overwrite)
                        stats.written_rows += 1
                        stats.output_files.append(str(row.output_path))
                    except Exception as exc:
                        failure = {
                            "row_index": row.row_index,
                            "file_path": str(row.file_path),
                            "output_path": str(row.output_path),
                            "error": str(exc),
                        }
                        stats.failures.append(failure)
                        stats.failed_rows += 1
                        if on_error == "fail":
                            raise
                        logger.warning("Skipping segmentation cutout row %s: %s", row.row_index, exc)
        except Exception as exc:
            if on_error == "fail":
                raise
            for row in pending_rows:
                failure = {
                    "row_index": row.row_index,
                    "file_path": str(row.file_path),
                    "output_path": str(row.output_path),
                    "error": str(exc),
                }
                stats.failures.append(failure)
                stats.failed_rows += 1
            logger.warning("Skipping segmentation map file %s: %s", file_path, exc)

    return stats
