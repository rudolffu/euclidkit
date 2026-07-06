"""Datalabs-local FITS spectra export to parquet.

This module reads SIR CombinedSpectra FITS files listed by catalog tables
(`datalabs_path`, `file_name`, `hdu_index`) and writes row-wise parquet parts
for downstream ML workflows. It intentionally does not use Euclid Datalink.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import sys
from typing import Any

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table


@dataclass
class RawParquetStats:
    requested_rows: int = 0
    exported_rows: int = 0
    skipped_rows: int = 0
    failed_rows: int = 0
    lambda_skipped_rows: int = 0
    file_missing_rows: int = 0
    output_files: list[str] = field(default_factory=list)
    manifest_path: str = ""
    failures_path: str = ""


@dataclass
class DithersParquetStats:
    objects_requested: int = 0
    objects_exported: int = 0
    skipped_rows: int = 0
    failed_rows: int = 0
    lambda_skipped_rows: int = 0
    file_missing_rows: int = 0
    objects_without_dithers: int = 0
    combined_rows: int = 0
    dither_rows: int = 0
    combined_output_files: list[str] = field(default_factory=list)
    dither_output_files: list[str] = field(default_factory=list)
    manifest_path: str = ""
    failures_path: str = ""
    raw_frame_table: str = ""
    raw_frame_rows: int = 0
    raw_frame_spectroimage_rows: int = 0
    dither_rows_with_raw_frame_metadata: int = 0
    dither_rows_missing_raw_frame_metadata: int = 0


def default_raw_parquet_workers() -> int:
    """Default process workers for Datalabs-local parquet export."""
    return min(os.cpu_count() or 1, 8)


def spectra_to_parquet(
    catalog_table: str,
    output_prefix: str,
    *,
    chunk_size: int = 2000,
    workers: int | None = None,
    lambda_range: str | None = None,
    overwrite: bool = False,
    on_error: str = "fail",
    show_progress: bool = False,
    limit: int | None = None,
) -> RawParquetStats:
    """Export catalog-listed combined spectra HDUs to parquet parts."""
    if chunk_size <= 0:
        raise ValueError("--chunk-size must be positive")
    if workers is None:
        workers = default_raw_parquet_workers()
    if workers <= 0:
        raise ValueError("--workers must be positive")
    if on_error not in {"fail", "skip"}:
        raise ValueError("--on-error must be 'fail' or 'skip'")
    if limit is not None and limit <= 0:
        raise ValueError("--limit must be a positive integer")

    target_lrange = str(lambda_range).upper() if lambda_range else None
    if target_lrange not in {None, "RGS", "BGS"}:
        raise ValueError("--lambda-range must be RGS or BGS")

    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    manifest_path = prefix.with_name(f"{prefix.name}_manifest.json")
    failures_path = prefix.with_name(f"{prefix.name}_failures.jsonl")
    _prepare_raw_output(prefix, manifest_path, failures_path, overwrite=overwrite)

    catalog = _read_catalog_table(Path(catalog_table))
    _validate_catalog(catalog)
    if limit is not None:
        catalog = catalog.iloc[:limit].copy()
    catalog = catalog.reset_index(drop=False).rename(columns={"index": "__row_index"})
    stats = RawParquetStats(requested_rows=int(len(catalog)), manifest_path=str(manifest_path))

    part_idx = 1
    failures: list[dict[str, Any]] = []

    for start in range(0, len(catalog), chunk_size):
        chunk = catalog.iloc[start : start + chunk_size].copy()
        tasks = [_catalog_row_to_task(row, target_lrange) for row in chunk.to_dict(orient="records")]
        results = _run_raw_tasks(tasks, workers)

        rows: list[dict[str, Any]] = []
        for result in sorted(results, key=lambda item: int(item["row_index"])):
            status = result["status"]
            if status == "ok":
                rows.append(result["row"])
            elif status == "lambda_skipped":
                stats.skipped_rows += 1
                stats.lambda_skipped_rows += 1
            elif status == "file_missing":
                failures.append(_failure_record(result))
                stats.skipped_rows += 1
                stats.file_missing_rows += 1
            else:
                failure = _failure_record(result)
                if on_error == "fail":
                    raise RuntimeError(
                        "failed to read catalog row {row_index} from {file_path} hdu={hdu_index}: {error}".format(
                            **failure
                        )
                    )
                failures.append(failure)
                stats.skipped_rows += 1
                stats.failed_rows += 1

        if rows:
            frame = pd.DataFrame(rows)
            part_path = prefix.with_name(f"{prefix.name}_part{part_idx:03d}.parquet")
            _write_parquet_atomic(frame, part_path)
            stats.output_files.append(str(part_path))
            stats.exported_rows += int(len(frame))
            part_idx += 1

        _progress(show_progress, min(start + len(chunk), len(catalog)), len(catalog), "spectra-to-parquet")

    if failures:
        _write_failures(failures_path, failures)
        stats.failures_path = str(failures_path)

    _write_manifest(
        manifest_path,
        {
            "created_at": _utc_now(),
            "export_type": "raw_spectra_parquet",
            "catalog_table": str(Path(catalog_table)),
            "requested_rows": stats.requested_rows,
            "exported_rows": stats.exported_rows,
            "skipped_rows": stats.skipped_rows,
            "failed_rows": stats.failed_rows,
            "lambda_skipped_rows": stats.lambda_skipped_rows,
            "file_missing_rows": stats.file_missing_rows,
            "chunk_size": int(chunk_size),
            "workers": int(workers),
            "lambda_range": target_lrange,
            "on_error": on_error,
            "limit": limit,
            "output_files": stats.output_files,
            "failures_path": stats.failures_path,
        },
    )
    if show_progress and stats.requested_rows:
        print(file=sys.stderr)
    return stats


def dithers_to_parquet(
    catalog_table: str,
    output_prefix: str,
    *,
    chunk_size: int = 2000,
    workers: int | None = None,
    lambda_range: str | None = None,
    include_combined: bool = True,
    overwrite: bool = False,
    on_error: str = "fail",
    show_progress: bool = False,
    raw_frame_table: str | None = None,
) -> DithersParquetStats:
    """Export combined and per-dither spectra for catalog-listed objects."""
    if chunk_size <= 0:
        raise ValueError("--chunk-size must be positive")
    if workers is None:
        workers = default_raw_parquet_workers()
    if workers <= 0:
        raise ValueError("--workers must be positive")
    if on_error not in {"fail", "skip"}:
        raise ValueError("--on-error must be 'fail' or 'skip'")

    target_lrange = str(lambda_range).upper() if lambda_range else None
    if target_lrange not in {None, "RGS", "BGS"}:
        raise ValueError("--lambda-range must be RGS or BGS")

    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    manifest_path = prefix.with_name(f"{prefix.name}_manifest.json")
    failures_path = prefix.with_name(f"{prefix.name}_failures.jsonl")
    _prepare_dither_output(prefix, manifest_path, failures_path, overwrite=overwrite)

    catalog = _read_catalog_table(Path(catalog_table))
    _validate_catalog(catalog)
    catalog = catalog.reset_index(drop=False).rename(columns={"index": "__row_index"})
    stats = DithersParquetStats(objects_requested=int(len(catalog)), manifest_path=str(manifest_path))
    raw_frame_lookup: dict[int, dict[str, Any]] | None = None
    if raw_frame_table is not None:
        raw_frame_lookup, raw_frame_rows, spectroimage_rows = _build_raw_frame_lookup(Path(raw_frame_table))
        stats.raw_frame_table = str(Path(raw_frame_table))
        stats.raw_frame_rows = raw_frame_rows
        stats.raw_frame_spectroimage_rows = spectroimage_rows

    failures: list[dict[str, Any]] = []
    combined_part_idx = 1
    dither_part_idx = 1

    for start in range(0, len(catalog), chunk_size):
        chunk = catalog.iloc[start : start + chunk_size].copy()
        tasks = [
            _catalog_row_to_task(row, target_lrange, include_combined=include_combined)
            for row in chunk.to_dict(orient="records")
        ]
        results = _run_dither_tasks(tasks, workers)

        combined_rows: list[dict[str, Any]] = []
        dither_rows: list[dict[str, Any]] = []
        for result in sorted(results, key=lambda item: int(item["row_index"])):
            status = result["status"]
            if status == "ok":
                stats.objects_exported += 1
                if result.get("combined_row") is not None:
                    combined_rows.append(result["combined_row"])
                object_dithers = result.get("dither_rows", [])
                dither_rows.extend(object_dithers)
                if not object_dithers:
                    stats.objects_without_dithers += 1
            elif status == "lambda_skipped":
                stats.skipped_rows += 1
                stats.lambda_skipped_rows += 1
            elif status == "file_missing":
                failures.append(_failure_record(result))
                stats.skipped_rows += 1
                stats.file_missing_rows += 1
            else:
                failure = _failure_record(result)
                if on_error == "fail":
                    raise RuntimeError(
                        "failed to read catalog row {row_index} from {file_path} hdu={hdu_index}: {error}".format(
                            **failure
                        )
                    )
                failures.append(failure)
                stats.skipped_rows += 1
                stats.failed_rows += 1

        if include_combined and combined_rows:
            frame = pd.DataFrame(combined_rows)
            part_path = prefix.with_name(f"{prefix.name}_combined_part{combined_part_idx:03d}.parquet")
            _write_parquet_atomic(frame, part_path)
            stats.combined_output_files.append(str(part_path))
            stats.combined_rows += int(len(frame))
            combined_part_idx += 1
        if dither_rows:
            if raw_frame_lookup is not None:
                (
                    dither_rows,
                    matched_metadata,
                    missing_metadata,
                ) = _enrich_dither_rows_with_raw_frame(dither_rows, raw_frame_lookup)
                stats.dither_rows_with_raw_frame_metadata += matched_metadata
                stats.dither_rows_missing_raw_frame_metadata += missing_metadata
            frame = pd.DataFrame(dither_rows)
            part_path = prefix.with_name(f"{prefix.name}_dithers_part{dither_part_idx:03d}.parquet")
            _write_parquet_atomic(frame, part_path)
            stats.dither_output_files.append(str(part_path))
            stats.dither_rows += int(len(frame))
            dither_part_idx += 1

        _progress(show_progress, min(start + len(chunk), len(catalog)), len(catalog), "dithers-to-parquet")

    if failures:
        _write_failures(failures_path, failures)
        stats.failures_path = str(failures_path)

    _write_manifest(
        manifest_path,
        {
            "created_at": _utc_now(),
            "export_type": "sir_dithers_parquet",
            "catalog_table": str(Path(catalog_table)),
            "objects_requested": stats.objects_requested,
            "objects_exported": stats.objects_exported,
            "skipped_rows": stats.skipped_rows,
            "failed_rows": stats.failed_rows,
            "lambda_skipped_rows": stats.lambda_skipped_rows,
            "file_missing_rows": stats.file_missing_rows,
            "objects_without_dithers": stats.objects_without_dithers,
            "combined_rows": stats.combined_rows,
            "dither_rows": stats.dither_rows,
            "chunk_size": int(chunk_size),
            "workers": int(workers),
            "lambda_range": target_lrange,
            "include_combined": bool(include_combined),
            "on_error": on_error,
            "combined_output_files": stats.combined_output_files,
            "dither_output_files": stats.dither_output_files,
            "failures_path": stats.failures_path,
            "raw_frame_table": stats.raw_frame_table,
            "raw_frame_rows": stats.raw_frame_rows,
            "raw_frame_spectroimage_rows": stats.raw_frame_spectroimage_rows,
            "dither_rows_with_raw_frame_metadata": stats.dither_rows_with_raw_frame_metadata,
            "dither_rows_missing_raw_frame_metadata": stats.dither_rows_missing_raw_frame_metadata,
        },
    )
    if show_progress and stats.objects_requested:
        print(file=sys.stderr)
    return stats


def parse_combined_extname(extname: str) -> str | None:
    match = re.match(r"^(.+?)_COMBINED1D_SIGNAL$", str(extname).strip(), flags=re.IGNORECASE)
    return match.group(1) if match else None


def parse_dither_extname(extname: str) -> tuple[str, int | None] | None:
    match = re.match(r"^(.+?)_DITH1D_(.+?)_SIGNAL$", str(extname).strip(), flags=re.IGNORECASE)
    if not match:
        return None
    try:
        dither_id = int(match.group(2))
    except ValueError:
        dither_id = None
    return match.group(1), dither_id


def _run_raw_tasks(tasks: list[dict[str, Any]], workers: int) -> list[dict[str, Any]]:
    if workers == 1:
        return [_read_raw_spectrum_task(task) for task in tasks]

    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(_read_raw_spectrum_task, task) for task in tasks]
        for future in as_completed(futures):
            results.append(future.result())
    return results


def _run_dither_tasks(tasks: list[dict[str, Any]], workers: int) -> list[dict[str, Any]]:
    if workers == 1:
        return [_read_dithers_task(task) for task in tasks]

    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(_read_dithers_task, task) for task in tasks]
        for future in as_completed(futures):
            results.append(future.result())
    return results


def _read_raw_spectrum_task(task: dict[str, Any]) -> dict[str, Any]:
    try:
        path = Path(str(task["file_path"]))
        hdu_index = int(task["hdu_index"])
        target_lrange = task.get("lambda_range")

        with fits.open(path, memmap=True) as hdul:
            primary = hdul[0].header if len(hdul) else {}
            lrange = _header_str(primary, "LRANGE")
            if target_lrange is not None and lrange != target_lrange:
                return {**task, "status": "lambda_skipped", "actual_lambda_range": lrange}

            if hdu_index <= 0 or hdu_index >= len(hdul):
                raise ValueError(f"hdu_index out of range for file with {len(hdul)} HDUs")

            hdu = hdul[hdu_index]
            row = {
                "object_id": int(task["object_id"]),
                "source_id": int(task["source_id"]),
                "ra": float(task["ra"]),
                "dec": float(task["dec"]),
                **_table_arrays(hdu),
            }

        return {**task, "status": "ok", "row": row}
    except FileNotFoundError as exc:
        return {**task, "status": "file_missing", "error": str(exc)}
    except Exception as exc:
        return {**task, "status": "failed", "error": str(exc)}


def _read_dithers_task(task: dict[str, Any]) -> dict[str, Any]:
    try:
        path = Path(str(task["file_path"]))
        hdu_index = int(task["hdu_index"])
        target_lrange = task.get("lambda_range")
        include_combined = bool(task["include_combined"])

        with fits.open(path, memmap=True) as hdul:
            primary = hdul[0].header if len(hdul) else {}
            lrange = _header_str(primary, "LRANGE")
            if target_lrange is not None and lrange != target_lrange:
                return {**task, "status": "lambda_skipped", "actual_lambda_range": lrange}

            if hdu_index <= 0 or hdu_index >= len(hdul):
                raise ValueError(f"hdu_index out of range for file with {len(hdul)} HDUs")

            combined_hdu = hdul[hdu_index]
            combined_extname = str(combined_hdu.header.get("EXTNAME", "")).strip()
            object_prefix = parse_combined_extname(combined_extname)
            if object_prefix is None:
                raise ValueError(f"anchor HDU is not COMBINED1D signal: {combined_extname}")

            combined_row = None
            if include_combined:
                combined_row = _base_row(task, path, lrange)
                combined_row.update(
                    {
                        "combined_hdu_index": hdu_index,
                        **_header_fields(combined_hdu.header),
                        **_table_arrays(combined_hdu),
                    }
                )

            dither_rows: list[dict[str, Any]] = []
            for idx, hdu in enumerate(hdul[1:], start=1):
                extname = str(hdu.header.get("EXTNAME", "")).strip()
                parsed = parse_dither_extname(extname)
                if parsed is None:
                    continue
                prefix, dither_id = parsed
                if prefix.upper() != object_prefix.upper():
                    continue
                header_dither_id = _header_int(hdu.header, "DITH_ID")
                row = _base_row(task, path, lrange)
                row.update(
                    {
                        "dither_hdu_index": idx,
                        "dither_id": header_dither_id if header_dither_id is not None else dither_id,
                        "ptgid": _header_int(hdu.header, "PTGID"),
                        "gwa_pos": _header_str(hdu.header, "GWA_POS"),
                        "gwa_tilt": _header_int(hdu.header, "GWA_TILT"),
                        "det_id": _header_int(hdu.header, "DET_ID"),
                        "det_id_2": _header_int(hdu.header, "DET_ID_2"),
                        **_header_fields(hdu.header),
                        **_table_arrays(hdu),
                    }
                )
                dither_rows.append(row)

        return {**task, "status": "ok", "combined_row": combined_row, "dither_rows": dither_rows}
    except FileNotFoundError as exc:
        return {**task, "status": "file_missing", "error": str(exc)}
    except Exception as exc:
        return {**task, "status": "failed", "error": str(exc)}


def _base_row(task: dict[str, Any], path: Path, lrange: str | None) -> dict[str, Any]:
    return {
        "object_id": int(task["object_id"]),
        "source_id": int(task["source_id"]),
        "ra": float(task["ra"]),
        "dec": float(task["dec"]),
        "input_fits": str(path),
        "catalog_row": int(task["row_index"]),
        "lambda_range": lrange,
    }


def _table_arrays(hdu: fits.BinTableHDU) -> dict[str, Any]:
    data = hdu.data
    if data is None:
        raise ValueError("requested HDU has no table data")
    names = list(data.columns.names or [])
    name_lookup = {name.upper(): name for name in names}
    if "SIGNAL" not in name_lookup:
        raise ValueError("requested HDU does not contain SIGNAL")

    row: dict[str, Any] = {}
    for name in names:
        out_name = str(name).lower()
        if out_name in {"err", "valid", "ivar"}:
            continue
        row[out_name] = _array_to_list(data[name])

    if "wavelength" not in row:
        wavelength = _wavelength_from_header(hdu.header)
        if wavelength is None:
            raise ValueError("missing WAVELENGTH column and WMIN/BINWIDTH/BINCOUNT header grid")
        row["wavelength"] = wavelength.astype(np.float64, copy=False).tolist()
    if "signal" not in row:
        raise ValueError("missing SIGNAL column")
    row["flux"] = row["signal"]
    return row


def _header_fields(header: Any) -> dict[str, Any]:
    return {
        "extname": str(header.get("EXTNAME", "")).strip(),
        "exptime": _header_float(header, "EXPTIME"),
        "wmin": _header_float(header, "WMIN"),
        "binwidth": _header_float(header, "BINWIDTH"),
        "bincount": _header_int(header, "BINCOUNT"),
        "lsf_sig": _header_float(header, "LSF_SIG"),
        "ext_prof": _header_str(header, "EXT_PROF"),
    }


def _read_catalog_table(catalog_path: Path) -> pd.DataFrame:
    name = catalog_path.name.lower()
    suffix = catalog_path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(catalog_path)
    if suffix == ".parquet":
        return pd.read_parquet(catalog_path)
    if suffix == ".feather":
        return pd.read_feather(catalog_path)
    if suffix in {".fits", ".fit"} or name.endswith((".fits.gz", ".fit.gz")):
        from astropy.table import Table

        return Table.read(catalog_path, format="fits", character_as_bytes=False).to_pandas()
    raise ValueError(f"Unsupported catalog table format: {catalog_path}")


def _read_raw_frame_table(raw_frame_path: Path) -> pd.DataFrame:
    name = raw_frame_path.name.lower()
    suffix = raw_frame_path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(raw_frame_path)
    if suffix == ".parquet":
        return pd.read_parquet(raw_frame_path)
    if suffix == ".feather":
        return pd.read_feather(raw_frame_path)
    if suffix in {".fits", ".fit"} or name.endswith((".fits.gz", ".fit.gz")):
        return Table.read(raw_frame_path, format="fits", character_as_bytes=False).to_pandas()
    if suffix in {".vot", ".xml"} or name.endswith((".vot.gz", ".xml.gz")):
        return Table.read(raw_frame_path, format="votable").to_pandas()
    try:
        return Table.read(raw_frame_path).to_pandas()
    except Exception as exc:
        raise ValueError(f"Unsupported raw-frame table format: {raw_frame_path}") from exc


def _build_raw_frame_lookup(raw_frame_path: Path) -> tuple[dict[int, dict[str, Any]], int, int]:
    frame = _read_raw_frame_table(raw_frame_path)
    row_count = int(len(frame))
    required = {"pointing_id", "technique", "obs_time_mjd", "obs_time_utc", "pa"}
    column_lookup = {str(col).lower(): col for col in frame.columns}
    missing = sorted(required.difference(column_lookup))
    if missing:
        raise ValueError(f"raw-frame table missing required columns: {missing}")

    raw = pd.DataFrame({name: frame[column_lookup[name]] for name in required})
    technique = raw["technique"].map(_as_text).str.strip().str.upper()
    raw = raw.loc[technique == "SPECTROIMAGE"].copy()
    spectroimage_count = int(len(raw))

    raw["pointing_id"] = pd.to_numeric(raw["pointing_id"], errors="coerce")
    raw["obs_time_mjd"] = pd.to_numeric(raw["obs_time_mjd"], errors="coerce")
    raw["pa"] = pd.to_numeric(raw["pa"], errors="coerce")
    raw["obs_time_utc"] = raw["obs_time_utc"].map(_clean_text_or_none)
    raw = raw.dropna(subset=["pointing_id"])

    lookup: dict[int, dict[str, Any]] = {}
    for row in raw.to_dict(orient="records"):
        pointing_id = _optional_int(row.get("pointing_id"))
        if pointing_id is None or pointing_id in lookup:
            continue
        lookup[pointing_id] = {
            "obs_time_mjd": _optional_float(row.get("obs_time_mjd")),
            "obs_time_utc": row.get("obs_time_utc"),
            "pa": _optional_float(row.get("pa")),
        }
    return lookup, row_count, spectroimage_count


def _enrich_dither_rows_with_raw_frame(
    dither_rows: list[dict[str, Any]],
    raw_frame_lookup: dict[int, dict[str, Any]],
) -> tuple[list[dict[str, Any]], int, int]:
    matched = 0
    missing = 0
    enriched: list[dict[str, Any]] = []
    for row in dither_rows:
        out = dict(row)
        pointing_id = _dither_pointing_id(out)
        metadata = raw_frame_lookup.get(pointing_id) if pointing_id is not None else None
        if metadata is not None:
            out.update(metadata)
            matched += 1
        else:
            out["obs_time_mjd"] = None
            out["obs_time_utc"] = None
            out["pa"] = None
            missing += 1
        enriched.append(out)
    return enriched, matched, missing


def _dither_pointing_id(row: dict[str, Any]) -> int | None:
    dither_id = _optional_int(row.get("dither_id"))
    if dither_id is not None:
        return dither_id
    return _optional_int(row.get("ptgid"))


def _validate_catalog(catalog: pd.DataFrame) -> None:
    required = {"datalabs_path", "file_name", "hdu_index", "source_id", "object_id", "ra_obj", "dec_obj"}
    missing = required.difference(set(catalog.columns))
    if missing:
        raise ValueError(f"catalog table missing required columns: {sorted(missing)}")


def _catalog_row_to_task(
    row: dict[str, Any],
    lambda_range: str | None,
    *,
    include_combined: bool | None = None,
) -> dict[str, Any]:
    file_name = _as_text(row["file_name"]).strip()
    datalabs_path = _as_text(row["datalabs_path"]).strip()
    if not file_name:
        raise ValueError("catalog row has empty file_name")
    if not datalabs_path:
        raise ValueError("catalog row has empty datalabs_path")
    task = {
        "row_index": int(row["__row_index"]),
        "file_path": str(Path(datalabs_path) / file_name),
        "hdu_index": int(row["hdu_index"]),
        "object_id": int(row["object_id"]),
        "source_id": int(row["source_id"]),
        "ra": float(row["ra_obj"]),
        "dec": float(row["dec_obj"]),
        "lambda_range": lambda_range,
    }
    if include_combined is not None:
        task["include_combined"] = include_combined
    return task


def _prepare_raw_output(prefix: Path, manifest_path: Path, failures_path: Path, *, overwrite: bool) -> None:
    existing = (
        list(prefix.parent.glob(f"{prefix.name}_part*.parquet"))
        + list(prefix.parent.glob(f"{prefix.name}_part*.parquet.tmp"))
        + [path for path in (manifest_path, failures_path) if path.exists()]
    )
    _clear_or_raise(existing, overwrite=overwrite)


def _prepare_dither_output(prefix: Path, manifest_path: Path, failures_path: Path, *, overwrite: bool) -> None:
    existing = (
        list(prefix.parent.glob(f"{prefix.name}_combined_part*.parquet"))
        + list(prefix.parent.glob(f"{prefix.name}_combined_part*.parquet.tmp"))
        + list(prefix.parent.glob(f"{prefix.name}_dithers_part*.parquet"))
        + list(prefix.parent.glob(f"{prefix.name}_dithers_part*.parquet.tmp"))
        + [path for path in (manifest_path, failures_path) if path.exists()]
    )
    _clear_or_raise(existing, overwrite=overwrite)


def _clear_or_raise(existing: list[Path], *, overwrite: bool) -> None:
    if existing and not overwrite:
        names = ", ".join(str(path) for path in existing[:5])
        if len(existing) > 5:
            names += ", ..."
        raise FileExistsError(f"output already exists; pass --overwrite to replace: {names}")
    for path in existing:
        path.unlink()


def _write_parquet_atomic(frame: pd.DataFrame, path: Path) -> None:
    tmp_path = path.with_name(f"{path.name}.tmp")
    if tmp_path.exists():
        tmp_path.unlink()
    frame.to_parquet(tmp_path, engine="pyarrow", compression="snappy", index=False)
    tmp_path.replace(path)


def _array_to_list(values: Any) -> list[Any]:
    arr = np.asarray(values)
    if arr.ndim == 2 and arr.shape[0] == 1:
        arr = np.asarray(arr[0])
    else:
        arr = np.asarray(arr).reshape(-1)
    return arr.tolist()


def _wavelength_from_header(header: Any) -> np.ndarray | None:
    wmin = _header_float(header, "WMIN")
    binwidth = _header_float(header, "BINWIDTH")
    bincount = _header_int(header, "BINCOUNT")
    if wmin is None or binwidth is None or bincount is None:
        return None
    return wmin + binwidth * np.arange(bincount, dtype=np.float64)


def _header_str(header: Any, key: str) -> str | None:
    try:
        value = header.get(key)
    except Exception:
        return None
    if value is None:
        return None
    return str(value).strip().upper()


def _header_float(header: Any, key: str) -> float | None:
    try:
        value = header.get(key)
    except Exception:
        return None
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


def _header_int(header: Any, key: str) -> int | None:
    try:
        value = header.get(key)
    except Exception:
        return None
    try:
        return None if value is None else int(value)
    except (TypeError, ValueError):
        return None


def _as_text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def _clean_text_or_none(value: Any) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = _as_text(value).strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return None
    return text


def _optional_float(value: Any) -> float | None:
    if value is None or pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _optional_int(value: Any) -> int | None:
    number = _optional_float(value)
    if number is None:
        return None
    try:
        return int(number)
    except (TypeError, ValueError, OverflowError):
        return None


def _failure_record(result: dict[str, Any]) -> dict[str, Any]:
    return {
        "row_index": int(result["row_index"]),
        "file_path": str(result["file_path"]),
        "hdu_index": int(result["hdu_index"]),
        "error": str(result.get("error", "")),
    }


def _write_failures(path: Path, failures: list[dict[str, Any]]) -> None:
    path.write_text(
        "".join(json.dumps(item, sort_keys=True) + "\n" for item in failures),
        encoding="utf-8",
    )


def _write_manifest(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _progress(enabled: bool, done: int, total: int, label: str) -> None:
    if not enabled or total <= 0:
        return
    print(f"\r[{label}] rows {done}/{total} ({(100.0 * done / total):5.1f}%)", end="", file=sys.stderr, flush=True)
