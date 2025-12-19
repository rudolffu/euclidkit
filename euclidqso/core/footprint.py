"""
Utilities for selecting catalog entries within Euclid survey footprints defined
by Multi-Order Coverage (MOC) maps.
"""

import logging
from importlib import resources
from pathlib import Path
from typing import Dict, Optional

import astropy.units as u
from astropy.coordinates import SkyCoord
from mocpy import MOC

from euclidqso.utils.io import load_table, save_table

logger = logging.getLogger(__name__)

# Map available MOCs by data release and survey
MOC_FILES: Dict[str, Dict[str, str]] = {
    "DR1": {
        "WIDE": "dr1_mer_wide_bins_o13_moc.fits",
        "DEEP": "dr1_mer_deep_o13_moc.fits",
        "BOTH": "dr1_mer_wide_deep_union_o13_moc.fits",
        "CGV_INPUT": "cgv_map_dr1input_o13_moc.fits",
    }
}


# Common aliases (normalized) for MOC selection keys
MOC_ALIASES: Dict[str, str] = {
    "ALL": "BOTH",
    "UNION": "BOTH",
    "WIDE_DEEP": "BOTH",
    "WIDE+DEEP": "BOTH",
    "CGV": "CGV_INPUT",
    "DR1INPUT": "CGV_INPUT",
}


def _normalize_key(value: str) -> str:
    return value.strip().upper().replace("-", "_").replace(" ", "_")


def list_available_surveys(data_release: str = "DR1") -> tuple[str, ...]:
    """Return available packaged footprint keys for a given data release."""
    release_key = _normalize_key(data_release)
    if release_key not in MOC_FILES:
        return tuple()
    return tuple(sorted(MOC_FILES[release_key].keys()))


def _resolve_column(table, target: str) -> str:
    """Return the actual column name matching target (case-insensitive)."""
    for col in table.colnames:
        if col.lower() == target.lower():
            return col
    raise KeyError(f"Column '{target}' not found in table (available: {table.colnames})")


def get_moc_path(survey: str, data_release: str = "DR1", moc_path: Optional[str] = None) -> Path:
    """
    Resolve the path to the requested MOC file.
    """
    if moc_path:
        return Path(moc_path)

    release_key = _normalize_key(data_release)
    survey_key = _normalize_key(survey)
    survey_key = MOC_ALIASES.get(survey_key, survey_key)

    if release_key not in MOC_FILES or survey_key not in MOC_FILES[release_key]:
        available = ", ".join(
            f"{rel}:{'/'.join(fields.keys())}" for rel, fields in MOC_FILES.items()
        )
        raise ValueError(
            f"No packaged MOC for data release '{data_release}' and survey '{survey}'. "
            f"Available combinations: {available}"
        )

    moc_filename = MOC_FILES[release_key][survey_key]
    moc_resource = resources.files("euclidqso").joinpath("data", "mocs", moc_filename)

    try:
        with resources.as_file(moc_resource) as resolved_path:
            return Path(resolved_path)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"MOC file '{moc_filename}' not found in packaged data."
        ) from exc


def filter_catalog_by_moc(
    input_catalog: str,
    output_catalog: str,
    survey: str = "WIDE",
    data_release: str = "DR1",
    ra_col: str = "ra",
    dec_col: str = "dec",
    moc_path: Optional[str] = None,
):
    """
    Filter catalog rows to those lying within the survey MOC footprint.

    Parameters
    ----------
    input_catalog : str
        Path to input catalog (FITS, CSV, VOTable, etc.).
    output_catalog : str
        Path to write the filtered catalog.
    survey : str, default "WIDE"
        Packaged footprint key to use (e.g., "WIDE", "DEEP", "BOTH", "CGV_INPUT").
    data_release : str, default "DR1"
        Data release version. Currently only DR1 MOCs are packaged.
    ra_col : str, default "ra"
        Column name for right ascension in degrees.
    dec_col : str, default "dec"
        Column name for declination in degrees.
    moc_path : str, optional
        Explicit MOC file path. Overrides packaged defaults.

    Returns
    -------
    astropy.table.Table
        Filtered subset of the catalog.
    """
    catalog = load_table(input_catalog)
    ra_name = _resolve_column(catalog, ra_col)
    dec_name = _resolve_column(catalog, dec_col)

    moc_file = get_moc_path(survey=survey, data_release=data_release, moc_path=moc_path)
    moc = MOC.from_fits(moc_file)

    coords = SkyCoord(
        ra=catalog[ra_name],
        dec=catalog[dec_name],
        unit=(u.deg, u.deg),
        frame="icrs",
    )
    # mocpy expects lon/lat inputs, not a SkyCoord instance
    mask = moc.contains(coords.ra, coords.dec)
    filtered = catalog[mask]

    logger.info(
        "Selected %d of %d sources inside %s %s footprint",
        len(filtered),
        len(catalog),
        survey.upper(),
        data_release.upper(),
    )

    save_table(filtered, output_catalog)
    return filtered
