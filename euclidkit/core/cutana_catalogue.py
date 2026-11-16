"""Deterministic Cutana catalogue assembly and alternative-product downloads."""

import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u


def validate_size(unit, value):
    value = float(value)
    limit = 10000 if unit == 'pixel' else 3600
    if not np.isfinite(value) or not 0 < value <= limit:
        raise ValueError(f'diameter_{unit} must be finite, positive and at most {limit}')
    if unit == 'pixel' and not value.is_integer():
        raise ValueError('diameter_pixel must be an integer')
    return int(value) if unit == 'pixel' else value


def _identify(frame, ra, dec):
    frame = frame.copy()
    for name, lower, upper in [(ra, 0, 360), (dec, -90, 90)]:
        frame[name] = pd.to_numeric(frame[name], errors='raise').astype(float)
        if not (np.isfinite(frame[name]) & frame[name].between(lower, upper)).all():
            raise ValueError(f'Invalid {name}: expected finite values in [{lower}, {upper}]')
    # Hash the complete coordinate pair, never rounded display strings.
    position = pd.util.hash_pandas_object(frame[[ra, dec]], index=False).astype(str)
    frame['SourceID'] = 'SRC_' + position
    if 'object_id' in frame:
        valid = frame.object_id.notna()
        frame.loc[valid, 'SourceID'] = 'OBJ_' + frame.loc[valid, 'object_id'].astype(str).str.replace(
            r'\.0$', '', regex=True)
        pairs = frame[['SourceID', ra, dec]].drop_duplicates()
        repeated = pairs.SourceID[pairs.SourceID.duplicated(keep=False)]
        mask = frame.SourceID.isin(repeated)
        frame.loc[mask, 'SourceID'] += '_' + position[mask]
    return frame


def write_catalogue(files, coordinates, output_file, ra, dec, bands, unit, size, keep_columns,
                    environment=None, data_release=None):
    """Choose nearest image centres per source/band; retain all alternate candidates.

    Group catalogues are written for heterogeneous band sets: each can be run with
    its own Cutana channel configuration. The requested output retains all sources
    with at least one matched band and acts as a combined inventory in that case.
    """
    output = Path(output_file)
    keys = [ra, dec] + (['object_id'] if 'object_id' in files and 'object_id' in coordinates else [])
    identities = _identify(pd.concat([coordinates[keys], files[keys]], ignore_index=True)
                           .drop_duplicates(keys), ra, dec)
    files = files.drop(columns=['SourceID'], errors='ignore').merge(
        identities[keys + ['SourceID']], on=keys, validate='many_to_one')
    files['filter_name'] = files.filter_name.astype(str).str.replace('-', '_', regex=False)
    files = files[files.filter_name.isin(bands)].copy()
    if files.empty:
        raise ValueError('No mosaic files matched the requested bands')
    if files[['datalabs_path', 'file_name']].isna().any().any():
        raise ValueError('Mosaic file paths must not be null')
    files['fits_file_paths'] = (files.datalabs_path.astype(str).str.rstrip('/')
                                + '/' + files.file_name.astype(str))
    files['centre_distance_deg'] = np.inf
    if {'image_ra', 'image_dec'} <= set(files.columns):
        ira = pd.to_numeric(files.image_ra, errors='coerce')
        idec = pd.to_numeric(files.image_dec, errors='coerce')
        valid = np.isfinite(ira) & np.isfinite(idec) & idec.between(-90, 90)
        files.loc[valid, 'centre_distance_deg'] = SkyCoord(
            files.loc[valid, ra].to_numpy() * u.deg,
            files.loc[valid, dec].to_numpy() * u.deg,
        ).separation(SkyCoord(ira[valid].to_numpy() * u.deg,
                              idec[valid].to_numpy() * u.deg)).deg
    # Multiple archive records for the same physical file are not alternatives.
    files['_oid'] = files.get('mosaic_product_oid', pd.Series('', index=files.index)).astype(str)
    files = files.sort_values(['SourceID', 'filter_name', 'centre_distance_deg',
                               'fits_file_paths', '_oid'], kind='stable')
    files = files.drop_duplicates(['SourceID', 'filter_name', 'fits_file_paths'])
    files['candidate_rank'] = files.groupby(['SourceID', 'filter_name']).cumcount() + 1
    files['selected'] = files.candidate_rank.eq(1)
    alternatives = files[files.duplicated(['SourceID', 'filter_name'], keep=False)].copy()
    if isinstance(environment, str):
        alternatives['archive_environment'] = environment
    if isinstance(data_release, str):
        alternatives['data_release'] = data_release
    selected = files[files.selected].copy()
    selected['_band_order'] = selected.filter_name.map({b: i for i, b in enumerate(bands)})
    selected = selected.sort_values(['SourceID', '_band_order'])
    identity = ['SourceID', ra, dec] + (['object_id'] if 'object_id' in files else [])
    result = selected.groupby('SourceID', sort=True).agg(
        fits_file_paths=('fits_file_paths', lambda x: json.dumps(list(x))),
        _bands=('filter_name', lambda x: tuple(x)),
    ).reset_index().merge(selected[identity].drop_duplicates(), on='SourceID', validate='one_to_one')
    if keep_columns:
        extra = [c for c in coordinates if c not in result and c not in (ra, dec)]
        if extra:
            # Join by object and position when possible, never just sky position.
            keys = [ra, dec] + (['object_id'] if 'object_id' in coordinates and 'object_id' in result else [])
            metadata = coordinates[keys + extra].drop_duplicates(keys)
            result = result.merge(metadata, on=keys, how='left', validate='many_to_one')
    result = result.rename(columns={ra: 'RA', dec: 'Dec'})
    result.insert(3, 'diameter_' + unit, size)
    groups = result.pop('_bands')
    leading = ['SourceID', 'RA', 'Dec', 'diameter_' + unit]
    if 'object_id' in result:
        leading.append('object_id')
    leading.append('fits_file_paths')
    result = result[leading + [c for c in result if c not in leading]]
    if result.SourceID.duplicated().any():
        raise ValueError('SourceID collision after catalogue assembly')

    output.parent.mkdir(parents=True, exist_ok=True)
    alternatives_path = output.with_name(output.stem + '.alternatives.csv')
    alternatives.drop(columns=['_oid']).to_csv(alternatives_path, index=False)
    present = selected.groupby('SourceID').filter_name.agg(list)
    coverage = identities.copy()
    coverage['available_bands'] = coverage.SourceID.map(present).map(
        lambda x: json.dumps(x) if isinstance(x, list) else '[]')
    coverage['missing_bands'] = coverage.available_bands.map(
        lambda x: json.dumps([b for b in bands if b not in json.loads(x)]))
    coverage['has_mosaic'] = coverage.SourceID.isin(result.SourceID)
    coverage_path = output.with_name(output.stem + '.coverage.csv')
    coverage.to_csv(coverage_path, index=False)
    result.to_csv(output, index=False)
    catalogue_groups = []
    unique_groups = sorted(set(groups))
    for subset in unique_groups:
        mask = groups.map(lambda x: x == subset)
        path = output if len(unique_groups) == 1 else output.with_name(
            output.stem + '.bands-' + '-'.join(subset) + '.csv')
        if path != output:
            result.loc[mask].to_csv(path, index=False)
        catalogue_groups.append({'catalogue': str(path.resolve()), 'bands': list(subset),
                                 'source_count': int(mask.sum()),
                                 'missing_bands': [b for b in bands if b not in subset]})
    report = {'requested_bands': bands, 'groups': catalogue_groups,
              'archive_environment': environment if isinstance(environment, str) else None,
              'data_release': data_release if isinstance(data_release, str) else None,
              'alternatives': str(alternatives_path.resolve()),
              'coverage': str(coverage_path.resolve()),
              'sources_without_mosaics': int((~coverage.has_mosaic).sum()),
              'selection': 'nearest image centre; ties by full path then product OID',
              'sources_with_alternatives': int(alternatives.SourceID.nunique()),
              'combined_catalogue_requires_separate_group_runs': len(catalogue_groups) > 1}
    manifest = output.with_name(output.stem + '.manifest.json')
    manifest.write_text(json.dumps(report, indent=2) + '\n')
    result.attrs['cutana_manifest'] = str(manifest)
    result.attrs['cutana_groups'] = catalogue_groups
    return result


def download_alternatives(archive, candidates_file, output_dir, schema=None):
    """Download distinct unselected files from an alternatives CSV; resume existing files."""
    from euclidkit.core.archive_service import ArchiveService, product_schema_for_client
    schema = product_schema_for_client(archive, schema, getattr(archive, 'data_release', None)
                                       if isinstance(getattr(archive, 'data_release', None), str) else None)
    candidates = pd.read_csv(candidates_file)
    if isinstance(getattr(archive, 'service', None), ArchiveService):
        for column, expected in [('archive_environment', archive.environment),
                                 ('data_release', archive.data_release)]:
            if column in candidates:
                values = set(candidates[column].dropna().astype(str))
                if values and values != {expected}:
                    raise ValueError(f'Alternatives belong to {column}={sorted(values)}; select the matching service/release')
    if not {'selected', 'file_name'} <= set(candidates):
        raise ValueError('Alternatives table requires selected and file_name columns')
    selected = candidates.selected.astype(str).str.lower()
    if not selected.isin(['true', 'false']).all():
        raise ValueError('selected must contain true/false values')
    names = sorted(set(candidates.loc[selected.eq('false'), 'file_name'].dropna().astype(str)))
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    archive._ensure_client()
    paths = []
    for name in names:
        if Path(name).name != name or name in ('', '.', '..'):
            raise ValueError(f'Expected an archive basename, got {name!r}')
        target = directory / name
        if not (target.exists() and target.stat().st_size > 0):
            with tempfile.TemporaryDirectory(dir=directory) as temporary:
                staging = Path(temporary) / name
                downloaded = archive.euclid.get_product(file_name=name, schema=schema,
                                                         output_file=str(staging))
                if not downloaded or not staging.exists() or staging.stat().st_size == 0:
                    raise RuntimeError(f'Failed to download alternative {name}')
                shutil.move(str(staging), str(target))
        paths.append(str(target))
    return paths
