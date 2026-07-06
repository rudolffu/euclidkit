# euclidkit

[![PyPI version](https://img.shields.io/pypi/v/euclidkit.svg)](https://pypi.org/project/euclidkit/)
[![Read the Docs](https://img.shields.io/readthedocs/euclidkit?label=docs)](https://euclidkit.readthedocs.io/en/latest/index.html)

A comprehensive Python package for Euclid archival data analysis, designed for use within the ESA Datalabs environment.

## Overview

`euclidkit` facilitates advanced data exploration and visualization for Euclid Q1/(I)DR1 archival releases, including:

- **Data Access**: Query and crossmatch sources with the Euclid MER catalogue
- **Segmentation Maps**: Resolve MER segmentation-map metadata from crossmatch results
- **Spectroscopic Analysis**: Access, download, and combine NISP spectra of archival sources
- **Unified Workflow**: Streamlined tools for researchers working with Euclid spectroscopic data

The package is designed for efficient archive querying and Euclid spectrum compilation workflows.

## Installation

### Requirements

- Python 3.11+
- Access to ESA Datalabs environment (for data volumes)
- COSMOS credentials for Euclid archive access

### Basic Installation

```bash
pip install euclidkit
```

### Development Installation

```bash
git clone https://github.com/rudolffu/euclidkit.git
cd euclidkit
pip install -e .
```

## Quick Start

### Setup Credentials

Store credentials in a private file under your home directory and restrict permissions:
```bash
mkdir -p ~/.euclidkit
touch ~/.euclidkit/.cred.txt
chmod 600 ~/.euclidkit/.cred.txt
```

Edit `~/.euclidkit/.cred.txt` manually with your preferred editor (do not put credentials in shell history).

Use two lines:
1. COSMOS username
2. COSMOS password

### Configuration

Create and edit the user config file:
```bash
euclidkit init-config --output ~/.euclidkit/euclidkit_config.yaml --template basic
```

Then edit `~/.euclidkit/euclidkit_config.yaml` and set the credential path.

Set the credential path in the config:

```yaml
data:
  credentials_file: /home/<user>/.euclidkit/.cred.txt
```

### Basic Usage

```python
# Note: the Python import path is currently still `euclidkit`.
from euclidkit.core.data_access import EuclidArchive

# Initialize archive connection
archive = EuclidArchive(environment='PDR')
archive.login()

# Crossmatch your sources with Euclid MER catalogue
results = archive.crossmatch_sources(
    user_table="my_sources.csv",
    radius=1.0,  # arcseconds
    output_file="crossmatch_results.fits"
)

# Query for available spectra
spectra_table = archive.query_spectra_sources(
    crossmatch_table=results,
    output_file="spectra_sources.fits"
)

# Export local Datalabs spectra rows to Parquet parts
from euclidkit.core.spectra_parquet import spectra_to_parquet

stats = spectra_to_parquet(
    catalog_table="spectra_sources.fits",
    output_prefix="./output/raw_spectra",
    lambda_range="RGS",
)
print(stats.output_files)
```

## Command Line Interface

### Crossmatching Sources

```bash
# Crossmatch user table with Euclid MER catalogue
euclidkit crossmatch \
    --input my_sources.csv \
    --output crossmatch_results.fits \
    --radius 1.0 \
    --verbose

# Optionally remove columns that are entirely null/missing from the final
# output table. Zero, False, and empty-string columns are retained.
euclidkit crossmatch \
    --input my_sources.csv \
    --output crossmatch_results.fits \
    --drop-empty-columns

# Submit the entire table as a single async job (no batching). The output file
# uses async TAP mode; for very large tables euclidkit splits into async chunks.
euclidkit crossmatch \
    --input my_sources.csv \
    --output crossmatch_results.fits \
    --full-async \
    --async-chunk-size 500000

# When using the IDR environment the command defaults to the WIDE field and
# writes results to wide_<filename>. Use --idr-field DEEP to query DEEP MER.
# DEEP defaults to the deep_survey partition (EDFN, EDFF, EDFS).
euclidkit crossmatch \
    --input my_sources.csv \
    --output crossmatch_results.fits \
    --environment IDR \
    --idr-field DEEP

# Query the deep_mode partition (CDFS, COSMOS), or use "both" for
# deep_survey followed by deep_mode.
euclidkit crossmatch \
    --input my_sources.csv \
    --output crossmatch_results.fits \
    --environment IDR \
    --idr-field DEEP \
    --idr-deep-partition mode

# Crossmatch an already-uploaded archive user table (no local upload needed)
euclidkit crossmatch \
    --user-table-name my_table \
    --output crossmatch_results.fits \
    --match-mode object-id \
    --environment IDR \
    --idr-field WIDE
```

`--full-async` behavior:

- For smaller inputs, euclidkit submits one async TAP job, downloads the result to the requested output file, and then removes the remote job.
- For large local input tables (`--input`), euclidkit splits the upload into async chunks, saves each chunk to `<output>_part_####.fits`, removes each remote job after the chunk is saved, writes `<output>.manifest.json`, and merges the chunk files into the requested final output.
- For large archive user tables (`--user-table-name`), euclidkit uses the same on-disk chunking pattern and final merge.

Matching mode recommendation:

- Prefer `--match-mode object-id` whenever the input already contains Euclid `object_id` values, or `source_id` values that should be joined to MER `object_id`. This avoids positional matching and is usually faster and more robust for large tables.

`--max-sources` vs `--async-chunk-size`:

- `--max-sources`: limits how many rows from the input table are processed in total.
- `--async-chunk-size`: controls rows per async TAP job when `--full-async` is enabled.
- `--drop-empty-columns`: drops columns where every result value is null or missing before saving the final `--output` table. Intermediate async part files are left unchanged.

### Uploading Tables

```bash
# Upload a FITS table to your Euclid TAP workspace
euclidkit upload-table \
    --input my_sources.fits \
    --table-name my_workspace_table \
    --description "Sources awaiting deep crossmatch" \
    --overwrite

# Upload CSV data as-is (format inferred automatically)
euclidkit upload-table \
    --input trimmed_sources.csv \
    --table-name trimmed_sources
```

### Querying Spectra

```bash
# Query spectra-source rows from an ID or coordinate table.
# Auto mode uses object_id when present, otherwise spatial RA/Dec matching.
euclidkit query-spectra \
    --crossmatch my_spectral_ids_or_coordinates.fits \
    --output spectra_sources.fits \
    --environment IDR \
    --idr-field WIDE \
    --verbose
```

`query-spectra` supports `--match-mode auto|object-id|spatial`. Object-ID mode
joins `spectra_source.source_id = object_id`; spatial mode matches input
coordinates to `ra_obj`/`dec_obj` with nearest-only output. Spatial mode accepts
common RA/Dec aliases such as `RA`/`DEC`, `right_ascension`/`declination`, and
`ra_deg`/`dec_deg` case-insensitively. A MER crossmatch table is optional; any
local table with the needed spectra-source IDs or coordinates can be used. The
result is the usual input to `compile-spectra`; see the Sphinx spectra
compilation guide for the full Datalabs and Datalink workflow details.

### Querying Segmentation Maps

```bash
# Resolve segmentation-map files and WCS metadata from MER crossmatch results.
# The input must contain SEGMENTATION_MAP_ID, object_id, and coordinates.
# query-segmap prefers ra/dec and falls back to mer_ra/mer_dec.
euclidkit query-segmap \
    --input crossmatch_results.fits \
    --output segmentation_maps.fits \
    --environment IDR

# Create one raw-label 10 arcsec FITS cutout per source row
euclidkit compile-segmap \
    --input segmentation_maps.fits \
    --output-dir ./segmap_cutouts
```

`query-segmap` computes `tile_index = floor(SEGMENTATION_MAP_ID / 1_000_000)`
locally, then joins to `q1.mer_segmentation_map` for PDR,
`dr1.mer_segmentation_map` for IDR, and `sedm.mer_segmentation_map` for OTF/REG.
Run `euclidkit crossmatch` first if your table does not yet contain
`SEGMENTATION_MAP_ID`. `compile-segmap` reads local segmentation-map FITS files
from `datalabs_path` + `file_name`, groups rows by tile file, and preserves the
raw segmentation-label pixels in each cutout. See the segmentation-map guide in
the Sphinx docs for required columns, output filenames, and error-handling
options.

### Building Cutana Input

```bash
# Build Cutana CSV from a source table with object_id or ra/dec columns
euclidkit query-cutana \
    --sources my_sources.fits \
    --output cutana_input.csv \
    --instrument VIS \
    --cutout-size arcsec \
    --cutout-size-value 15

# NISP example with explicit filters
euclidkit query-cutana \
    --sources my_sources.fits \
    --output cutana_input_nisp.csv \
    --instrument NISP \
    --nisp-filters NIR_Y,NIR_H \
    --environment IDR \
    --idr-field DEEP \
    --idr-deep-partition both \
    --cutout-size arcsec \
    --cutout-size-value 15
```

### Compiling Spectra

```bash
# Default local Datalabs mode: export raw spectra to parquet parts
euclidkit compile-spectra \
    --spectra-table spectra_sources.fits \
    --output-dir ./output \
    --prefix raw_spectra \
    --chunk-size 2000 \
    --workers 8 \
    -L RGS

# Export both arms into separate raw_spectra_rgs / raw_spectra_bgs parquet families
euclidkit compile-spectra \
    --spectra-table spectra_sources.fits \
    --output-dir ./output \
    --prefix raw_spectra \
    -L BOTH

# Legacy local FITS mode remains available explicitly
euclidkit compile-spectra \
    --spectra-table spectra_sources.fits \
    --output-dir ./output \
    --prefix compiled_spectra \
    --output-format fits \
    --max-extensions 1000

# Datalink mode is unchanged and writes FITS outputs
euclidkit compile-spectra \
    --spectra-table spectra_sources.fits \
    --output-dir ./output \
    --prefix compiled_dl \
    --use-datalink \
    --environment IDR \
    --schema sedm \
    -L BOTH

# Export per-dither spectra from local Datalabs FITS files
euclidkit dithers-to-parquet \
    --catalog-table spectra_sources.fits \
    --output-prefix ./output/raw_sir \
    --workers 8 \
    --lambda-range RGS \
    --environment IDR
```

Note: non-Datalink `compile-spectra` now defaults to Parquet and reads `LRANGE` from FITS headers for `-L/--lambda-range` filtering. Datalink remains FITS-only; `RGS`/`BGS` map to corresponding retrieval types, and `BOTH` runs two passes and writes separate `_rgs` and `_bgs` FITS files. `--retrieval-type` is kept for backward compatibility.

`dithers-to-parquet` automatically annotates per-dither rows with
`obs_time_mjd`, `obs_time_utc`, and `pa` from the archive raw-frame table
filtered to `technique = 'SPECTROIMAGE'`. PDR/Q1 uses `q1.raw_frame`; IDR/DR1
uses `dr1.raw_frame`; OTF/REG use `sedm.raw_frame`. It matches
`raw_frame.pointing_id` to the dither HDU `ptgid`, falling back to parsed
`dither_id` only when `ptgid` is unavailable. Internally, it uses a TAP upload
join rather than a long literal ID-list clause:

```sql
SELECT
  r.pointing_id,
  r.obs_time_mjd,
  r.obs_time_utc,
  r.pa
FROM dr1.raw_frame AS r
JOIN TAP_UPLOAD.dither_pointings AS p
  ON r.pointing_id = p.pointing_id
WHERE r.technique = 'SPECTROIMAGE'
```

## Key Features

### Data Archive Integration

- **Multiple Environments**: Support for PDR, IDR, OTF, and REG archive environments
- **Efficient Queries**: Batch processing with TAP table uploads for large datasets
- **Crossmatching**: Position-based matching with configurable search radius

### Spectroscopic Tools

- **Spectrum Access**: Direct access to Euclid data volumes on ESA Datalabs
- **Parquet Spectra Export**: Export local Datalabs spectra to raw Parquet parts by default
- **FITS Compatibility**: Keep legacy multi-extension FITS compilation and Datalink FITS outputs available

### Analysis Pipeline

- **Quality Control**: Spectrum validation and quality assessment

## Data Environment

### ESA Datalabs Integration

This package is optimized for the ESA Datalabs environment with direct access to:

- **Euclid Q1 Data**: `/data/euclid_q1/` (35 TB volume)

## API Reference

### Core Classes

#### `EuclidArchive`

Main interface to the Euclid science archive.

```python
archive = EuclidArchive(environment='PDR')
archive.login(credentials_file='~/.euclidkit/.cred.txt')

# Crossmatch sources
results = archive.crossmatch_sources(
    user_table="sources.csv",
    radius=1.0,
    output_file="results.fits"
)

# Query spectra
spectra = archive.query_spectra_sources(
    crossmatch_table=results,
    output_file="spectra.fits"
)

# Get individual spectrum
spectrum_hdu = archive.get_individual_spectrum(
    datalabs_path="/data/euclid_q1/path",
    file_name="spectrum_file.fits", 
    hdu_index=42
)

# Export queried local Datalabs spectra to raw Parquet parts
from euclidkit.core.spectra_parquet import spectra_to_parquet

stats = spectra_to_parquet(
    catalog_table="spectra.fits",
    output_prefix="./output/raw_spectra",
    chunk_size=2000,
    lambda_range="RGS",
)
print(stats.manifest_path)
```

#### Spectra Parquet Export

Default local Datalabs spectra export writes raw Parquet parts directly from
catalog rows containing `datalabs_path`, `file_name`, `hdu_index`, `source_id`,
`object_id`, `ra_obj`, and `dec_obj`.

```python
from euclidkit.core.spectra_parquet import spectra_to_parquet, dithers_to_parquet

raw_stats = spectra_to_parquet(
    catalog_table="spectra_sources.fits",
    output_prefix="./output/raw_spectra",
    chunk_size=2000,
    workers=8,
    lambda_range="RGS",
    on_error="skip",
)

dither_stats = dithers_to_parquet(
    catalog_table="spectra_sources.fits",
    output_prefix="./output/raw_sir",
    chunk_size=2000,
    workers=8,
    lambda_range="RGS",
    include_combined=True,
    environment="IDR",
)
```

`SpectrumCompiler` remains available for legacy local multi-extension FITS
compilation and Datalink FITS outputs.

### Workflow Examples

#### Complete Spectroscopic Analysis Pipeline

```python
from euclidkit.core.data_access import EuclidArchive
from euclidkit.core.spectra_parquet import spectra_to_parquet, dithers_to_parquet
import pandas as pd

# 1. Initialize archive
archive = EuclidArchive(environment='PDR')
archive.login()

# 2. Load your QSO candidates
qso_candidates = pd.read_csv('qso_candidates.csv')

# 3. Crossmatch with Euclid MER catalogue
crossmatches = archive.crossmatch_sources(
    user_table=qso_candidates,
    radius=2.0,  # 2 arcsecond radius
    output_file='qso_crossmatches.fits'
)

# 4. Find available spectra
spectra_sources = archive.query_spectra_sources(
    crossmatch_table=crossmatches,
    output_file='qso_spectra_sources.fits'
)

print(f"Found {len(spectra_sources)} spectra for {len(crossmatches)} crossmatches")

# 5. Export raw spectra to Parquet parts for downstream analysis
raw_stats = spectra_to_parquet(
    catalog_table='qso_spectra_sources.fits',
    output_prefix='./spectra_parquet/qso_raw',
    chunk_size=2000,
    lambda_range='RGS',
    on_error='skip',
)
print(f"Wrote {len(raw_stats.output_files)} raw spectra parquet parts")

# 6. Optionally export per-dither spectra for the same catalog rows
dither_stats = dithers_to_parquet(
    catalog_table='qso_spectra_sources.fits',
    output_prefix='./spectra_parquet/qso_sir',
    lambda_range='RGS',
    environment='IDR',
)
print(f"Wrote {len(dither_stats.dither_output_files)} dither parquet parts")

archive.logout()
```

## Diagnostics

Check your installation and environment:

```bash
# Check all components
euclidkit diagnostics

# Check specific components
euclidkit diagnostics --check-deps --check-data
```

## Archive Environments

Use ``--environment`` (CLI) or ``environment=...`` (Python API) to select the
archive backend:

- **PDR**: Public Data Release archive.
- **IDR**: Internal Data Release archive (consortium access).
- **OTF**: On-the-fly archive environment.
- **REG**: Regression/testing archive environment.

For **IDR**, you can also select the field with ``--idr-field``:

- **WIDE**: Uses the IDR WIDE MER catalogue. This queries
  ``catalogue.mer_catalogue_wide_survey`` first, then queries
  ``catalogue.mer_catalogue_wide_mode`` only for sources not matched in the
  survey table.
- **DEEP**: Uses IDR DEEP MER partitions. By default this queries
  ``catalogue.mer_catalogue_deep_survey`` (EDFN, EDFF, EDFS). Use
  ``--idr-deep-partition mode`` for ``catalogue.mer_catalogue_deep_mode``
  (CDFS, COSMOS), or ``--idr-deep-partition both`` to query survey first and
  mode second.

Examples:

```bash
# IDR WIDE (default IDR field)
euclidkit crossmatch \
  --input my_sources.fits \
  --output xmatch_wide.fits \
  --environment IDR \
  --idr-field WIDE

# IDR DEEP
euclidkit crossmatch \
  --input my_sources.fits \
  --output xmatch_deep.fits \
  --environment IDR \
  --idr-field DEEP

# IDR DEEP mode partition
euclidkit crossmatch \
  --input my_sources.fits \
  --output xmatch_deep_mode.fits \
  --environment IDR \
  --idr-field DEEP \
  --idr-deep-partition mode
```

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Documentation

For detailed documentation and examples, visit:
- [euclidkit docs](https://euclidkit.readthedocs.io/en/latest/index.html)
- [Package Documentation](https://github.com/rudolffu/euclidkit/docs)
- [Euclid Science Archive](https://s2e2.cosmos.esa.int/www/euclid_iscience/Public_User_Guide.html)
- [astroquery.esa.euclid](https://astroquery.readthedocs.io/en/latest/esa/euclid/euclid.html)

## Support

- **Issues**: [GitHub Issues](https://github.com/rudolffu/euclidkit/issues)
- **Discussions**: [GitHub Discussions](https://github.com/rudolffu/euclidkit/discussions)
- **Email**: fuympku@outlook.com

## Author

**Yuming Fu** ([@rudolffu](https://github.com/rudolffu))
- Email: fuympku@outlook.com
- GitHub: https://github.com/rudolffu/euclidkit

## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- ESA Euclid Mission and Euclid Consortium
- ESA Datalabs and Euclid Data Space infrastructure team
- Astropy and astroquery communities

## Changelog

### Latest Changes

- **Spectroscopic Pipeline**: Complete pipeline for accessing and combining Euclid spectra
- **CLI Integration**: Added `--combine-output` option to `query-spectra` command
- **TAP Upload**: Improved query performance using TAP table uploads
- **Parquet Spectra Export**: Default local Datalabs spectra export to Parquet
- **Error Handling**: Robust handling of long filenames and missing data

See [CHANGELOG.md](CHANGELOG.md) for detailed version history.
