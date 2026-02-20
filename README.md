# euclidkit

A comprehensive Python package for Euclid archival data analysis, designed for use within the ESA Datalabs environment.

## Overview

`euclidkit` facilitates advanced data exploration and visualization for Euclid Q1/(I)DR1 archival releases, including:

- **Data Access**: Query and crossmatch sources with the Euclid MER catalogue
- **Spectroscopic Analysis**: Access, download, and combine NISP spectra of archival sources
- **Template Generation**: Run redshift determination workflows and generate composite spectra
- **Multi-Survey Integration**: Interface with DESI, SDSS, and other survey data
- **Unified Workflow**: Streamlined tools for researchers working with Euclid spectroscopic data

The package integrates with other specialized tools (`specbox` and `qsofitmore`) to provide a complete analysis pipeline.

## Installation

### Requirements

- Python 3.9+
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

### Optional Dependencies

For DESI integration:
```bash
pip install euclidkit[desi]
```

## Quick Start

### Setup Credentials

Create a credentials file in your workspace:
```bash
# Create /media/user/cred.txt with your COSMOS credentials
echo "your_username" > /media/user/cred.txt
echo "your_password" >> /media/user/cred.txt
```

### Configuration

Generate a configuration file:
```bash
euclidkit init-config --output my_config.yaml --template basic
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

# Combine spectra into a single FITS file
combined_file = archive.combine_spectra_to_fits(
    spectra_table=spectra_table,
    output_file="my_combined_spectra.fits"
)
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

# Submit the entire table as a single async job (no batching). The output file
# will contain TAP job metadata instead of immediate crossmatch results.
euclidkit crossmatch \
    --input my_sources.csv \
    --output crossmatch_results.fits \
    --full-async

# When using the IDR environment the command defaults to the WIDE field and
# writes results to wide_<filename>. Use --idr-field DEEP to query the deep stack:
euclidkit crossmatch \
    --input my_sources.csv \
    --output crossmatch_results.fits \
    --environment IDR \
    --idr-field DEEP
```

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
# Query spectra from crossmatch results
euclidkit query-spectra \
    --crossmatch crossmatch_results.fits \
    --output spectra_sources.fits \
    --verbose

# Query spectra by object IDs and auto-combine
euclidkit query-spectra \
    --object-ids 123456,789012,345678 \
    --output spectra_sources.fits \
    --combine-output my_spectra.fits \
    --max-spectra 100 \
    --verbose
```

### Building Cutana Input

```bash
# Build Cutana CSV from a source table with object_id or ra/dec columns
euclidkit query-cutana \
    --sources my_sources.fits \
    --output cutana_input.csv \
    --instrument VIS \
    --cutout-size pixel \
    --cutout-size-value 50

# NISP example with explicit filters
euclidkit query-cutana \
    --sources my_sources.fits \
    --output cutana_input_nisp.csv \
    --instrument NISP \
    --nisp-filters NIR_Y,NIR_H \
    --cutout-size arcsec \
    --cutout-size-value 15
```

### Compiling Spectra

```bash
# Compile individual spectra into chunked FITS files
euclidkit compile-spectra \
    --spectra-table spectra_sources.fits \
    --output-dir ./output \
    --prefix compiled_spectra \
    --max-extensions 5000 \
    --verbose
```

## Key Features

### Data Archive Integration

- **Multiple Environments**: Support for PDR, IDR, OTF, and REG archive environments
- **Efficient Queries**: Batch processing with TAP table uploads for large datasets
- **Crossmatching**: Position-based matching with configurable search radius

### Spectroscopic Tools

- **Spectrum Access**: Direct access to Euclid data volumes on ESA Datalabs
- **FITS Compilation**: Combine individual spectra into multi-extension FITS files
- **Metadata Preservation**: Maintain source IDs, coordinates, and provenance information

### Analysis Pipeline

- **Template Generation**: Tools for creating composite spectra and templates
- **Redshift Analysis**: Cross-correlation redshift determination
- **Quality Control**: Spectrum validation and quality assessment

## Data Environment

### ESA Datalabs Integration

This package is optimized for the ESA Datalabs environment with direct access to:

- **Euclid Q1 Data**: `/data/euclid_q1/` (35 TB volume)
- **Euclid ERO Data**: `/data/euc_ero_data_01/` (600 GB volume)

### File Naming Conventions

The package handles standard Euclid data products:

- **VIS Images**: `EUC_VIS_*_TILE*.fits`
- **NIR Images**: `EUC_NIR_*_{Y,J,H}_*.fits`
- **MER Mosaics**: `EUC_MER_BGSUB-MOSAIC-{VIS,NIR-Y,NIR-J,NIR-H}_TILE*.fits`
- **Spectra**: `EUC_SIR_W-COMBSPEC_*.fits`

## API Reference

### Core Classes

#### `EuclidArchive`

Main interface to the Euclid science archive.

```python
archive = EuclidArchive(environment='PDR')
archive.login(credentials_file='/media/user/cred.txt')

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

# Combine spectra
combined = archive.combine_spectra_to_fits(
    spectra_table=spectra,
    output_file="combined.fits",
    max_spectra=1000
)
```

#### `SpectrumCompiler`

Advanced spectrum compilation with chunking support.

```python
from euclidkit.core.spectra import SpectrumCompiler

compiler = SpectrumCompiler(max_extensions=5000)

# Compile into chunked files
output_files = compiler.compile_spectra(
    spectra_table=spectra_table,
    output_dir="./output",
    output_prefix="compiled_spectra"
)

# Create single FITS file
single_file = compiler.compile_single_fits(
    spectra_table=spectra_table,
    output_file="all_spectra.fits"
)

# Generate metadata table
metadata = compiler.create_metadata_table(
    spectra_table=spectra_table,
    output_files=output_files,
    output_dir="./output"
)
```

### Workflow Examples

#### Complete Spectroscopic Analysis Pipeline

```python
from euclidkit.core.data_access import EuclidArchive
from euclidkit.core.spectra import SpectrumCompiler
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

# 5. Create combined FITS file (for small samples)
if len(spectra_sources) <= 1000:
    combined_spectra = archive.combine_spectra_to_fits(
        spectra_table=spectra_sources,
        output_file='qso_combined_spectra.fits'
    )
    print(f"Combined spectra saved to: {combined_spectra}")

# 6. Or use chunked compilation for large samples
else:
    compiler = SpectrumCompiler(max_extensions=2000)
    output_files = compiler.compile_spectra(
        spectra_table=spectra_sources,
        output_dir='./spectra_chunks',
        output_prefix='qso_spectra'
    )
    print(f"Created {len(output_files)} chunked files")

archive.logout()
```

## Diagnostics

Check your installation and environment:

```bash
# Check all components
euclidkit diagnostics

# Check specific components
euclidkit diagnostics --check-deps --check-data --check-desi
```

## Integration with Other Tools

### DESI Spectra Access

```python
# Access DESI spectra in SDSS-compatible format
from euclidkit.external.desi import get_desi_spectra

desi_spectra = get_desi_spectra(
    ra=150.0, dec=12.5, 
    dr='DESI-DR1',
    output_format='sdss'
)
```

### Template Generation

```python
from euclidkit.analysis.templates import TemplateGenerator

generator = TemplateGenerator()
template = generator.create_composite_template(
    spectra_list=['spec1.fits', 'spec2.fits', 'spec3.fits'],
    redshift_range=(2.0, 4.0),
    output_file='qso_template.fits'
)
```

## Archive Environments

- **PDR**: Public Data Release (`catalogue.mer_catalogue`)
- **IDR**: Internal Data Release (`catalogue.mer_catalogue`)  
- **OTF**: On-the-fly processing (`catalogue.mer_catalogue`)
- **REG**: Regression testing (`catalogue.mer_final_catalog_fits_file_regreproc1_r2`)

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Documentation

For detailed documentation and examples, visit:
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

This project is licensed under the GNU General Public License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- ESA Euclid Mission and Consortium
- ESA Datalabs infrastructure team
- Astropy and astroquery communities
- DESI Collaboration for SPARCL integration

## Changelog

### Latest Changes

- **Spectroscopic Pipeline**: Complete pipeline for accessing and combining Euclid spectra
- **CLI Integration**: Added `--combine-output` option to `query-spectra` command
- **TAP Upload**: Improved query performance using TAP table uploads
- **FITS Compilation**: Efficient multi-extension FITS file creation
- **Error Handling**: Robust handling of long filenames and missing data

See [CHANGELOG.md](CHANGELOG.md) for detailed version history.
