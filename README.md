# euclidkit

[![PyPI version](https://img.shields.io/pypi/v/euclidkit.svg)](https://pypi.org/project/euclidkit/)
[![Read the Docs](https://img.shields.io/readthedocs/euclidkit?label=docs)](https://euclidkit.readthedocs.io/en/latest/index.html)

`euclidkit` takes you from a list of targets to Euclid data for those targets:
MER photometry, morphology, spectroscopic redshifts, NISP slitless spectra,
image cutouts, and segmentation maps.

You bring a table of sky positions or Euclid object IDs. `euclidkit` writes the
archive queries, handles uploads and large asynchronous jobs, and saves the
results as FITS, CSV, or Parquet tables that you can open with Astropy or pandas.

It works with the public Euclid releases (DR1 by default, and Q1) and with the
internal consortium archive. Catalogue queries run from any computer. Bulk
access to spectra and images reads the Euclid data volumes on
[ESA Datalabs](https://datalabs.esa.int/).

## Example

Find Euclid counterparts for your sources, then collect their NISP spectra:

```bash
# 1. Match your positions to the Euclid MER catalogue (1 arcsec radius)
euclidkit crossmatch -i my_sources.csv -o matches.fits --radius 1.0

# 2. Find which of the matched objects have spectra
euclidkit query-spectra -x matches.fits -o spectra_sources.fits

# 3. On ESA Datalabs: export the red-grism spectra to Parquet
euclidkit compile-spectra -s spectra_sources.fits -o ./spectra -L RGS
```

`my_sources.csv` needs `ra` and `dec` columns in degrees (ICRS). The
[quickstart](https://euclidkit.readthedocs.io/en/latest/quickstart.html) walks
through each step and its output.

## What you can do

| Goal | Command | Runs outside Datalabs? |
|---|---|---|
| Keep only sources inside the Euclid footprint | `select-footprint` | Yes |
| Get MER photometry and flags for your sources | `crossmatch` | Yes |
| Get morphology (CAS, Gini, M20, Sérsic) | `query-morph` | Yes |
| Get spectroscopic redshift candidates | `query-zspe` | Yes |
| Find and export NISP spectra | `query-spectra`, `compile-spectra` | Query: yes. Export: Datalabs, or Datalink from anywhere |
| Export individual dither spectra | `dithers-to-parquet` | No |
| Prepare image cutouts with [Cutana](https://github.com/esa/Cutana) | `query-cutana` | Query: yes |
| Cut out segmentation maps | `query-segmap`, `compile-segmap` | Query: yes. Cutouts: no |

Every command is also available from Python through
`euclidkit.core.data_access.EuclidArchive`.

## Installation

```bash
pip install euclidkit
```

Requires Python 3.11 or later. See the
[installation guide](https://euclidkit.readthedocs.io/en/latest/installation.html)
for a development install.

## Archive credentials

Archive queries need your ESA Cosmos account. Store the username and password
in a private two-line file, which `euclidkit` finds automatically:

```bash
mkdir -p ~/.euclidkit
touch ~/.euclidkit/.cred.txt
chmod 600 ~/.euclidkit/.cred.txt
```

Open `~/.euclidkit/.cred.txt` in an editor. Put your username on the first line
and your password on the second. Do not type the password on the command line,
because the shell records it in its history. Other options are in the
[configuration guide](https://euclidkit.readthedocs.io/en/latest/configuration.html).

## Documentation

Full documentation: <https://euclidkit.readthedocs.io>

- [Quickstart](https://euclidkit.readthedocs.io/en/latest/quickstart.html): a complete worked example
- [Data releases and fields](https://euclidkit.readthedocs.io/en/latest/environments.html): DR1, Q1, WIDE, and DEEP
- [Catalogues](https://euclidkit.readthedocs.io/en/latest/catalogues.html), [spectra](https://euclidkit.readthedocs.io/en/latest/spectra_compilation.html), [image cutouts](https://euclidkit.readthedocs.io/en/latest/cutouts.html), and [segmentation maps](https://euclidkit.readthedocs.io/en/latest/segmentation_maps.html)
- [Command reference](https://euclidkit.readthedocs.io/en/latest/cli.html) and [Python API](https://euclidkit.readthedocs.io/en/latest/api.html)
- [Glossary](https://euclidkit.readthedocs.io/en/latest/glossary.html) of Euclid terms
- [Changelog](CHANGELOG.md)

Related resources:

- [Euclid Science Archive](https://eas.esac.esa.int/sas/)
- [Euclid DR1 Data Product Description](https://euclid.esac.esa.int/dr/dr1/dpdd/), for the definition and units of each column
- [astroquery.esa.euclid](https://astroquery.readthedocs.io/en/latest/esa/euclid/euclid.html), which `euclidkit` builds on

## Support

Report bugs and ask questions on the
[GitHub issue tracker](https://github.com/rudolffu/euclidkit/issues). You can
also email the author, Yuming Fu ([@rudolffu](https://github.com/rudolffu)), at
fuympku@outlook.com.

## Contributing

Contributions are welcome. Open an issue first if the change is large.

```bash
git clone https://github.com/rudolffu/euclidkit.git
cd euclidkit
pip install -e ".[dev]"
pytest -q
```

The tests mock the archive and the data volumes, so they run without
credentials. Add tests for new behaviour, and send a pull request against
`main`.

## License

BSD 3-Clause. See [LICENSE](LICENSE).

## Acknowledgments

- The ESA Euclid mission and the Euclid Consortium
- The ESA Datalabs and Euclid archive teams
- Kristin Anett Remmelgas and Héctor Cánovas Cabrera, for the spectra-source query approach
- The Astropy and astroquery communities
