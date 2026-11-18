# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0] - 2026-11-18

### Migration
- Python 3.11 or later is required.
- PDR now defaults to public DR1. Add `--release Q1` or
  `EuclidArchive(data_release="Q1")` to retain the previous public dataset.
- IDR remains an independent internal service. Omit public-release flags when
  selecting IDR; it is not redirected to PDR.
- Use `--field` and `--deep-partition` for DR1/IDR field selection. The legacy
  `--idr-field` and `--idr-deep-partition` aliases and Python keywords still work.

### Added
- Public DR1 routing for MER, spectra sources, SPE candidates, morphology,
  segmentation maps, mosaics and raw-frame queries, with explicit public Q1 selection.
- `query-morph` and `EuclidArchive.query_morphology` for compact MER morphology
  measurements, including WIDE fallback, DEEP partitions and asynchronous queries.
- Deterministic per-band Cutana mosaic selection, alternative-product and coverage
  reports, and `download-cutana-alternatives` for later retrieval.
- Band-group Cutana catalogues that retain sources with partial band coverage.
  Process each group separately with its available channels.
- Service and release provenance in export manifests and async recovery records.

### Changed
- Cutana preparation validates sizes and coordinates, guarantees unique source IDs,
  and checks source centres against archive mosaic footprints.
- Documentation rewritten for astronomers with task-based guides, a command
  reference and glossary, and the Furo theme for readable wide tables.
- Removed tests for retired APIs and refreshed current archive, spectra and I/O
  tests with isolated credentials and temporary files.

### Fixed
- Local FITS compilation copies reused HDUs before assigning row-specific headers,
  preserving identifiers and coordinates for repeated HDU references.
- Compiled FITS `NSPECTRA` counts written extensions after missing inputs are skipped.
- Removed `euclidkit-spectra`, `euclidkit-composite` and `euclidkit-pipeline`
  entry points, which referenced modules that were not shipped.
- Diagnostics use the configured credentials location and report available data
  volumes. Standalone `--check-data` and `--check-desi` checks no longer crash.
- Package license metadata now matches BSD-3-Clause, and `requirements.txt`
  includes `pyarrow`.

## Earlier releases (0.2.1–0.2.4)

The following consolidated notes were previously listed under Unreleased;
they describe changes already present in v0.2.4.

### Added
- IDR DEEP MER partition selection via `idr_deep_partition` and
  `--idr-deep-partition` for MER crossmatch and Cutana workflows.
- Optional empty-column pruning for crossmatch outputs via
  `drop_empty_columns` and `euclidkit crossmatch --drop-empty-columns`.
- Local Datalabs spectra Parquet export via `spectra_to_parquet`, the default
  non-Datalink `euclidkit compile-spectra` path, and new
  `euclidkit dithers-to-parquet` command.
- `euclidkit query-segmap` and `EuclidArchive.query_segmentation_maps` for
  resolving MER segmentation-map metadata from crossmatch outputs.
- `euclidkit compile-segmap` and `compile_segmap_cutouts` for local raw-label
  FITS cutouts from MER segmentation-map tiles.

### Changed
- IDR DEEP MER queries now use `catalogue.mer_catalogue_deep_survey` by
  default, with optional `deep_mode` or `both` partition selection.
- Non-Datalink `compile-spectra` now defaults to Parquet output with automatic
  multiprocessing; legacy local FITS output is available via
  `--output-format fits`, and Datalink remains FITS-only.
- Cutana segmentation-map shortcut now computes
  `floor(segmentation_map_id / 1_000_000)` locally and joins directly on
  `tile_index`, matching `query-segmap`.

## [0.2.0] - 2026-03-02

### Added
- `query-cutana` CLI command for generating Cutana input files from object IDs or coordinates.
- Canonical spectra compilation support for chunk-level parallel workers via `compile-spectra --workers`.
- Datalink spectra compilation mode with `compile-spectra --use-datalink`.
- Sphinx + Read the Docs documentation scaffold and user guide pages.

### Changed
- Default spectra chunk size changed to `1000` extensions per output FITS.
- Spectra compilation now resumes by default and can continue correctly across previous runs with different chunk sizes.
- Crossmatch now includes `point_like_prob` and `extended_prob` in outputs.
- Credential setup guidance updated to use secure user-local files under `~/.euclidkit`.
- Project license updated to BSD 3-Clause for public release/distribution.

### Fixed
- Crossmatch async mode now attempts result download and writes table output when available.
- Crossmatch object-id mode logging now reports mode-specific messaging (no radius confusion).
- Metadata table generation in spectra compilation fixed for mixed chunk/resume scenarios.

## [0.2.0rc.4] - 2026-02-23

### Changed
- Resume logic for canonical spectra compilation now inspects existing chunk FITS files and skips rows based on actual compiled extension counts.
- Resume now supports differing historical chunk sizes by continuing from already-compiled row count rather than assuming the current chunk size.
- Metadata generation now maps chunk/extension assignments from real output files, improving correctness for mixed-size resume cases.

## [0.2.0rc.3] - 2026-02-22

### Changed
- README cleaned for public release: removed unimplemented feature claims and examples.
- Added a project changelog file and release entries.

## [0.2.0rc.2] - 2026-02-22

### Added
- `compile-spectra --workers` option for chunk-level parallel compilation in canonical (local FITS) mode.
- Resume behavior for compilation by default when output chunks already exist.
- Datalink spectrum retrieval path using private astroquery client (`__eucliddata.load_data`) for improved robustness.
- `compile-spectra --limit` option to process only the first N rows for quick tests.
- Additional CLI tests for datalink mode, resume behavior, and worker forwarding.

### Changed
- Spectrum compilation now logs grouping statistics for FITS-path grouping per chunk.
- Compiled FITS extension header key standardized to `SOURC_ID` for consistency with archive datalink FITS headers.
- Metadata output typing hardened to avoid object-dtype FITS write failures.
- `compile-spectra` output directory handling now prefers in-place resume and only prompts on explicit overwrite requests.

### Fixed
- Datalink handling for empty/invalid downloads after server-side 500 responses.
- Metadata table generation error: `chunk_file` object dtype caused FITS serialization failures.

## [0.2.0rc.1] - 2026-02-21

### Added
- New `query-cutana` CLI command to build Cutana input tables from source tables.
- Support for IDR field selection (`WIDE`/`DEEP`) in Cutana-related archive lookups.
- Improved async behavior for large crossmatch queries.
- Initial PyPI publishing workflow via GitHub Actions.

### Changed
- Cutout defaults updated to `15` arcsec.
- Query logic updated to rely on `segmentation_map_id` workflows for mosaic product retrieval.
- Environment-dependent table prefix handling refined (`q1`/`dr1` behavior across PDR/REG/IDR).

### Fixed
- Removed legacy distance-based matching fields and unnecessary `DISTANCE(...)` projections in cutout query paths.
- Corrected join strategy to use object-ID joins against MER catalog where applicable.

[Unreleased]: https://github.com/rudolffu/euclidkit/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/rudolffu/euclidkit/compare/v0.2.0rc.4...v0.2.0
[0.2.0rc.4]: https://github.com/rudolffu/euclidkit/compare/v0.2.0rc.3...v0.2.0rc.4
[0.2.0rc.3]: https://github.com/rudolffu/euclidkit/compare/v0.2.0rc.2...v0.2.0rc.3
[0.2.0rc.2]: https://github.com/rudolffu/euclidkit/compare/v0.2.0rc.1...v0.2.0rc.2
[0.2.0rc.1]: https://github.com/rudolffu/euclidkit/releases/tag/v0.2.0rc.1
