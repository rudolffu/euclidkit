# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/rudolffu/euclidkit/compare/v0.2.0rc.3...HEAD
[0.2.0rc.3]: https://github.com/rudolffu/euclidkit/compare/v0.2.0rc.2...v0.2.0rc.3
[0.2.0rc.2]: https://github.com/rudolffu/euclidkit/compare/v0.2.0rc.1...v0.2.0rc.2
[0.2.0rc.1]: https://github.com/rudolffu/euclidkit/releases/tag/v0.2.0rc.1
