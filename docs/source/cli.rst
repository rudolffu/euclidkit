CLI Reference
=============

Top-level commands
------------------

- ``euclidkit init-config``
- ``euclidkit diagnostics``
- ``euclidkit crossmatch``
- ``euclidkit query-spectra``
- ``euclidkit query-zspe``
- ``euclidkit compile-spectra``
- ``euclidkit dithers-to-parquet``
- ``euclidkit query-cutana``
- ``euclidkit cutouts``
- ``euclidkit query-segmap``
- ``euclidkit compile-segmap``
- ``euclidkit upload-table``
- ``euclidkit select-footprint``

Use ``--help`` on any command for full options:

.. code-block:: bash

   euclidkit query-spectra --help

Environment options
-------------------

Commands that access the Euclid archive support:

- ``--environment``: ``PDR``, ``IDR``, ``OTF``, ``REG``
- ``--idr-field`` (IDR-only commands): ``WIDE`` or ``DEEP``
- ``--idr-deep-partition`` for MER-based IDR DEEP workflows:
  ``survey`` (default), ``mode``, or ``both``

crossmatch
----------

Crossmatch a user table against Euclid MER sources.

Example:

.. code-block:: bash

   euclidkit crossmatch \
     --input my_sources.fits \
     --output crossmatch_results.fits \
     --radius 1.0 \
     --environment IDR \
     --idr-field WIDE

IDR DEEP partition example:

.. code-block:: bash

   euclidkit crossmatch \
     --input my_sources.fits \
     --output deep_mode_crossmatch.fits \
     --match-mode object-id \
     --environment IDR \
     --idr-field DEEP \
     --idr-deep-partition mode

Archive user-table example (no re-upload):

.. code-block:: bash

   euclidkit crossmatch \
     --user-table-name my_table \
     --output crossmatch_results.fits \
     --match-mode object-id \
     --environment IDR \
     --idr-field WIDE

Large-table async example:

.. code-block:: bash

   euclidkit crossmatch \
     --input huge_sources.fits \
     --output huge_crossmatch.fits \
     --match-mode object-id \
     --full-async \
     --async-chunk-size 500000

Drop empty result columns:

.. code-block:: bash

   euclidkit crossmatch \
     --input my_sources.fits \
     --output crossmatch_results.fits \
     --drop-empty-columns

Option semantics:

- ``--max-sources`` limits total processed rows from the input table.
- ``--async-chunk-size`` controls rows per async TAP job in ``--full-async`` mode.
- ``--drop-empty-columns`` removes columns where every final result value is
  null or missing before saving ``--output``. Zero, ``False``, and empty-string
  columns are retained. Async chunk part files remain raw.
- ``--idr-deep-partition`` is used only with ``--environment IDR --idr-field DEEP``
  for MER-based commands. ``survey`` queries
  ``catalogue.mer_catalogue_deep_survey`` (EDFN, EDFF, EDFS), ``mode`` queries
  ``catalogue.mer_catalogue_deep_mode`` (CDFS, COSMOS), and ``both`` queries
  survey first and mode second. The default is ``survey``.

``--full-async`` behavior:

- For smaller inputs, ``euclidkit`` submits one async TAP job, downloads the result to the requested output file, and removes the remote job after saving.
- For large local input tables supplied with ``--input``, ``euclidkit`` splits the upload into async chunks, writes chunk files named ``<output>_part_####.fits``, writes ``<output>.manifest.json``, removes each remote job after its chunk is saved, and merges the chunk files into the requested final output.
- For large archive user tables supplied with ``--user-table-name``, ``euclidkit`` uses the same chunk-file plus manifest workflow before producing the final merged FITS file.

Matching mode recommendation:

- Prefer ``--match-mode object-id`` whenever the input already contains Euclid ``object_id`` values, or ``source_id`` values that should be joined to MER ``object_id``. This avoids positional matching and is usually faster and more robust for large tables.

query-spectra
-------------

Query spectra-source rows for objects in an object-ID or coordinate table.
``--match-mode`` accepts the same values as ``crossmatch``:
``auto`` (default), ``object-id``, or ``spatial``.

.. code-block:: bash

   euclidkit query-spectra \
     --crossmatch my_spectral_ids.fits \
     --output spectra_sources.fits \
     --environment IDR \
     --idr-field WIDE

Spatial example:

.. code-block:: bash

   euclidkit query-spectra \
     --crossmatch my_coordinates.fits \
     --output spectra_sources.fits \
     --match-mode spatial \
     --radius 1.0 \
     --ra-col ra \
     --dec-col dec \
     --environment IDR \
     --idr-field WIDE

Matching behavior:

- ``auto`` uses object-ID mode when ``object_id`` or ``object_id_euclid`` is
  present; otherwise it uses spatial mode when RA/Dec columns can be resolved.
- Object-ID mode uploads unique IDs as a temporary ``object_id`` column and
  joins ``spectra_source.source_id = uploaded.object_id``.
- A MER crossmatch table is not required. Any local table with the needed
  spectra-source IDs or coordinates can be used. If your IDs are in a
  ``source_id`` column, rename or copy that column to ``object_id`` before
  running ``query-spectra``.
- Spatial mode matches ``ra_obj``/``dec_obj`` to the input coordinates within
  ``--radius`` arcsec and keeps only the nearest spectra-source row per input
  row.
- Spatial coordinate columns are resolved using exact ``--ra-col``/``--dec-col``
  names first, then case-insensitive matches and common aliases such as
  ``RA``/``DEC``, ``right_ascension``/``declination``,
  ``ra_deg``/``dec_deg``, ``RAJ2000``/``DEJ2000``, and Euclid-specific
  ``ra_euclid``/``dec_euclid`` names.
- Rows with missing ``datalabs_path`` are excluded because local Datalabs
  compilation needs a file path.

Output columns include ``source_id``, ``ra_obj``, ``dec_obj``,
``datalabs_path``, ``file_name``, ``hdu_index``, and the uploaded
``object_id``. Spatial mode also includes ``input_row_id``, ``input_ra``,
``input_dec``, and ``separation_arcsec``; its ``object_id`` is set from the
matched ``source_id`` for downstream ``compile-spectra`` compatibility. For
IDR, ``--idr-field WIDE`` queries ``catalogue.spectra_source_wide`` and
``--idr-field DEEP`` queries ``catalogue.spectra_source_deep``. Other
environments use the corresponding ``spectra_source`` table for the selected
release.

The resulting table is the usual input to ``compile-spectra``. See
:doc:`spectra_compilation` for the full spectra workflow.

query-zspe
----------

Query IDR SPE redshift candidates by ``object_id`` and join the matched SPE
columns back to the input table.

.. code-block:: bash

   euclidkit query-zspe \
     --crossmatch crossmatch_results.fits \
     --output zspe_matches.fits \
     --object-type qso \
     --idr-field WIDE

``query-zspe`` supports ``--object-type qso|galaxy``. WIDE queries search the
``wide_survey`` candidate table first and then use ``wide`` as the fallback for
remaining objects. Use ``--full-async`` for one server-side async query with the
same WIDE fallback behavior.

compile-spectra
---------------

Compile or export spectra rows returned by ``query-spectra``. Non-Datalink mode
is local Datalabs file access and defaults to raw Parquet output. Datalink mode
retrieves spectra through the archive service and remains FITS-only.

Default local Datalabs Parquet example:

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix raw_spectra \
     --chunk-size 2000 \
     -L RGS

Legacy local FITS example:

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix compiled_spectra \
     --output-format fits

Datalink dual-arm FITS example:

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix compiled_dl \
     --use-datalink \
     --environment IDR \
     --schema sedm \
     -L BOTH

Option semantics:

- Non-Datalink ``compile-spectra`` reads local ``datalabs_path`` + ``file_name``
  FITS files and writes Parquet by default.
- ``--output-format fits`` selects the legacy local multi-extension FITS compiler.
- ``--use-datalink`` retrieves spectra from the archive by source ID, remains
  FITS-only, and ignores ``--output-format parquet``.
- Parquet workers default to ``min(os.cpu_count(), 8)`` when ``--workers`` is omitted.
- FITS compatibility and Datalink paths default to one worker.
- ``--on-error fail`` stops on the first unreadable Parquet row; ``--on-error skip`` continues and writes ``<prefix>_failures.jsonl``.

See :doc:`spectra_compilation` for the detailed Datalabs, Datalink, arm
selection, and output-file behavior.

dithers-to-parquet
------------------

Export local Datalabs combined and per-dither SIR spectra to Parquet parts.
This command reads local ``datalabs_path`` + ``file_name`` entries and does not
use Datalink.

.. code-block:: bash

   euclidkit dithers-to-parquet \
     --catalog-table spectra_sources.fits \
     --output-prefix ./output/raw_sir \
     --lambda-range RGS \
     --workers 8 \
     --environment IDR

Per-dither parquet rows are automatically annotated with ``obs_time_mjd``,
``obs_time_utc``, and ``pa`` from the environment raw-frame table
(``q1.raw_frame`` for PDR/Q1, ``dr1.raw_frame`` for IDR/DR1, and
``sedm.raw_frame`` for OTF/REG), joined on ``pointing_id`` and
``grism_wheel_pos = gwa_pos``.

query-cutana
------------

Build Cutana input CSV from source rows containing object IDs or coordinates.

Example:

.. code-block:: bash

   euclidkit query-cutana \
     --sources my_sources.fits \
     --output cutana_input.csv \
     --instrument VIS \
     --cutout-size arcsec \
     --cutout-size-value 15

IDR DEEP example:

.. code-block:: bash

   euclidkit query-cutana \
     --sources my_sources.fits \
     --output cutana_deep.csv \
     --instrument NISP \
     --nisp-filters NIR_Y,NIR_H \
     --environment IDR \
     --idr-field DEEP \
     --idr-deep-partition both \
     --cutout-size arcsec \
     --cutout-size-value 15

cutouts
-------

Generate Euclid image cutouts from catalogue positions or object metadata.

.. code-block:: bash

   euclidkit cutouts --help

query-segmap
------------

Query MER segmentation-map metadata for rows that already contain
``SEGMENTATION_MAP_ID`` from MER crossmatch output.

Example:

.. code-block:: bash

   euclidkit query-segmap \
     --input crossmatch_results.fits \
     --output segmentation_maps.fits \
     --environment IDR

``query-segmap`` computes ``tile_index`` locally as
``floor(SEGMENTATION_MAP_ID / 1_000_000)`` and joins that value to the
environment-specific segmentation-map table. See :doc:`segmentation_maps` for
required columns, environment mapping, and downstream cutout compilation.

compile-segmap
--------------

Create local FITS cutouts from ``query-segmap`` output rows.

Example:

.. code-block:: bash

   euclidkit compile-segmap \
     --input segmentation_maps.fits \
     --output-dir ./segmap_cutouts \
     --size-arcsec 10

``compile-segmap`` opens local ``datalabs_path`` + ``file_name`` FITS files and
writes one raw-label FITS cutout per source row. See :doc:`segmentation_maps`
for required columns, output semantics, and error-handling options.

upload-table
------------

Upload a local table into Euclid TAP user workspace.

Example:

.. code-block:: bash

   euclidkit upload-table \
     --input my_sources.fits \
     --table-name my_sources_work

select-footprint
----------------

Filter an input catalogue to rows that fall within a packaged Euclid footprint
MOC.

.. code-block:: bash

   euclidkit select-footprint \
     --input my_sources.fits \
     --output in_footprint.fits \
     --survey mer-wide
