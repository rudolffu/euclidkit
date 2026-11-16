CLI Reference
=============

Top-level commands
------------------

- ``euclidkit init-config``
- ``euclidkit diagnostics``
- ``euclidkit crossmatch``
- ``euclidkit query-morph``
- ``euclidkit query-spectra``
- ``euclidkit query-zspe``
- ``euclidkit compile-spectra``
- ``euclidkit dithers-to-parquet``
- ``euclidkit query-cutana``
- ``euclidkit download-cutana-alternatives``
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
- ``--release``: ``DR1`` (default) or ``Q1``, on PDR only
- ``--field`` (DR1/IDR catalogue commands): ``WIDE`` or ``DEEP``
- ``--deep-partition`` for MER-based DR1 DEEP workflows:
  ``survey`` (default), ``mode``, or ``both``

.. _query-morph:

Public release selection
------------------------

Archive commands accept ``--release DR1|Q1`` on PDR; DR1 is the default.
Use ``--environment IDR`` without ``--release`` for the independent internal
service. Legacy ``--idr-field`` and ``--idr-deep-partition`` names remain aliases
for ``--field`` and ``--deep-partition``. See :doc:`environments` for migration
and the public table mappings. Public-release-only operations unavailable in
Q1 fail explicitly rather than reading DR1.

query-morph
-----------

Query public DR1 MER morphology by object ID and inner join measurements onto matched
input rows. A prior MER crossmatch is unnecessary if Euclid ``object_id`` values
are already available. The key is case-insensitive; returned morphology names
are lowercase. Repeated input IDs retain separate rows. Existing columns with
the selected measurement names are replaced.

.. code-block:: bash

   euclidkit query-morph -i sources.fits -o morphology.fits
   euclidkit query-morph -i sources.fits -o deep_morphology.parquet \
       --field DEEP --deep-partition both --full-async

This command defaults to public DR1 and also supports explicit IDR access. Q1 is unsupported. WIDE (default) queries
``dr1.mer_morphology_wide_survey`` first, then
``dr1.mer_morphology_wide_mode`` for unmatched IDs. DEEP defaults to
``dr1.mer_morphology_deep_survey``; ``--deep-partition mode`` selects
``dr1.mer_morphology_deep_mode``. ``both`` queries both DEEP partitions
directly, keeping the survey measurement if an ID occurs in both.

The selection contains each of these columns and its ``_err`` column:

- ``concentration``, ``asymmetry``, ``smoothness`` (CAS)
- ``moment_20`` (M20), ``gini``
- ``sersic_sersic_vis_radius`` (arcsec)
- ``sersic_sersic_vis_axis_ratio``, ``sersic_sersic_vis_index``
- ``sersic_angle`` (degrees)
- ``flux_radius`` (pixels)

Other measurements are dimensionless. The DR1 FITS schema specifies 64-bit
integer IDs and 32-bit floating-point measurements; nulls are preserved.
See the `DR1 morphology schema <http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmdr1/merdpd/dpcards/mer_finalcatalog.html#morphology-catalog-fits-file>`_.

Sync queries batch 1,000 unique IDs. ``--full-async`` uses server-side WIDE
fallback and splits inputs above ``--async-chunk-size`` (default 500,000 unique
IDs). Raw chunk parts retain internal partition provenance; the final output
contains input columns and the 20 measurements/errors only. Failed downloads
write recovery metadata (``.job.json`` or ``.manifest.json``) without writing a
partial final table. Empty matches produce an empty table with the output schema.
Use ``--credentials/-c`` for a credentials file and ``--verbose/-v`` for logging.

crossmatch
----------

Crossmatch a user table against Euclid MER sources.

Example:

.. code-block:: bash

   euclidkit crossmatch \
     --input my_sources.fits \
     --output crossmatch_results.fits \
     --radius 1.0 \
     --environment PDR --release DR1 \
     --field WIDE

DR1 DEEP partition example:

.. code-block:: bash

   euclidkit crossmatch \
     --input my_sources.fits \
     --output deep_mode_crossmatch.fits \
     --match-mode object-id \
     --environment PDR --release DR1 \
     --field DEEP \
     --deep-partition mode

Archive user-table example (no re-upload):

.. code-block:: bash

   euclidkit crossmatch \
     --user-table-name my_table \
     --output crossmatch_results.fits \
     --match-mode object-id \
     --environment PDR --release DR1 \
     --field WIDE

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
- ``--deep-partition`` is used only with ``--environment PDR --release DR1 --field DEEP``
  for MER-based commands. ``survey`` queries
  ``dr1.mer_catalogue_deep_survey`` (EDFN, EDFF, EDFS), ``mode`` queries
  ``dr1.mer_catalogue_deep_mode`` (CDFS, COSMOS), and ``both`` queries
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
     --environment PDR --release DR1 \
     --field WIDE

Spatial example:

.. code-block:: bash

   euclidkit query-spectra \
     --crossmatch my_coordinates.fits \
     --output spectra_sources.fits \
     --match-mode spatial \
     --radius 1.0 \
     --ra-col ra \
     --dec-col dec \
     --environment PDR --release DR1 \
     --field WIDE

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
public DR1, ``--field WIDE`` queries ``dr1.spectra_source_wide`` and
``--field DEEP`` queries ``dr1.spectra_source_deep``. Other
environments use the corresponding ``spectra_source`` table for the selected
release.

The resulting table is the usual input to ``compile-spectra``. See
:doc:`spectra_compilation` for the full spectra workflow.

query-zspe
----------

Query public DR1 SPE redshift candidates by ``object_id`` and join the matched SPE
columns back to the input table.

.. code-block:: bash

   euclidkit query-zspe \
     --crossmatch crossmatch_results.fits \
     --output zspe_matches.fits \
     --object-type qso \
     --field WIDE

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
     --environment PDR --release DR1 \
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
     --environment PDR --release DR1

Per-dither parquet rows are automatically annotated with ``obs_time_mjd``,
``obs_time_utc``, and ``pa`` from the environment raw-frame table
(``dr1.raw_frame`` for PDR/DR1, ``q1.raw_frame`` for PDR/Q1, and
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

DR1 DEEP example:

.. code-block:: bash

   euclidkit query-cutana \
     --sources my_sources.fits \
     --output cutana_deep.csv \
     --instrument NISP \
     --nisp-filters NIR_Y,NIR_H \
     --environment PDR --release DR1 \
     --field DEEP \
     --deep-partition both \
     --cutout-size arcsec \
     --cutout-size-value 15

Cutana preparation products
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``query-cutana`` chooses one mosaic per source and band. It selects the nearest
image centre using spherical separation, then breaks ties by full file path
and product OID. Candidates without image-centre coordinates rank last. This
rule is deterministic, not an exposure-depth or image-quality ranking. Exact
duplicate file records are collapsed. Requested NISP band order is preserved;
the default is ``NIR_Y,NIR_J,NIR_H``. Hyphens in filter names are accepted.

The archive query requires the source centre to lie inside the mosaic footprint.
Keep Cutana's FITS checks enabled: archive footprints do not replace FITS/WCS
validation, and catalogue paths must be readable where Cutana runs.

For ``cutana_input.csv`` the command also writes:

* ``cutana_input.alternatives.csv``: all candidates for ambiguous source/band
  pairs, including the selected candidate, rank, distance, archive filename,
  paths and product OID when available.
* ``cutana_input.coverage.csv``: resolved sources, available/missing bands and
  whether a mosaic was found. Sources with no mosaic cannot be extracted.
* ``cutana_input.manifest.json``: band groups, counts and companion-file paths.

Sources missing some bands are retained. If available band sets differ, the
main CSV is a combined inventory and separate ``.bands-*.csv`` files are written.
**Run each manifest group separately in Cutana**, configuring
``selected_extensions`` and ``channel_weights`` for that group's available
bands. Cutana does not promise to fill missing bands: a fixed multiband channel
configuration can reject a source missing a configured channel. When all
sources share one band set, the manifest points to the original CSV.

Diameters must be finite and positive, at most 10000 pixels or 3600 arcseconds;
pixel diameters must be integers. Identical source/position rows collapse to
one output source; repeated object IDs at distinct positions receive distinct
stable IDs. Extra metadata uses the first input row for each object/position.
Coverage reports describe sources resolved through MER, not unresolved input IDs.

download-cutana-alternatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the unselected candidates later, without repeating the source query:

.. code-block:: bash

   euclidkit download-cutana-alternatives \
     --candidates cutana_input.alternatives.csv \
     --output-dir alternative_mosaics \
     --environment PDR --release DR1

Use the original query's archive environment, release and product service's schema
(the default follows the selected release; ``--schema`` overrides it). Each distinct alternative filename is downloaded once with
Astroquery's ``get_product(file_name=...)``. Existing nonempty files are reused;
new downloads are staged so interrupted transfers are not treated as complete.
The alternatives table can be filtered before downloading. This command does
not modify the prepared catalogue or automatically download mosaics during
``query-cutana``.

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
     --environment PDR --release DR1

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
