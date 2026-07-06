Spectra Compilation
===================

Typical workflow
----------------

The spectra workflow starts with an input table that contains either Euclid
spectral source IDs in an ``object_id`` column or sky coordinates. A MER
crossmatch table is one valid input, but it is not required.

If your table already contains an ``object_id`` column whose values correspond
to ``spectra_source.source_id``, query by object ID:

.. code-block:: bash

   euclidkit query-spectra \
     --crossmatch my_spectral_ids.fits \
     --output spectra_sources.fits \
     --match-mode object-id \
     --environment IDR \
     --idr-field WIDE

If your table only has positions, query by coordinates:

.. code-block:: bash

   euclidkit query-spectra \
     --crossmatch my_coordinates.fits \
     --output spectra_sources.fits \
     --match-mode spatial \
     --radius 1.0 \
     --environment IDR \
     --idr-field WIDE

Then compile the returned spectra-source catalogue:

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix raw_spectra \
     -L RGS

``query-spectra`` is the bridge between an ID/coordinate table and spectrum
file access. It queries the archive spectra-source table by either an uploaded
``object_id`` key joined to ``spectra_source.source_id`` or by matching input
coordinates to ``ra_obj``/``dec_obj``, then writes a spectra-source catalogue
with the local Datalabs file information needed by ``compile-spectra``. If you
do not already have this spectra-source catalogue, generate it with the
`query-spectra CLI
<https://euclidkit.readthedocs.io/en/latest/cli.html#query-spectra>`__ before
running ``compile-spectra``.

``query-spectra`` matching
--------------------------

``query-spectra`` supports ``--match-mode auto|object-id|spatial``. Auto mode
uses object-ID matching when ``object_id`` or ``object_id_euclid`` exists, and
otherwise uses spatial matching when RA/Dec columns can be resolved.

- Object-ID mode uploads unique IDs as a temporary ``object_id`` column and
  joins ``spectra_source.source_id = uploaded.object_id``.
- If your input column is named ``source_id`` rather than ``object_id``, rename
  or copy it to ``object_id`` before running ``query-spectra``; in this context
  ``object_id`` is the uploaded key used to join against
  ``spectra_source.source_id``.
- Spatial mode matches input coordinates to ``spectra_source.ra_obj`` and
  ``spectra_source.dec_obj`` within ``--radius`` arcsec, default ``1.0``.
- Spatial mode keeps only the nearest spectra-source row per input row.
- Spatial coordinate columns are resolved using exact ``--ra-col``/``--dec-col``
  names first, then case-insensitive matches and common aliases such as
  ``RA``/``DEC``, ``right_ascension``/``declination``,
  ``ra_deg``/``dec_deg``, ``RAJ2000``/``DEJ2000``, and
  ``ra_euclid``/``dec_euclid``.
- Rows without ``datalabs_path`` are excluded from the output.

The output contains the columns needed by local and Datalink compilation,
including ``source_id``, ``object_id``, ``ra_obj``, ``dec_obj``,
``datalabs_path``, ``file_name``, and ``hdu_index``. Spatial mode also includes
``input_row_id``, ``input_ra``, ``input_dec``, and ``separation_arcsec``; its
``object_id`` is set from the matched ``source_id`` so the result can be passed
directly to ``compile-spectra``.

Default Parquet mode (local Datalabs FITS paths)
------------------------------------------------

Non-Datalink ``compile-spectra`` is a local Datalabs workflow. It reads each row
from the spectra-source catalogue, combines ``datalabs_path`` and ``file_name``
to locate the FITS file, opens the requested ``hdu_index``, and writes the raw
spectrum arrays to chunked Parquet files. This is the default output mode.

If you do not already have ``spectra_sources.fits``, create it first from an
ID table or coordinate table:

.. code-block:: bash

   euclidkit query-spectra \
     --crossmatch my_spectral_ids_or_coordinates.fits \
     --output spectra_sources.fits \
     --environment IDR \
     --idr-field WIDE

Then export the local Datalabs spectra:

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix raw_spectra \
     --chunk-size 2000 \
     --workers 8 \
     -L RGS

This writes:

- ``raw_spectra_part001.parquet``, ...
- ``raw_spectra_manifest.json``
- ``raw_spectra_failures.jsonl`` only when ``--on-error skip`` records failures

Parquet rows include ``object_id``, ``source_id``, ``ra``, ``dec``, raw array
columns from the FITS HDU, ``wavelength``, and ``flux`` as an alias of
``signal``. Derived ``err``, ``valid``, and ``ivar`` columns are intentionally
not written.

The Parquet exporter is intended for bulk machine-learning or table-processing
workflows where direct array columns are easier to consume than multi-extension
FITS files. Worker count defaults to ``min(os.cpu_count(), 8)`` when
``--workers`` is omitted.

Arm selection
-------------

In Parquet mode, ``-L/--lambda-range`` filters using the FITS primary ``LRANGE``
header. ``-L BOTH`` runs separate RGS and BGS exports:

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix raw_spectra \
     -L BOTH

This writes separate ``raw_spectra_rgs_*`` and ``raw_spectra_bgs_*`` Parquet
families. XML LambdaRange annotation is not used in Parquet mode; the FITS
headers are the source of truth.

Legacy local FITS mode
----------------------

Use ``--output-format fits`` to keep the previous local multi-extension FITS
compiler. This path still reads local ``datalabs_path`` + ``file_name`` FITS
files, but writes FITS chunks instead of Parquet parts.

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix compiled_spectra \
     --output-format fits \
     --max-extensions 1000

The legacy FITS path supports IDR DEEP XML LambdaRange annotation and FITS
metadata output. Resume behavior applies when output chunk files already exist
and ``--overwrite`` is not provided. FITS compatibility mode defaults to one
worker.

Datalink mode
-------------

Use ``--use-datalink`` to retrieve spectra through the Euclid archive Datalink
service instead of reading local Datalabs files. This path uses the
spectra-source rows and source identifiers to request spectra from the archive
client, then writes multi-extension FITS chunks. It remains FITS-only and
ignores the local ``--output-format parquet`` default.

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix compiled_dl \
     --use-datalink \
     -L BOTH \
     --environment IDR \
     --schema sedm

In Datalink mode, ``-L RGS`` maps to ``SPECTRA_RGS`` retrieval,
``-L BGS`` maps to ``SPECTRA_BGS`` retrieval, and ``-L BOTH`` runs the two arms
separately. Dual-arm output writes separate FITS families:

- ``<prefix>_rgs_chunk_###.fits``
- ``<prefix>_bgs_chunk_###.fits``

``--retrieval-type`` is retained for compatibility, but ``-L/--lambda-range`` is
the preferred arm-selection interface. Datalink mode runs with one worker.

Dithers Parquet export
----------------------

``dithers-to-parquet`` exports the combined spectrum and matching
``*_DITH1D_*_SIGNAL`` spectra for each catalog object. It is local-Datalabs only
by data access model: it reads catalog ``datalabs_path`` + ``file_name`` entries
and does not use Datalink.

.. code-block:: bash

   euclidkit dithers-to-parquet \
     --catalog-table spectra_sources.fits \
     --output-prefix ./output/raw_sir \
     --chunk-size 2000 \
     --workers 8 \
     --lambda-range RGS \
     --raw-frame-table rawframe_by_pointing_id.vot.gz

This writes:

- ``raw_sir_combined_part001.parquet``, ...
- ``raw_sir_dithers_part001.parquet``, ...
- ``raw_sir_manifest.json``

Use ``--dithers-only`` when combined spectra were already exported with
``compile-spectra`` Parquet mode.

Raw-frame observation metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Per-dither parquet rows can be annotated with ``obs_time_mjd``,
``obs_time_utc``, and ``pa`` by passing ``--raw-frame-table``. The table must
contain ``pointing_id``, ``technique``, ``obs_time_mjd``, ``obs_time_utc``, and
``pa``. ``dithers-to-parquet`` filters it to ``technique = 'SPECTROIMAGE'`` and
matches ``pointing_id`` to the dither ID parsed from ``*_DITH1D_<n>_SIGNAL``;
if that dither ID is unavailable, it falls back to ``ptgid``.

For many pointings, avoid a long literal ID-list clause. Upload a distinct
one-column pointing table and join it to the archive raw-frame table:

.. code-block:: sql

   SELECT
     r.pointing_id,
     r.obs_time_mjd,
     r.obs_time_utc,
     r.pa
   FROM dr1.raw_frame AS r
   JOIN TAP_UPLOAD.dither_pointings AS p
     ON r.pointing_id = p.pointing_id
   WHERE r.technique = 'SPECTROIMAGE'

Use the raw-frame table/schema for the relevant release if you are not querying
DR1.

Notes
-----

- Use ``--limit N`` with ``compile-spectra`` for quick Parquet or FITS tests.
- Existing Parquet part, manifest, or failures files require ``--overwrite``.
- ``--on-error fail`` stops on the first unreadable Parquet row.
- ``--on-error skip`` continues past unreadable rows and writes a failures JSONL.
- Local Datalabs modes require valid ``datalabs_path``, ``file_name``, and
  ``hdu_index`` values in the spectra-source catalogue.
