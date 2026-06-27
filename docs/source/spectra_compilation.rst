Spectra Compilation
===================

Typical workflow
----------------

A common spectra workflow is:

.. code-block:: bash

   euclidkit crossmatch \
     --input my_sources.fits \
     --output crossmatch_results.fits \
     --environment IDR \
     --idr-field WIDE

   euclidkit query-spectra \
     --crossmatch crossmatch_results.fits \
     --output spectra_sources.fits \
     --environment IDR \
     --idr-field WIDE

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix raw_spectra \
     -L RGS

``query-spectra`` is the bridge between catalogue matching and spectrum file
access. It reads Euclid ``object_id`` values from the crossmatch table, joins
them to the archive spectra-source table through ``source_id``, and writes a
spectra-source catalogue with the local Datalabs file information needed by
``compile-spectra``.

``query-spectra`` matching
--------------------------

``query-spectra`` matches by object ID only. It does not perform positional
matching with RA/Dec.

- The input table must contain ``object_id``. If that column is absent,
  ``object_id_euclid`` is accepted as a fallback.
- Unique IDs are uploaded as a temporary ``object_id`` column.
- The archive query joins ``spectra_source.source_id = uploaded.object_id``.
- Rows without ``datalabs_path`` are excluded from the output.

The output contains the columns needed by local and Datalink compilation,
including ``source_id``, ``object_id``, ``ra_obj``, ``dec_obj``,
``datalabs_path``, ``file_name``, and ``hdu_index``.

Default Parquet mode (local Datalabs FITS paths)
------------------------------------------------

Non-Datalink ``compile-spectra`` is a local Datalabs workflow. It reads each row
from the spectra-source catalogue, combines ``datalabs_path`` and ``file_name``
to locate the FITS file, opens the requested ``hdu_index``, and writes the raw
spectrum arrays to chunked Parquet files. This is the default output mode.

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
     --lambda-range RGS

This writes:

- ``raw_sir_combined_part001.parquet``, ...
- ``raw_sir_dithers_part001.parquet``, ...
- ``raw_sir_manifest.json``

Use ``--dithers-only`` when combined spectra were already exported with
``compile-spectra`` Parquet mode.

Notes
-----

- Use ``--limit N`` with ``compile-spectra`` for quick Parquet or FITS tests.
- Existing Parquet part, manifest, or failures files require ``--overwrite``.
- ``--on-error fail`` stops on the first unreadable Parquet row.
- ``--on-error skip`` continues past unreadable rows and writes a failures JSONL.
- Local Datalabs modes require valid ``datalabs_path``, ``file_name``, and
  ``hdu_index`` values in the spectra-source catalogue.
