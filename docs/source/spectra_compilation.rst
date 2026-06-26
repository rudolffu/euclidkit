Spectra Compilation
===================

Default Parquet mode (local Datalabs FITS paths)
------------------------------------------------

``compile-spectra`` reads catalog rows with ``datalabs_path``, ``file_name``,
``hdu_index``, ``source_id``, ``object_id``, ``ra_obj``, and ``dec_obj``. By
default it reads the requested FITS HDU directly and writes raw spectra Parquet
parts.

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
columns from the FITS HDU, ``wavelength``, and ``flux`` as a ``signal`` alias.
Derived ``err``, ``valid``, and ``ivar`` columns are intentionally not written.

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
families. XML LambdaRange annotation is not used in Parquet mode.

Legacy local FITS mode
----------------------

Use ``--output-format fits`` to keep the previous multi-extension FITS compiler.
This path still supports IDR DEEP XML LambdaRange annotation and FITS metadata
output.

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix compiled_spectra \
     --output-format fits \
     --max-extensions 1000

Resume behavior applies to this FITS compatibility mode when output chunk files
already exist and ``--overwrite`` is not provided.

Datalink mode
-------------

Use ``--use-datalink`` to retrieve spectra by source ID instead of local
``datalabs_path`` file access. Datalink remains FITS-only and ignores the local
``--output-format`` default.

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --use-datalink \
     -L BOTH \
     --environment IDR \
     --schema sedm

In datalink ``-L BOTH`` mode, euclidkit writes separate outputs:

- ``<prefix>_rgs_chunk_###.fits``
- ``<prefix>_bgs_chunk_###.fits``

Dithers Parquet Export
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

- Parquet workers default to ``min(os.cpu_count(), 8)``.
- FITS compatibility and Datalink modes default to one worker.
- Use ``--limit N`` with ``compile-spectra`` for quick Parquet or FITS tests.
- Use ``--on-error skip`` to continue past unreadable rows and write a failures JSONL.
