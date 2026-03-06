Spectra Compilation
===================

Canonical mode (local FITS paths on Datalabs)
---------------------------------------------

Canonical compilation reads source FITS files from ``datalabs_path`` + ``file_name``
columns and copies selected ``hdu_index`` rows into chunked outputs. This requires 
access to the Datalabs volumes.

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix compiled_spectra \
     --max-extensions 1000

For IDR DEEP canonical tables containing both BGS and RGS rows, select the arm
with ``-L/--lambda-range``. In ``BOTH`` mode, euclidkit writes separate outputs
with ``_rgs`` and ``_bgs`` suffixes.

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix compiled_deep \
     --environment IDR \
     --idr-field DEEP \
     -L BOTH

Resume behavior
---------------

Resume is the default when output directory already has chunk files and
``--overwrite`` is not provided.

During resume, euclidkit:

1. Inspects existing contiguous chunk files.
2. Counts compiled spectra from actual extension counts.
3. Skips already compiled input rows.
4. Continues with remaining rows, even if previous runs used a different chunk size.

Parallel workers (not recommended on Datalabs)
----------------------------------------------

Use ``--workers`` for chunk-level parallelism in canonical mode:

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --workers 2

On shared Datalabs storage, ``--workers 2`` is often not faster due to I/O
contention. Benchmark on your setup; ``--workers 1`` is usually the safest default.

Datalink mode
-------------

Use ``--use-datalink`` to retrieve spectra by source ID instead of local
``datalabs_path`` file access.

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --use-datalink \
     --environment IDR \
     --schema sedm \
     --retrieval-type SPECTRA_RGS

Notes:

- ``--workers`` currently applies to canonical mode only.
- For quick tests, use ``--limit N``.
- ``-L/--lambda-range`` is the unified arm selector for both canonical and
  datalink modes. In datalink mode, euclidkit maps it to
  ``--retrieval-type`` (``RGS -> SPECTRA_RGS``, ``BGS -> SPECTRA_BGS``,
  ``BOTH -> ALL``). ``--retrieval-type`` is kept for backward compatibility.
