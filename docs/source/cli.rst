CLI Reference
=============

Top-level commands
------------------

- ``euclidkit init-config``
- ``euclidkit diagnostics``
- ``euclidkit crossmatch``
- ``euclidkit query-spectra``
- ``euclidkit query-cutana``
- ``euclidkit compile-spectra``
- ``euclidkit upload-table``
- ``euclidkit cutouts``
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

Option semantics:

- ``--max-sources`` limits total processed rows from the input table.
- ``--async-chunk-size`` controls rows per async TAP job in ``--full-async`` mode.
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

Query spectra-source rows for objects in a crossmatch table.

Example:

.. code-block:: bash

   euclidkit query-spectra \
     --crossmatch crossmatch_results.fits \
     --output spectra_sources.fits \
     --environment IDR \
     --idr-field WIDE

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

compile-spectra
---------------

Compile spectra into chunked multi-extension FITS outputs.

Example:

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix compiled_spectra

IDR DEEP arm-selection example:

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix compiled_deep \
     --environment IDR \
     --idr-field DEEP \
     --lambda-range BOTH

Datalink dual-arm example:

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix compiled_dl \
     --use-datalink \
     --environment IDR \
     --schema sedm \
     -L BOTH

upload-table
------------

Upload a local table into Euclid TAP user workspace.

Example:

.. code-block:: bash

   euclidkit upload-table \
     --input my_sources.fits \
     --table-name my_sources_work
