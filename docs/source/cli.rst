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
