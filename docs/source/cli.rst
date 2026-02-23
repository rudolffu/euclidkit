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

query-spectra
-------------

Query spectra-source rows for objects in a crossmatch table.

Example:

.. code-block:: bash

   euclidkit query-spectra \
     --crossmatch crossmatch_results.fits \
     --output spectra_sources.fits \
     --environment IDR

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

upload-table
------------

Upload a local table into Euclid TAP user workspace.

Example:

.. code-block:: bash

   euclidkit upload-table \
     --input my_sources.fits \
     --table-name my_sources_work

