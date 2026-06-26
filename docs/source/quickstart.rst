Quickstart
==========

Basic CLI workflow
------------------

Crossmatch sources:

.. code-block:: bash

   euclidkit crossmatch \
     --input my_sources.csv \
     --output crossmatch_results.fits \
     --radius 1.0 \
     --environment PDR

For IDR DEEP MER data, the default partition is ``survey``. Use
``--idr-deep-partition mode`` for CDFS/COSMOS, or ``both`` to query survey
first and mode second:

.. code-block:: bash

   euclidkit crossmatch \
     --input my_deep_sources.fits \
     --output deep_crossmatch_results.fits \
     --environment IDR \
     --idr-field DEEP \
     --idr-deep-partition both

Query spectra for crossmatched objects:

.. code-block:: bash

   euclidkit query-spectra \
     --crossmatch crossmatch_results.fits \
     --output spectra_sources.fits

Export local Datalabs spectra to raw Parquet parts. This is the default
non-Datalink ``compile-spectra`` output; Datalink still writes FITS files.

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix raw_spectra \
     --chunk-size 2000 \
     -L RGS

Export per-dither spectra from the same local Datalabs catalog rows:

.. code-block:: bash

   euclidkit dithers-to-parquet \
     --catalog-table spectra_sources.fits \
     --output-prefix ./output/raw_sir \
     --lambda-range RGS

Compile both Euclid arms via datalink (RGS + BGS):

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix compiled_dl \
     --use-datalink \
     --environment IDR \
     -L BOTH

Python API example
------------------

.. code-block:: python

   from euclidkit.core.data_access import EuclidArchive

   archive = EuclidArchive(environment="PDR")
   archive.login(credentials_file="~/.euclidkit/.cred.txt")

   xmatch = archive.crossmatch_sources(
       user_table="my_sources.csv",
       radius=1.0,
       output_file="crossmatch_results.fits",
   )

   deep_xmatch = archive.crossmatch_sources(
       user_table="my_deep_sources.fits",
       output_file="deep_crossmatch_results.fits",
       idr_field="DEEP",
       idr_deep_partition="survey",
   )

   spectra = archive.query_spectra_sources(
       crossmatch_table=xmatch,
       output_file="spectra_sources.fits",
   )
