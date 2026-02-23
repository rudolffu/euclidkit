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

Query spectra for crossmatched objects:

.. code-block:: bash

   euclidkit query-spectra \
     --crossmatch crossmatch_results.fits \
     --output spectra_sources.fits

Compile spectra into chunked FITS (on Datalabs):

.. code-block:: bash

   euclidkit compile-spectra \
     --spectra-table spectra_sources.fits \
     --output-dir ./output \
     --prefix compiled_spectra \
     --max-extensions 1000

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

   spectra = archive.query_spectra_sources(
       crossmatch_table=xmatch,
       output_file="spectra_sources.fits",
   )
