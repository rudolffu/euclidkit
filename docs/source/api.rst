API Overview
============

Core classes
------------

EuclidArchive
~~~~~~~~~~~~~

Main interface for archive login, TAP operations, crossmatch, and spectra query.

Typical usage:

.. code-block:: python

   from euclidkit.core.data_access import EuclidArchive

   archive = EuclidArchive(environment="PDR")
   archive.login(credentials_file="~/.euclidkit/.cred.txt")

   results = archive.crossmatch_sources(
       user_table="sources.csv",
       radius=1.0,
       output_file="results.fits",
   )

IDR DEEP MER partition selection:

.. code-block:: python

   archive = EuclidArchive(environment="IDR")
   archive.login()

   # Default DEEP partition is "survey" (EDFN, EDFF, EDFS).
   deep_survey = archive.crossmatch_sources(
       user_table="deep_sources.fits",
       output_file="deep_survey_crossmatch.fits",
       idr_field="DEEP",
       idr_deep_partition="survey",
   )

   # Use "mode" for CDFS/COSMOS, or "both" to query survey first then mode.
   deep_both = archive.crossmatch_sources(
       user_table="deep_sources.fits",
       output_file="deep_both_crossmatch.fits",
       idr_field="DEEP",
       idr_deep_partition="both",
   )

The same ``idr_deep_partition`` argument is available on
``crossmatch_user_table`` and Cutana input generation paths that resolve MER
metadata. It applies only to MER catalogue selection, not spectra-source or
SPE redshift candidate table selection.

SpectrumCompiler
~~~~~~~~~~~~~~~~

Compiles queried spectra into chunked FITS products and metadata tables.

.. code-block:: python

   from euclidkit.core.spectra import SpectrumCompiler

   compiler = SpectrumCompiler(max_extensions=1000)
   output_files = compiler.compile_spectra(
       spectra_table=spectra_table,
       output_dir="./output",
       output_prefix="compiled_spectra",
       workers=1,
   )

   metadata_file = compiler.create_metadata_table(
       spectra_table=spectra_table,
       output_files=output_files,
       output_dir="./output",
   )
