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
       drop_empty_columns=True,
   )

``drop_empty_columns=True`` removes columns where every final crossmatch result
value is null or missing before saving and returning the table. Zero, ``False``,
and empty strings are preserved. For async crossmatch workflows, intermediate
part files are left unchanged and only the final merged output is pruned.

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

Spectra Parquet Export
~~~~~~~~~~~~~~~~~~~~~~

The default local Datalabs spectra workflow writes raw Parquet parts directly
from catalog rows with ``datalabs_path``, ``file_name``, and ``hdu_index``.

.. code-block:: python

   from euclidkit.core.spectra_parquet import spectra_to_parquet

   stats = spectra_to_parquet(
       catalog_table="spectra_sources.fits",
       output_prefix="./output/raw_spectra",
       chunk_size=2000,
       workers=8,
       lambda_range="RGS",
   )

Use ``dithers_to_parquet`` from the same module to export combined and per-dither
SIR spectra. ``SpectrumCompiler`` remains available for legacy multi-extension
FITS compilation and Datalink FITS outputs.
