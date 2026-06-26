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
from catalog rows. Required catalog columns are ``datalabs_path``, ``file_name``,
``hdu_index``, ``source_id``, ``object_id``, ``ra_obj``, and ``dec_obj``.

.. code-block:: python

   from euclidkit.core.spectra_parquet import spectra_to_parquet

   stats = spectra_to_parquet(
       catalog_table="spectra_sources.fits",
       output_prefix="./output/raw_spectra",
       chunk_size=2000,
       workers=8,
       lambda_range="RGS",
       on_error="skip",
   )

``spectra_to_parquet`` returns ``RawParquetStats`` with row counts, output part
paths, manifest path, and optional failures JSONL path.

Use ``dithers_to_parquet`` to export the combined spectrum and matching
``*_DITH1D_*_SIGNAL`` spectra for each object:

.. code-block:: python

   from euclidkit.core.spectra_parquet import dithers_to_parquet

   dither_stats = dithers_to_parquet(
       catalog_table="spectra_sources.fits",
       output_prefix="./output/raw_sir",
       chunk_size=2000,
       workers=8,
       lambda_range="RGS",
       include_combined=True,
   )

``dithers_to_parquet`` returns ``DithersParquetStats`` with object counts,
combined/dither row counts, output part paths, manifest path, and optional
failures JSONL path. ``SpectrumCompiler`` remains available for legacy
multi-extension FITS compilation and Datalink FITS outputs.
