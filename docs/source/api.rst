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

Segmentation-map metadata lookup:

.. code-block:: python

   segmaps = archive.query_segmentation_maps(
       source_table=results,
       output_file="segmentation_maps.fits",
   )

``query_segmentation_maps`` requires ``SEGMENTATION_MAP_ID``, ``object_id``,
and source coordinates. It prefers ``ra``/``dec`` columns and falls back to
``mer_ra``/``mer_dec`` when the former are absent. It computes ``tile_index`` as
``floor(SEGMENTATION_MAP_ID / 1_000_000)`` before joining to the environment
segmentation-map table. See :doc:`segmentation_maps` for table mapping and CLI
examples.

Compile local segmentation-map cutouts from those query results:

.. code-block:: python

   from euclidkit.core.segmap import compile_segmap_cutouts

   stats = compile_segmap_cutouts(
       input_table="segmentation_maps.fits",
       output_dir="./segmap_cutouts",
       size_arcsec=10.0,
   )

``compile_segmap_cutouts`` opens local ``datalabs_path`` + ``file_name`` FITS
files with memmap, writes one raw-label FITS cutout per input row, and returns
``SegmapCutoutStats`` with requested, written, failed, skipped-existing, and
output-file counts. See :doc:`segmentation_maps` for required columns and output
semantics.

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
       environment="IDR",
       credentials_file="~/.euclidkit/.cred.txt",
   )

``dithers_to_parquet`` returns ``DithersParquetStats`` with object counts,
combined/dither row counts, raw-frame metadata match counts, output part paths,
manifest path, and optional failures JSONL path. It automatically queries
``q1.raw_frame`` for PDR/Q1, ``dr1.raw_frame`` for IDR/DR1, or
``sedm.raw_frame`` for OTF/REG to annotate per-dither rows with observation
metadata. ``SpectrumCompiler`` remains available for legacy multi-extension
FITS compilation and Datalink FITS outputs.
