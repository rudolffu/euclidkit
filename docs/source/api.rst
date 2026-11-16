Python API
==========

Everything that the commands do is available from Python. This page lists the
main entry points. The guides explain the behaviour behind each one.

Archive queries return ``astropy.table.Table`` objects. Call
``table.to_pandas()`` if you prefer a pandas DataFrame.

Connect to the archive
----------------------

.. code-block:: python

   from euclidkit.core.data_access import EuclidArchive

   archive = EuclidArchive()      # public DR1
   archive.login()                # reads ~/.euclidkit/.cred.txt
   try:
       ...                        # your queries
   finally:
       archive.logout()

``EuclidArchive(environment="PDR", data_release=None)``
   ``data_release`` is ``"DR1"`` (default) or ``"Q1"`` for the public archive.
   Leave it unset for ``environment="IDR"``, ``"OTF"``, or ``"REG"``. See
   :doc:`environments`.

``login(credentials_file=None, user=None, password=None, use_keyring=True, prompt=False, store_in_keyring=False)``
   See :doc:`configuration` for the order in which credentials are looked up.

Many query methods take ``idr_field`` (``"WIDE"`` or ``"DEEP"``) and
``idr_deep_partition`` (``"survey"``, ``"mode"``, or ``"both"``). Despite the
names, they apply to public DR1 as well. They correspond to ``--field`` and
``--deep-partition``.

Catalogue queries
-----------------

These are methods of ``EuclidArchive``. Guide: :doc:`catalogues`.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Method
     - Purpose
   * - ``crossmatch_sources``
     - Match a local table, DataFrame, or file to the MER catalogue
   * - ``upload_user_table``
     - Upload a table to your archive user space
   * - ``crossmatch_user_table``
     - Match a table that is already in your archive user space
   * - ``query_morphology``
     - Add MER morphology by object ID
   * - ``query_zspe_candidates``
     - Add SPE redshift candidates by object ID
   * - ``query_spectra_sources``
     - Find the spectra of your sources
   * - ``query_segmentation_maps``
     - Find the segmentation-map file for each source

.. code-block:: python

   from astropy.table import Table

   matches = archive.crossmatch_sources(
       user_table="my_sources.csv",   # a path, a Table, or a DataFrame
       radius=1.0,                    # arcsec
       output_file="matches.fits",
       drop_empty_columns=True,
   )

   deep = archive.crossmatch_sources(
       user_table="deep_sources.fits",
       idr_field="DEEP",
       idr_deep_partition="both",
   )

   morphology = archive.query_morphology(
       Table.read("matches.fits"),
       output_file="morphology.fits",
   )

   redshifts = archive.query_zspe_candidates(
       crossmatch_table=matches,
       object_type="qso",
       output_file="redshifts.fits",
   )

   spectra_sources = archive.query_spectra_sources(
       crossmatch_table=matches,
       match_mode="auto",             # or "object-id", "spatial"
       output_file="spectra_sources.fits",
   )

Things to know:

- ``crossmatch_sources`` matches by position unless the table has an ID
  column. Pass ``use_object_id=True`` or ``False`` to force one.
- With ``full_async=True``, ``crossmatch_sources``, ``query_morphology``, and
  ``query_zspe_candidates`` run asynchronous jobs. ``query_morphology`` then
  requires ``output_file``, and returns job or manifest information instead of
  a table.
- ``query_morphology`` and ``query_zspe_candidates`` raise ``ValueError`` for
  Q1.

Footprint selection
-------------------

.. code-block:: python

   from euclidkit.core.footprint import filter_catalog_by_moc

   in_wide = filter_catalog_by_moc(
       "my_sources.fits", "in_wide.fits", survey="WIDE",
   )

The function writes the filtered table to the second argument and also returns
it. ``list_available_surveys()`` in the same module returns the packaged
footprint names.

Spectra export
--------------

Guide: :doc:`spectra_compilation`. These functions read the Datalabs volumes.

.. code-block:: python

   from euclidkit.core.spectra_parquet import spectra_to_parquet, dithers_to_parquet

   stats = spectra_to_parquet(
       catalog_table="spectra_sources.fits",
       output_prefix="./output/raw_spectra",
       chunk_size=2000,
       workers=8,
       lambda_range="RGS",
       on_error="skip",
   )
   print(stats.exported_rows, stats.output_files, stats.manifest_path)

   dither_stats = dithers_to_parquet(
       catalog_table="spectra_sources.fits",
       output_prefix="./output/raw_sir",
       lambda_range="RGS",
       include_combined=True,
       environment="PDR",
       data_release="DR1",
   )
   print(dither_stats.dither_rows, dither_stats.dither_output_files)

``spectra_to_parquet`` returns ``RawParquetStats``. ``dithers_to_parquet``
returns ``DithersParquetStats``, which also reports how many dither rows were
matched to a raw frame.

The catalogue table must have the columns ``datalabs_path``, ``file_name``,
``hdu_index``, ``source_id``, ``object_id``, ``ra_obj``, and ``dec_obj``. The
output of ``query_spectra_sources`` has them.

For FITS output and Datalink downloads, use
``euclidkit.core.spectra.SpectrumCompiler`` (``compile_spectra`` and
``compile_spectra_datalink``). To read one spectrum as a FITS HDU, use
``EuclidArchive.get_individual_spectrum(datalabs_path, file_name, hdu_index)``.

Image cutouts
-------------

Guide: :doc:`cutouts`.

.. code-block:: python

   from euclidkit.core.cutouts import CutoutGenerator

   cutouts = CutoutGenerator(archive=archive)

   # Prepare a Cutana catalogue
   cutouts.generate_cutana_input(
       sources="my_sources.fits",
       output_file="cutana_input.csv",
       instrument_name="VIS",
       cutout_size="arcsec",
       cutout_size_value=15.0,
   )

   # Or write cutouts directly (on Datalabs)
   n_written = cutouts.make_cutouts(
       {"ra": 150.0, "dec": 2.5},
       "./cutouts",
       cutout_size=10,                # arcsec
       instrument_name="VIS",
   )

Segmentation-map cutouts
------------------------

Guide: :doc:`segmentation_maps`.

.. code-block:: python

   from euclidkit.core.segmap import compile_segmap_cutouts

   segmaps = archive.query_segmentation_maps(
       source_table=matches,
       output_file="segmentation_maps.fits",
   )
   stats = compile_segmap_cutouts(
       input_table="segmentation_maps.fits",
       output_dir="./segmap_cutouts",
       size_arcsec=10.0,
   )

``compile_segmap_cutouts`` returns ``SegmapCutoutStats`` with the counts of
requested, written, failed, and skipped rows.

Reading and writing tables
--------------------------

.. code-block:: python

   from euclidkit.utils.io import load_table, save_table

   table = load_table("matches.fits")      # FITS, CSV, VOTable, ECSV, Parquet
   save_table(table, "matches.parquet")    # format follows the extension
