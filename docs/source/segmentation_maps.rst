Segmentation-Map Cutouts
========================

A MER segmentation map labels each pixel of a tile with the object that it
belongs to. A cutout of that map around your source shows the footprint of the
source and of its neighbours. You can use it to mask neighbours or to check
deblending.

The workflow has two steps:

1. ``query-segmap`` asks the archive which segmentation-map file covers each
   source. It runs on any computer.
2. ``compile-segmap`` cuts the maps out of the local files. It runs on ESA
   Datalabs and makes no archive queries.

Before you start
----------------

You need a ``crossmatch`` output table. ``query-segmap`` reads these columns
from it (names are matched without regard to case):

- ``segmentation_map_id``
- ``object_id``
- ``ra`` and ``dec``, or ``mer_ra`` and ``mer_dec`` if those are absent

If your table has no ``segmentation_map_id`` column, run
:ref:`crossmatch <crossmatch-guide>` first.

Step 1: Find the segmentation-map files
---------------------------------------

.. code-block:: bash

   euclidkit query-segmap -i matches.fits -o segmentation_maps.fits

The output adds the file location (``datalabs_path`` and ``file_name``) of the
tile for each source.

``--field`` and ``--deep-partition`` do not apply here, because the
segmentation-map table is not split by field. Use ``--release Q1`` for Q1
sources.

Step 2: Make the cutouts (on Datalabs)
--------------------------------------

.. code-block:: bash

   euclidkit compile-segmap -i segmentation_maps.fits -o ./segmap_cutouts \
     --size-arcsec 10

This writes one FITS file per source row, named:

.. code-block:: text

   segmap_cutout_row######_obj*_seg*.fits

Options:

- ``--size-arcsec`` is the side of the square cutout in arcsec (default 10).
- ``--overwrite`` replaces existing cutouts. By default they are skipped.
- ``--on-error skip`` continues past missing files and bad rows. The default,
  ``fail``, stops at the first error.
- ``--no-progress`` hides the progress display.

Read the cutouts
----------------

The primary HDU of each cutout holds the **raw segmentation labels**, not a
binary mask. Pixels of neighbouring objects keep their own labels. Where a
cutout extends past the edge of the tile, the missing pixels are zero.

For how a label value relates to an object, see the MER segmentation-map
product description in the
`Euclid DR1 Data Product Description <https://euclid.esac.esa.int/dr/dr1/dpdd/>`_.

.. code-block:: python

   from astropy.io import fits

   labels = fits.getdata("segmap_cutouts/<cutout file>.fits")
   print(labels.shape, labels.max())

Python API
----------

.. code-block:: python

   from euclidkit.core.data_access import EuclidArchive
   from euclidkit.core.segmap import compile_segmap_cutouts
   from euclidkit.utils.io import load_table

   archive = EuclidArchive()
   archive.login()
   archive.query_segmentation_maps(
       source_table=load_table("matches.fits"),
       output_file="segmentation_maps.fits",
   )
   archive.logout()

   stats = compile_segmap_cutouts(
       input_table="segmentation_maps.fits",
       output_dir="./segmap_cutouts",
       size_arcsec=10.0,
       on_error="skip",
   )
   print(stats.written_rows, "written,", stats.failed_rows, "failed")

How the tile is identified
--------------------------

``euclidkit`` computes the tile of each source from its segmentation-map ID:

.. code-block:: text

   tile_index = floor(segmentation_map_id / 1,000,000)

It then looks up that tile in the ``mer_segmentation_map`` table of the
selected release. :doc:`environments` lists the table for each release.
