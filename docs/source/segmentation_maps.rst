Segmentation Map Cutouts
========================

``euclidkit`` provides a two-step workflow for going from MER crossmatch rows to
local MER segmentation-map FITS cutouts:

1. ``query-segmap`` resolves the segmentation-map tile metadata from the archive.
2. ``compile-segmap`` opens the local tile FITS files and writes one cutout per
   source row.

This workflow is local-file oriented after ``query-segmap``: ``compile-segmap``
does not query TAP or use Datalink.

Workflow
--------

Start from a MER crossmatch table. The table must contain ``SEGMENTATION_MAP_ID``;
if it does not, rerun ``euclidkit crossmatch`` first.

.. code-block:: bash

   euclidkit query-segmap \
     --input crossmatch_results.fits \
     --output segmentation_maps.fits \
     --environment IDR

Then create cutouts from the local segmentation-map FITS files:

.. code-block:: bash

   euclidkit compile-segmap \
     --input segmentation_maps.fits \
     --output-dir ./segmap_cutouts \
     --size-arcsec 10

``compile-segmap`` groups rows by ``datalabs_path`` + ``file_name`` and opens
each large segmentation-map tile once with FITS memmap.

Required columns
----------------

``query-segmap`` input columns are matched case-insensitively:

- ``SEGMENTATION_MAP_ID``
- ``object_id``
- source coordinates: ``ra``/``dec``, or ``mer_ra``/``mer_dec`` if ``ra``/``dec`` are absent

``compile-segmap`` input is normally the output of ``query-segmap`` and requires:

- ``datalabs_path``
- ``file_name``
- ``segmentation_map_id``
- ``object_id``
- source coordinates: ``ra``/``dec``, or ``mer_ra``/``mer_dec`` if ``ra``/``dec`` are absent

Segmentation tile lookup
------------------------

``query-segmap`` computes the segmentation tile locally using the authoritative
formula:

.. code-block:: text

   tile_index = floor(SEGMENTATION_MAP_ID / 1_000_000)

It then joins the uploaded tile indices to the environment-specific
``mer_segmentation_map`` table:

- ``PDR``: ``q1.mer_segmentation_map``
- ``IDR``: ``dr1.mer_segmentation_map``
- ``OTF`` and ``REG``: ``sedm.mer_segmentation_map``

``--idr-field`` and ``--idr-deep-partition`` are not used here because
``mer_segmentation_map`` is not split by those MER catalogue partitions.

Cutout output
-------------

By default, ``compile-segmap`` writes one 10 arcsec x 10 arcsec FITS file per
source row. The output filename pattern is:

.. code-block:: text

   segmap_cutout_row######_obj*_seg*.fits

The output FITS primary HDU contains raw segmentation-label pixels from the MER
segmentation map. Pixels are not converted into a source-only binary mask, and
neighboring object labels inside the cutout are preserved.

Boundary sources use partial cutouts with zero-valued padding. Existing output
files are skipped by default; pass ``--overwrite`` to replace them. Use
``--on-error skip`` to continue past missing files or bad rows, otherwise the
command fails on the first error.

Useful options
--------------

.. code-block:: bash

   euclidkit compile-segmap \
     --input segmentation_maps.fits \
     --output-dir ./segmap_cutouts \
     --size-arcsec 10 \
     --overwrite \
     --on-error skip \
     --no-progress

- ``--size-arcsec`` sets the square cutout size; default is ``10.0``.
- ``--overwrite`` replaces existing cutout FITS files.
- ``--on-error fail|skip`` controls whether missing files or bad rows abort the run.
- ``--progress``/``--no-progress`` controls tile-processing progress output.

Python API
----------

.. code-block:: python

   from euclidkit.core.segmap import compile_segmap_cutouts

   stats = compile_segmap_cutouts(
       input_table="segmentation_maps.fits",
       output_dir="./segmap_cutouts",
       size_arcsec=10.0,
       on_error="skip",
       show_progress=True,
   )

   print(stats.written_rows, stats.failed_rows)

``compile_segmap_cutouts`` returns ``SegmapCutoutStats`` with requested, written,
failed, skipped-existing, and output-file counts.
