Image Cutouts
=============

``euclidkit`` offers two ways to get image cutouts of your sources:

- **For many sources**, prepare an input catalogue for
  `Cutana <https://github.com/esa/Cutana>`_, ESA's tool for making cutouts in
  bulk from large collections of FITS tiles. ``query-cutana`` finds the right
  mosaic for each source and band. Cutana then extracts the cutouts.
- **For a few sources**, ``euclidkit cutouts`` writes FITS cutouts directly.

Prepare a Cutana catalogue
--------------------------

``query-cutana`` takes a table with Euclid ``object_id`` values or with
``ra``/``dec`` positions in degrees.

VIS cutouts, 15 arcsec across:

.. code-block:: bash

   euclidkit query-cutana -s my_sources.fits -o cutana_input.csv \
     --instrument VIS --cutout-size arcsec --cutout-size-value 15

NISP cutouts in two bands, for DEEP-field sources:

.. code-block:: bash

   euclidkit query-cutana -s my_sources.fits -o cutana_nisp.csv \
     --instrument NISP --nisp-filters NIR_Y,NIR_H \
     --field DEEP --deep-partition both

Options:

- ``--cutout-size`` is the unit of the cutout diameter, ``arcsec`` (default)
  or ``pixel``. ``--cutout-size-value`` is the diameter (default 15). It must
  be positive, at most 3600 arcsec or 10000 pixels. Pixel diameters must be
  integers.
- ``--nisp-filters`` defaults to ``NIR_Y,NIR_J,NIR_H``. The order that you give
  is kept.
- ``--keep-noncutana-cols`` copies your other input columns into the
  catalogue. By default they are dropped.

Output files
~~~~~~~~~~~~

For ``-o cutana_input.csv``, the command writes four kinds of file:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - File
     - Content
   * - ``cutana_input.csv``
     - The Cutana catalogue: one row per source, with the mosaic file for each
       band
   * - ``cutana_input.coverage.csv``
     - For each source, which bands were found and which are missing
   * - ``cutana_input.alternatives.csv``
     - For sources covered by more than one mosaic, every candidate and its
       rank
   * - ``cutana_input.manifest.json``
     - The band groups, the source counts, and the paths of the files above

How the mosaic is chosen
~~~~~~~~~~~~~~~~~~~~~~~~

Euclid mosaics are tiles that overlap, so a source can lie in several. For
each source and band, ``euclidkit``:

1. Keeps the mosaics whose archive footprint contains the source.
2. Chooses the mosaic whose image centre is nearest to the source.
3. Breaks ties by file path, so that the choice is the same on every run.

.. note::

   This rule makes the result reproducible. It does not select the deepest or
   the best-quality image. If that matters for your science, inspect
   ``.alternatives.csv``.

Keep Cutana's own FITS checks switched on. The archive footprint is not a
substitute for checking the image WCS. The catalogue contains Datalabs file
paths, so run Cutana where those paths are readable.

Sources with missing bands
~~~~~~~~~~~~~~~~~~~~~~~~~~

A source that lacks one of the requested bands is kept. When sources differ in
which bands they have, ``euclidkit`` also writes one catalogue for each set of
bands, named ``cutana_input.bands-*.csv``. The manifest lists them.

In that case, **run Cutana once for each group**, and configure Cutana's
channels for the bands of that group. A single multi-band Cutana configuration
can reject sources that lack one of its channels.

A source with no mosaic in any band is listed in ``.coverage.csv`` and cannot
be extracted.

Duplicate sources
~~~~~~~~~~~~~~~~~

Rows with the same source and position are merged into one. The same object ID
at different positions gives separate rows with distinct IDs.

Download the alternative mosaics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To compare the chosen mosaic with the other candidates, download them later.
The source query is not repeated:

.. code-block:: bash

   euclidkit download-cutana-alternatives \
     --candidates cutana_input.alternatives.csv --output-dir alternatives

Use the same ``--environment`` and ``--release`` as for ``query-cutana``. Each
file is downloaded once. Files that already exist and are not empty are
skipped, so you can rerun the command after an interruption. To download only
some candidates, filter the CSV file first.

Make cutouts directly
---------------------

``euclidkit cutouts`` writes one FITS cutout per source and band. It reads the
image files from the Datalabs volumes, so run it on ESA Datalabs. The
``euclidkit-cutouts`` command is the same tool.

By object ID:

.. code-block:: bash

   euclidkit cutouts --object-id 123456789 --output ./cutouts --size 10

By position, in degrees:

.. code-block:: bash

   euclidkit cutouts --ra 150.0 --dec 2.5 --output ./cutouts \
     --instrument NISP --nisp-filters NIR_Y,NIR_J,NIR_H

For a list of sources:

.. code-block:: bash

   euclidkit cutouts --sources-file sources.csv --output ./cutouts

Options:

- ``--size`` is the cutout size in arcsec (default 10).
- ``--instrument`` is ``VIS`` (default) or ``NISP``.
- ``--product-type`` is ``mosaic`` (default), ``calibrated_frame``, or
  ``stacked_frame``.
- ``--sources-file`` accepts a CSV file, a JSON file, or a text file (``.txt``
  or ``.dat``) with one object ID or one ``ra dec`` pair per line.

Each output file is named after the source image, followed by
``-CUTOUT_<ra>_<dec>.fits``.
