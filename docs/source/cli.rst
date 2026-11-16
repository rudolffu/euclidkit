Command Reference
=================

This page lists every ``euclidkit`` command, what it needs, and its main
options. For the complete and current option list of a command, run it with
``--help``:

.. code-block:: bash

   euclidkit crossmatch --help

All commands
------------

"Login" means that the command queries the archive and needs your
:doc:`credentials <configuration>`. "Datalabs" means that it reads files from
the Euclid data volumes on ESA Datalabs.

.. list-table::
   :header-rows: 1
   :widths: 30 44 10 16

   * - Command
     - Purpose
     - Login
     - Datalabs
   * - :ref:`select-footprint <cmd-select-footprint>`
     - Keep sources inside a Euclid footprint
     - No
     - No
   * - :ref:`crossmatch <cmd-crossmatch>`
     - Match sources to the MER catalogue
     - Yes
     - No
   * - :ref:`upload-table <cmd-upload-table>`
     - Upload a table to your archive user space
     - Yes
     - No
   * - :ref:`query-morph <query-morph>`
     - Add MER morphology by object ID
     - Yes
     - No
   * - :ref:`query-zspe <cmd-query-zspe>`
     - Add SPE redshift candidates by object ID
     - Yes
     - No
   * - :ref:`query-spectra <cmd-query-spectra>`
     - Find the spectra of your sources
     - Yes
     - No
   * - :ref:`compile-spectra <cmd-compile-spectra>`
     - Export spectra to Parquet or FITS
     - Datalink mode only
     - Unless Datalink mode
   * - :ref:`dithers-to-parquet <cmd-dithers-to-parquet>`
     - Export combined and per-dither spectra
     - Yes
     - Yes
   * - :ref:`query-cutana <cmd-query-cutana>`
     - Build a Cutana input catalogue
     - Yes
     - No
   * - :ref:`download-cutana-alternatives <cmd-download-cutana-alternatives>`
     - Download the mosaics that ``query-cutana`` did not select
     - Yes
     - No
   * - :ref:`cutouts <cmd-cutouts>`
     - Write image cutouts directly
     - Yes
     - Yes
   * - :ref:`query-segmap <cmd-query-segmap>`
     - Find segmentation-map files
     - Yes
     - No
   * - :ref:`compile-segmap <cmd-compile-segmap>`
     - Write segmentation-map cutouts
     - No
     - Yes
   * - :ref:`init-config <cmd-init-config>`
     - Write a configuration file template
     - No
     - No
   * - :ref:`diagnostics <cmd-diagnostics>`
     - Check the installation
     - No
     - No

Options shared by archive commands
----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Option
     - Meaning
   * - ``-e``, ``--environment``
     - Archive service: ``PDR`` (default, public), ``IDR``, ``OTF``, or ``REG``
   * - ``--release``
     - Public release: ``DR1`` (default) or ``Q1``. Only with ``PDR``.
   * - ``--field``
     - ``WIDE`` (default) or ``DEEP``. For DR1 and IDR. Alias: ``--idr-field``.
   * - ``--deep-partition``
     - ``survey`` (default), ``mode``, or ``both``. Only with ``--field DEEP``.
       Alias: ``--idr-deep-partition``.
   * - ``-c``, ``--credentials``
     - Path of a credentials file
   * - ``-v``, ``--verbose``
     - Print more detail

:doc:`environments` explains these choices. Not every command accepts every
option. For example, only the MER-based commands accept ``--deep-partition``.

Tables can be read and written as FITS (``.fits``), CSV (``.csv``), VOTable
(``.vot``, ``.xml``), ECSV (``.ecsv``), and Parquet (``.parquet``). The format
follows the file extension.

Catalogue commands
------------------

.. _cmd-select-footprint:

select-footprint
~~~~~~~~~~~~~~~~

Keep the rows of a table that lie inside a packaged Euclid DR1 footprint.
Guide: :ref:`select-footprint-guide`.

.. code-block:: bash

   euclidkit select-footprint -i sources.fits -o in_footprint.fits --survey wide

- ``-i``, ``--input`` and ``-o``, ``--output``: input and output tables
  (required)
- ``-s``, ``--survey``: ``wide`` (default), ``deep``, ``both``, ``cgv``, or
  ``cgv_input``
- ``-d``, ``--data-release``: ``dr1``, the only release available
- ``--ra-col``, ``--dec-col``: coordinate columns in degrees (default ``ra``,
  ``dec``)
- ``--moc-path``: use your own MOC file instead

.. _cmd-crossmatch:

crossmatch
~~~~~~~~~~

Match a table to the MER catalogue by position or by object ID.
Guide: :ref:`crossmatch-guide`.

.. code-block:: bash

   euclidkit crossmatch -i sources.fits -o matches.fits --radius 1.0

- ``-i``, ``--input``: local input table. Give this or ``--user-table-name``.
- ``--user-table-name``: name of a table already in your archive user space
- ``-o``, ``--output``: output table (required)
- ``-r``, ``--radius``: match radius in arcsec (default 1.0)
- ``--ra-col``, ``--dec-col``: coordinate columns in degrees (default ``ra``,
  ``dec``)
- ``--match-mode``: ``auto`` (default), ``object-id``, or ``spatial``
- ``--max-sources``: process only the first N rows
- ``--full-async``: use asynchronous archive jobs
- ``--async-chunk-size``: rows per asynchronous job (default 500000)
- ``--drop-empty-columns``: remove output columns that are entirely missing

.. _cmd-upload-table:

upload-table
~~~~~~~~~~~~

Upload a local table to your archive user space.
Guide: :ref:`crossmatch-guide`.

.. code-block:: bash

   euclidkit upload-table -i sources.fits -t my_sources

- ``-i``, ``--input``: local table (required)
- ``-t``, ``--table-name``: name of the table in the archive (required)
- ``-d``, ``--description``: description stored with the table
- ``-f``, ``--format``: input format, if the file extension is not enough
- ``--overwrite``: replace an existing table of the same name

.. _query-morph:

query-morph
~~~~~~~~~~~

Add MER morphology measurements to a table with ``object_id``. DR1 and IDR
only. Guide: :ref:`morphology-guide`.

.. code-block:: bash

   euclidkit query-morph -i sources.fits -o morphology.fits

- ``-i``, ``--input`` and ``-o``, ``--output``: input and output tables
  (required)
- ``--full-async``, ``--async-chunk-size``: as for ``crossmatch``

Without ``--full-async``, the query runs in batches of 1,000 unique IDs.

.. _cmd-query-zspe:

query-zspe
~~~~~~~~~~

Add SPE redshift candidates to a table with ``object_id``. DR1 and IDR only.
Guide: :ref:`redshift-guide`.

.. code-block:: bash

   euclidkit query-zspe -x matches.fits -o redshifts.fits --object-type qso

- ``-x``, ``--crossmatch`` and ``-o``, ``--output``: input and output tables
  (required)
- ``--object-type``: ``qso`` (default) or ``galaxy``
- ``--full-async``, ``--async-chunk-size``: as for ``crossmatch``

Spectra commands
----------------

Guide for all three commands: :doc:`spectra_compilation`.

.. _cmd-query-spectra:

query-spectra
~~~~~~~~~~~~~

Find the spectra of the sources in a table of object IDs or positions.

.. code-block:: bash

   euclidkit query-spectra -x matches.fits -o spectra_sources.fits

- ``-x``, ``--crossmatch``: input table (required). It does not have to be a
  ``crossmatch`` output.
- ``-o``, ``--output``: output table (required)
- ``--match-mode``: ``auto`` (default), ``object-id``, or ``spatial``
- ``-r``, ``--radius``: match radius in arcsec for spatial mode (default 1.0)
- ``--ra-col``, ``--dec-col``: coordinate columns for spatial mode
- ``--combine-output``: also write the spectra into this FITS file
- ``--max-spectra``: limit the number of spectra in ``--combine-output``

.. _cmd-compile-spectra:

compile-spectra
~~~~~~~~~~~~~~~

Export the spectra listed by ``query-spectra``.

.. code-block:: bash

   euclidkit compile-spectra -s spectra_sources.fits -o ./output -L RGS

- ``-s``, ``--spectra-table`` and ``-o``, ``--output-dir``: input table and
  output directory (required)
- ``--prefix``: start of the output filenames (default ``compiled_spectra``)
- ``-L``, ``--lambda-range``: ``RGS`` (default), ``BGS``, or ``BOTH``
- ``--output-format``: ``parquet`` (default) or ``fits``, for the local modes
- ``--use-datalink``: download through the archive. Always writes FITS.
- ``--chunk-size``: spectra per Parquet file (default 2000)
- ``--max-extensions``: spectra per FITS file
- ``--workers``: parallel workers. Default: up to 8 for Parquet, 1 for FITS
  and Datalink.
- ``--limit``: process only the first N rows
- ``--on-error``: ``fail`` (default) or ``skip``, for Parquet
- ``--overwrite``: replace existing output with the same prefix
- ``--schema``: archive schema for Datalink. The default follows the release.
- ``--retrieval-type``: older Datalink grism selector. Prefer ``-L``.

.. _cmd-dithers-to-parquet:

dithers-to-parquet
~~~~~~~~~~~~~~~~~~

Export the combined spectrum and each single-dither spectrum to Parquet.

.. code-block:: bash

   euclidkit dithers-to-parquet --catalog-table spectra_sources.fits \
     --output-prefix ./output/raw_sir --lambda-range RGS

- ``--catalog-table``: ``query-spectra`` output (required)
- ``--output-prefix``: path and start of the output filenames (required)
- ``--lambda-range``: ``RGS`` or ``BGS``. Without it, both are exported.
- ``--dithers-only``: do not write the combined spectra
- ``--chunk-size``, ``--workers``, ``--on-error``, ``--overwrite``: as for
  ``compile-spectra``

Image commands
--------------

Guide for the three image commands: :doc:`cutouts`.

.. _cmd-query-cutana:

query-cutana
~~~~~~~~~~~~

Build a Cutana input catalogue from object IDs or positions.

.. code-block:: bash

   euclidkit query-cutana -s sources.fits -o cutana_input.csv --instrument VIS

- ``-s``, ``--sources`` and ``-o``, ``--output``: input table and output CSV
  (required)
- ``--instrument``: ``VIS`` (default) or ``NISP``
- ``--nisp-filters``: comma-separated list, for example ``NIR_Y,NIR_H``
- ``--cutout-size``: ``arcsec`` (default) or ``pixel``
- ``--cutout-size-value``: cutout diameter (default 15)
- ``--keep-noncutana-cols``: keep your other input columns

.. _cmd-download-cutana-alternatives:

download-cutana-alternatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the candidate mosaics that ``query-cutana`` did not select.

.. code-block:: bash

   euclidkit download-cutana-alternatives \
     --candidates cutana_input.alternatives.csv --output-dir alternatives

- ``--candidates``: the ``.alternatives.csv`` file (required)
- ``--output-dir``: download directory (required)
- ``--schema``: archive product schema. The default follows the release.

.. _cmd-cutouts:

cutouts
~~~~~~~

Write image cutouts directly from the Datalabs volumes. ``euclidkit-cutouts``
is the same command.

.. code-block:: bash

   euclidkit cutouts --object-id 123456789 --output ./cutouts --size 10

- ``--object-id``, ``--ra`` with ``--dec``, ``-s``/``--source``, or
  ``-f``/``--sources-file``: the source or sources
- ``-o``, ``--output``: output directory (required)
- ``--size``: cutout size in arcsec (default 10)
- ``--instrument``: ``VIS`` (default) or ``NISP``
- ``--nisp-filters``: comma-separated list of ``NIR_Y``, ``NIR_J``, ``NIR_H``
- ``--product-type``: ``mosaic`` (default), ``calibrated_frame``, or
  ``stacked_frame``

Segmentation-map commands
-------------------------

Guide for both commands: :doc:`segmentation_maps`.

.. _cmd-query-segmap:

query-segmap
~~~~~~~~~~~~

Find the segmentation-map file for each row of a ``crossmatch`` output.

.. code-block:: bash

   euclidkit query-segmap -i matches.fits -o segmentation_maps.fits

- ``-i``, ``--input`` and ``-o``, ``--output``: input and output tables
  (required)

.. _cmd-compile-segmap:

compile-segmap
~~~~~~~~~~~~~~

Write one segmentation-map cutout per row of a ``query-segmap`` output.

.. code-block:: bash

   euclidkit compile-segmap -i segmentation_maps.fits -o ./segmap_cutouts

- ``-i``, ``--input`` and ``-o``, ``--output-dir``: input table and output
  directory (required)
- ``--size-arcsec``: side of the square cutout in arcsec (default 10)
- ``--overwrite``: replace existing cutouts
- ``--on-error``: ``fail`` (default) or ``skip``

Setup commands
--------------

.. _cmd-init-config:

init-config
~~~~~~~~~~~

Write a configuration file template. See :doc:`configuration`.

.. code-block:: bash

   euclidkit init-config

- ``-o``, ``--output``: where to write the file (default
  ``~/.euclidkit/euclidkit_config.yaml``)
- ``-t``, ``--template``: ``basic`` (default), ``advanced``, or ``pipeline``

.. _cmd-diagnostics:

diagnostics
~~~~~~~~~~~

Check the installation. See :doc:`installation` for what ``--check-data``
covers.

.. code-block:: bash

   euclidkit diagnostics --check-deps

- ``--check-deps``: check that the dependencies import
- ``--check-data``: check for the Q1 and ERO data volumes, and for a
  credentials file in the configured location
- ``--check-desi``: check for the optional ``sparcl-client`` package

Without options, all three checks run.
