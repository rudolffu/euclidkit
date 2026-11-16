Spectra
=======

This guide shows how to collect the NISP slitless spectra of your sources. It
takes two steps:

1. ``query-spectra`` asks the archive which of your sources have spectra and
   where the files are.
2. ``compile-spectra`` reads those spectra and writes them into a few large
   files.

.. code-block:: bash

   euclidkit query-spectra -x my_sources.fits -o spectra_sources.fits
   euclidkit compile-spectra -s spectra_sources.fits -o ./output \
     --prefix raw_spectra -L RGS

Choose how to read the spectra
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 24 22 22 32

   * - Mode
     - Where it runs
     - Output
     - Use it when
   * - Local Parquet (default)
     - ESA Datalabs
     - Parquet tables
     - You want many spectra quickly for analysis or machine learning
   * - Local FITS (``--output-format fits``)
     - ESA Datalabs
     - Multi-extension FITS
     - Your tools expect one FITS extension per spectrum
   * - Datalink (``--use-datalink``)
     - Any computer
     - Multi-extension FITS
     - You are not on Datalabs. It is slower, because each spectrum is
       downloaded from the archive.

Step 1: Find the spectra
------------------------

``query-spectra`` accepts any table with Euclid object IDs or with positions.
A ``crossmatch`` output works, but you do not need one.

By object ID:

.. code-block:: bash

   euclidkit query-spectra -x my_ids.fits -o spectra_sources.fits \
     --match-mode object-id

By position, within 1 arcsec:

.. code-block:: bash

   euclidkit query-spectra -x my_positions.fits -o spectra_sources.fits \
     --match-mode spatial --radius 1.0

With the default ``--match-mode auto``, ``euclidkit`` uses object IDs if the
table has an ``object_id`` or ``object_id_euclid`` column, and positions
otherwise.

Things to know:

- **Object-ID mode** reads the ``object_id`` column. If your IDs are in a
  column called ``source_id``, rename it to ``object_id`` first.
- **Spatial mode** keeps only the nearest spectrum for each input row. It
  finds coordinate columns named ``ra``/``dec``, and also common variants such
  as ``RA``/``DEC``, ``right_ascension``/``declination``,
  ``ra_deg``/``dec_deg``, and ``RAJ2000``/``DEJ2000``. Use ``--ra-col`` and
  ``--dec-col`` for anything else.
- **WIDE or DEEP**: for DR1, ``--field WIDE`` (default) and ``--field DEEP``
  query different spectra tables. See :doc:`environments`.
- Spectra that have no file on the Datalabs volumes are left out.

The output has one row per spectrum:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Columns
     - Content
   * - ``object_id``, ``source_id``
     - Euclid identifier of the object
   * - ``ra_obj``, ``dec_obj``
     - Position of the spectrum's source, in degrees
   * - ``datalabs_path``, ``file_name``, ``hdu_index``
     - The FITS file on the Datalabs volumes, and the extension that holds
       this spectrum
   * - ``input_row_id``, ``input_ra``, ``input_dec``, ``separation_arcsec``
     - Spatial mode only: the input row that matched, and the match distance

To also write the spectra into one FITS file in the same call, add
``--combine-output combined.fits``. ``--max-spectra N`` limits how many are
included. This needs the Datalabs volumes.

Step 2: Export to Parquet (default)
-----------------------------------

Run this on ESA Datalabs:

.. code-block:: bash

   euclidkit compile-spectra -s spectra_sources.fits -o ./output \
     --prefix raw_spectra --chunk-size 2000 --workers 8 -L RGS

This writes:

- ``raw_spectra_part001.parquet``, ``raw_spectra_part002.parquet``, and so on,
  with ``--chunk-size`` spectra per file (default 2000)
- ``raw_spectra_manifest.json``, which records the settings, the release, the
  row counts, and the list of output files
- ``raw_spectra_failures.jsonl``, only if you use ``--on-error skip`` and some
  rows fail

Each Parquet row is one spectrum, with these columns:

- ``object_id``, ``source_id``, ``ra``, ``dec``
- ``wavelength`` and ``flux`` arrays. ``flux`` is a copy of the ``signal``
  column.
- every other array column of the spectrum's FITS extension, unchanged

The exporter writes the archive values as they are. It does not compute an
error array, an inverse variance, or a validity mask. For the meaning and
units of each column, see the SIR combined-spectra product description in the
`Euclid DR1 Data Product Description <https://euclid.esac.esa.int/dr/dr1/dpdd/>`_.

Useful options:

- ``--workers N`` sets the number of parallel readers. The default is the
  number of CPU cores, up to 8.
- ``--limit N`` processes only the first ``N`` rows, for a quick test.
- ``--on-error skip`` continues past unreadable spectra. The default,
  ``fail``, stops at the first one.
- ``--overwrite`` is required if output files with the same prefix exist.

Choose the grism
----------------

``-L`` (``--lambda-range``) selects the grism:

- ``-L RGS``: red grism (default)
- ``-L BGS``: blue grism
- ``-L BOTH``: both, written as two separate sets of files,
  ``<prefix>_rgs_*`` and ``<prefix>_bgs_*``

In the local modes, the grism of each spectrum is read from the ``LRANGE``
keyword in the FITS primary header.

Export to FITS on Datalabs
--------------------------

.. code-block:: bash

   euclidkit compile-spectra -s spectra_sources.fits -o ./output \
     --prefix compiled_spectra --output-format fits --max-extensions 1000

This writes ``compiled_spectra_chunk_001.fits`` and so on. Each file holds up
to ``--max-extensions`` spectra, one per extension, and each extension header
carries the source ID in ``SOURC_ID``. A metadata table lists which spectrum
is in which file and extension.

If chunk files already exist, the run resumes after the last compiled
spectrum. Add ``--overwrite`` to start again. This mode uses one worker unless
you set ``--workers``.

Download through Datalink
-------------------------

Outside Datalabs, add ``--use-datalink``. ``euclidkit`` then requests each
spectrum from the archive by source ID. This needs your credentials.

.. code-block:: bash

   euclidkit compile-spectra -s spectra_sources.fits -o ./output \
     --prefix compiled_dl --use-datalink -L BOTH

The output is always FITS:

- ``compiled_dl_rgs_chunk_###.fits``
- ``compiled_dl_bgs_chunk_###.fits``

The archive schema follows the release, for example ``dr1`` for public DR1.
Override it with ``--schema`` if you need to. The older ``--retrieval-type``
option still works, but ``-L`` is the preferred way to select the grism.

.. _dithers-guide:

Export the individual dither spectra
------------------------------------

The archive's combined spectrum is a stack of several dithered exposures.
``dithers-to-parquet`` exports the combined spectrum together with each
single-dither spectrum. Use it to study variability between exposures, or to
check whether a feature appears in every dither.

Run this on ESA Datalabs. It also needs archive credentials, to look up the
observation time of each exposure.

.. code-block:: bash

   euclidkit dithers-to-parquet --catalog-table spectra_sources.fits \
     --output-prefix ./output/raw_sir --lambda-range RGS --workers 8

This writes:

- ``raw_sir_combined_part001.parquet``, and so on: the combined spectra
- ``raw_sir_dithers_part001.parquet``, and so on: one row per dither spectrum
- ``raw_sir_manifest.json``

Add ``--dithers-only`` to skip the combined spectra, for example if you have
already exported them with ``compile-spectra``.

Each dither row includes these columns from the exposure:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Column
     - Content
   * - ``dither_id``, ``ptgid``
     - Dither number and pointing ID
   * - ``gwa_pos``, ``gwa_tilt``
     - Grism wheel position and tilt
   * - ``obs_time_mjd``, ``obs_time_utc``
     - Observation time of the exposure
   * - ``pa``
     - Position angle of the exposure

The last three come from the archive's ``raw_frame`` table. ``euclidkit``
matches each dither to its raw frame by pointing ID and grism wheel position.
The manifest reports how many dither rows found a match.

Troubleshooting
---------------

``compile-spectra`` reports missing files
   The local modes need the Euclid data volumes. Run them on ESA Datalabs, or
   use ``--use-datalink``.

``compile-spectra`` stops because output files exist
   Choose another ``--prefix``, or add ``--overwrite``.

The export is empty for ``-L BGS``
   The table may contain only red-grism spectra. Check the
   ``lambda_skipped_rows`` count in the manifest.
