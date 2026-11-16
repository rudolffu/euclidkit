Data Releases and Fields
========================

Most commands need to know three things: which archive to query, which data
release, and which survey field. The defaults are the public archive, DR1, and
the WIDE field. This page explains when and how to change them.

Quick reference
---------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - To query...
     - Use
   * - Public DR1, WIDE field
     - No options (the defaults)
   * - Public DR1, DEEP fields
     - ``--field DEEP``
   * - Public Q1
     - ``--release Q1``
   * - The internal consortium archive
     - ``--environment IDR``

Archive and release
-------------------

``--environment`` selects the archive service. ``--release`` selects a data
release on the public service.

- ``PDR`` (default) is the public archive. It serves ``--release DR1`` (default)
  and ``--release Q1``.
- ``IDR`` is the internal archive for Euclid Consortium members. It is a
  separate service, not another name for public DR1. Do not pass ``--release``
  with it.
- ``OTF`` and ``REG`` are development and testing services. Do not pass
  ``--release`` with them.

Passing ``--release`` together with ``IDR``, ``OTF``, or ``REG`` is an error.

In Python, the same choices are arguments of ``EuclidArchive``:

.. code-block:: python

   from euclidkit.core.data_access import EuclidArchive

   dr1 = EuclidArchive()                         # public DR1
   q1 = EuclidArchive(data_release="Q1")         # public Q1
   internal = EuclidArchive(environment="IDR")   # internal archive

WIDE and DEEP fields
--------------------

DR1 and the internal archive store the WIDE and DEEP surveys in separate
tables. Select one with ``--field WIDE`` (default) or ``--field DEEP``.

Each field is further split into two partitions, ``survey`` and ``mode``:

- **WIDE**: you do not choose. ``euclidkit`` queries the ``survey`` table
  first, then queries the ``mode`` table only for sources that were not found.
- **DEEP**: choose with ``--deep-partition``.

  - ``survey`` (default) covers the Euclid Deep Fields: EDF-N, EDF-F, and EDF-S.
  - ``mode`` covers CDFS and COSMOS.
  - ``both`` queries ``survey`` and then ``mode``. If an object is in both,
    the ``survey`` row is kept.

.. code-block:: bash

   # A COSMOS sample
   euclidkit crossmatch -i cosmos_sources.fits -o cosmos_matches.fits \
     --field DEEP --deep-partition mode

In Python, the corresponding keywords are ``idr_field`` and
``idr_deep_partition``. The names are historical. They apply to public DR1 as
well as to the internal archive.

.. code-block:: python

   matches = archive.crossmatch_sources(
       user_table="cosmos_sources.fits",
       idr_field="DEEP",
       idr_deep_partition="mode",
   )

What Q1 does not support
------------------------

Q1 has a single MER catalogue. With ``--release Q1``:

- ``--field DEEP`` is an error.
- ``query-morph`` and ``query-zspe`` are not available.

These commands stop with an error. They never fall back to DR1 data.

Archive tables behind each choice
---------------------------------

You do not need these names to use ``euclidkit``. They are listed so that you
can reproduce a query by hand in ADQL, or cite the tables that you used.

.. list-table::
   :header-rows: 1
   :widths: 22 39 39

   * - Data
     - Public DR1
     - Public Q1
   * - MER catalogue
     - ``dr1.mer_catalogue_wide_survey``, ``_wide_mode``, ``_deep_survey``,
       ``_deep_mode``
     - ``q1.mer_catalogue``
   * - Morphology
     - ``dr1.mer_morphology_*`` with the same four suffixes
     - Not available
   * - Spectra index
     - ``dr1.spectra_source_wide``, ``dr1.spectra_source_deep``
     - ``q1.spectra_source``
   * - Redshift candidates
     - ``dr1.spectro_z_spe_<type>_candidates_wide_survey``, ``_wide``, ``_deep``
     - Not available
   * - Segmentation maps, mosaics, raw frames
     - ``dr1.mer_segmentation_map``, ``dr1.raw_frame``, and other ``dr1.*``
       product tables
     - The same names under ``q1.*``

The internal archive (``IDR``) uses the same DR1 table layout under the
``catalogue`` schema for catalogues, ``dr1.raw_frame`` and
``dr1.mer_segmentation_map`` for those two product tables, and ``sedm`` for
other products.

Every manifest and recovery file that ``euclidkit`` writes records the archive
and the release, so you can check later where a result came from.

Migrating from earlier versions
-------------------------------

- **PDR now defaults to DR1.** Commands that relied on PDR meaning Q1 must add
  ``--release Q1``.
- **Public DR1 work that used ``--environment IDR``** should now use the
  default public archive. Remove the option.
- ``--idr-field`` and ``--idr-deep-partition`` still work. They are aliases
  for ``--field`` and ``--deep-partition``.
- ``crossmatch`` on the public archive writes to exactly the output filename
  that you give. Only ``--environment IDR`` still adds a ``wide_`` or ``deep_``
  prefix to the filename.
