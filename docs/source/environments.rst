Archive Environments and IDR Fields
===================================

euclidkit supports multiple Euclid archive environments via ``--environment``
in CLI commands (or ``EuclidArchive(environment=...)`` in Python).

Supported environments
----------------------

- ``PDR``: Public Data Release archive.
- ``IDR``: Internal Data Release archive (consortium access).
- ``OTF``: On-the-fly archive environment.
- ``REG``: Regression/testing archive environment.

IDR field selection
-------------------

When using ``IDR``, select the field with ``--idr-field``:

- ``WIDE``: queries IDR WIDE MER catalogues. MER crossmatch queries
  ``catalogue.mer_catalogue_wide_survey`` first, then queries
  ``catalogue.mer_catalogue_wide_mode`` only for sources not matched in the
  survey table.
- ``DEEP``: queries IDR DEEP MER partitions. MER crossmatch defaults to
  ``catalogue.mer_catalogue_deep_survey``.

For commands where ``--idr-field`` is available, the default is ``WIDE``.

IDR DEEP MER partition selection
--------------------------------

MER-based commands that support ``--idr-field DEEP`` also support
``--idr-deep-partition``:

- ``survey`` (default): query ``catalogue.mer_catalogue_deep_survey`` for EDFN,
  EDFF, and EDFS.
- ``mode``: query ``catalogue.mer_catalogue_deep_mode`` for CDFS and COSMOS.
- ``both``: query ``deep_survey`` first, then ``deep_mode``. If an unexpected
  duplicate ``object_id`` appears in both partitions, the survey row is kept.

This selector applies to MER workflows such as ``crossmatch`` and
``query-cutana``. It does not change spectra-source or SPE redshift table
selection.

Examples
--------

Crossmatch in IDR WIDE:

.. code-block:: bash

   euclidkit crossmatch \
     --input my_sources.fits \
     --output xmatch_wide.fits \
     --environment IDR \
     --idr-field WIDE

Crossmatch in IDR DEEP:

.. code-block:: bash

   euclidkit crossmatch \
     --input my_sources.fits \
     --output xmatch_deep.fits \
     --environment IDR \
     --idr-field DEEP

Crossmatch in the IDR DEEP mode partition:

.. code-block:: bash

   euclidkit crossmatch \
     --input my_sources.fits \
     --output xmatch_deep_mode.fits \
     --environment IDR \
     --idr-field DEEP \
     --idr-deep-partition mode

Cutana query in IDR DEEP:

.. code-block:: bash

   euclidkit query-cutana \
     --sources my_sources.fits \
     --output cutana_deep.csv \
     --environment IDR \
     --idr-field DEEP \
     --idr-deep-partition both
