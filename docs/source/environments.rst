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

- ``WIDE``: queries IDR WIDE MER catalogue.
- ``DEEP``: queries IDR DEEP MER catalogue.

For commands where ``--idr-field`` is available, the default is ``WIDE``.

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

Cutana query in IDR DEEP:

.. code-block:: bash

   euclidkit query-cutana \
     --sources my_sources.fits \
     --output cutana_deep.csv \
     --environment IDR \
     --idr-field DEEP

