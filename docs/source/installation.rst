Installation
============

Requirements
------------

- Python 3.11 or later.
- An ESA Cosmos account, for archive queries. See :doc:`configuration`.
- An `ESA Datalabs <https://datalabs.esa.int/>`_ session, only for the commands
  that read spectra and images directly from the Euclid data volumes. The table
  in :doc:`cli` lists which commands need it.

Install from PyPI
-----------------

.. code-block:: bash

   pip install euclidkit

This installs all dependencies, including ``astroquery`` for archive access
and ``pyarrow`` for Parquet output.

Check the installation
----------------------

.. code-block:: bash

   euclidkit --version
   euclidkit --help

``euclidkit --help`` lists every command. Each command has its own help, for
example ``euclidkit crossmatch --help``.

To check that all dependencies import correctly:

.. code-block:: bash

   euclidkit diagnostics --check-deps

To check that ``euclidkit`` can find your credentials file:

.. code-block:: bash

   euclidkit diagnostics --check-data

.. note::

   ``--check-data`` also looks for the Q1 and ERO data volumes
   (``/data/euclid_q1/`` and ``/data/euc_ero_data_01/``). That check fails
   outside ESA Datalabs, which is expected. It does not test DR1 file access.

Install for development
-----------------------

.. code-block:: bash

   git clone https://github.com/rudolffu/euclidkit.git
   cd euclidkit
   pip install -e ".[dev]"
   pytest -q

The tests mock the archive and the data volumes, so they run without
credentials and outside Datalabs.

Next step
---------

Set up your archive credentials: :doc:`configuration`.
