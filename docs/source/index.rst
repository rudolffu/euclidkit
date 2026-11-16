euclidkit
=========

.. image:: https://img.shields.io/badge/GitHub-rudolffu%2Feuclidkit-black?logo=github
   :target: https://github.com/rudolffu/euclidkit
   :alt: GitHub repository

.. image:: https://img.shields.io/pypi/v/euclidkit.svg
   :target: https://pypi.org/project/euclidkit/
   :alt: PyPI version

``euclidkit`` takes you from a list of targets to Euclid data for those
targets: MER photometry, morphology, spectroscopic redshifts, NISP slitless
spectra, image cutouts, and segmentation maps.

You bring a table of sky positions or Euclid object IDs. ``euclidkit`` writes
the archive queries, handles uploads and large asynchronous jobs, and saves the
results as FITS, CSV, or Parquet tables. You can use it from the command line
or from Python.

It works with the public Euclid releases (DR1 by default, and Q1) and with the
internal consortium archive. Catalogue queries run from any computer. Bulk
access to spectra and images reads the Euclid data volumes on
`ESA Datalabs <https://datalabs.esa.int/>`_.

Where to start
--------------

New to ``euclidkit``? Read these three pages in order:

1. :doc:`installation`
2. :doc:`configuration`, to set up your archive credentials
3. :doc:`quickstart`, a complete worked example

Then go to the page for your task:

.. list-table::
   :header-rows: 1
   :widths: 55 45

   * - I want to...
     - Read
   * - Choose between DR1 and Q1, or between the WIDE and DEEP fields
     - :doc:`environments`
   * - Check which of my sources are inside the Euclid footprint
     - :ref:`select-footprint-guide`
   * - Get Euclid photometry and quality flags for my sources
     - :ref:`crossmatch-guide`
   * - Get morphology or spectroscopic redshifts
     - :ref:`morphology-guide` and :ref:`redshift-guide`
   * - Get NISP spectra, combined or per dither
     - :doc:`spectra_compilation`
   * - Make image cutouts
     - :doc:`cutouts`
   * - Make segmentation-map cutouts
     - :doc:`segmentation_maps`
   * - Look up a command, an option, or a Python function
     - :doc:`cli` and :doc:`api`
   * - Understand a Euclid term such as MER, SIR, or RGS
     - :doc:`glossary`

Getting help
------------

Report bugs and ask questions on the
`GitHub issue tracker <https://github.com/rudolffu/euclidkit/issues>`_.
Release history is in the
`changelog <https://github.com/rudolffu/euclidkit/blob/main/CHANGELOG.md>`_.

``euclidkit`` is distributed under the BSD 3-Clause License.

About this documentation
------------------------

This documentation is revised with
`nicedocs <https://github.com/rudolffu/nicedocs>`_.

.. toctree::
   :maxdepth: 2
   :caption: Get started
   :hidden:

   installation
   configuration
   quickstart
   environments

.. toctree::
   :maxdepth: 2
   :caption: Guides
   :hidden:

   catalogues
   spectra_compilation
   cutouts
   segmentation_maps

.. toctree::
   :maxdepth: 2
   :caption: Reference
   :hidden:

   cli
   api
   glossary
