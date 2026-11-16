Quickstart
==========

This example starts from a list of quasar candidates with sky positions and
ends with their Euclid photometry and NISP spectra. It uses the public DR1
release and the WIDE field, which are the defaults.

Before you start, :doc:`install euclidkit <installation>` and
:doc:`set up your credentials <configuration>`.

Steps 1 to 3 run on any computer. Step 4 reads the Euclid data volumes, so run
it in an ESA Datalabs session.

Prepare the input table
-----------------------

Create ``my_sources.csv`` with one row per target. Coordinates are in degrees
(ICRS):

.. code-block:: text

   name,ra,dec
   target_a,51.8975,-35.1766
   target_b,152.1781,68.8275

These two example positions lie inside the DR1 WIDE footprint. Replace them
with your own targets. FITS, VOTable, ECSV, and Parquet tables also work. If your columns have other
names, pass ``--ra-col`` and ``--dec-col``.

Step 1: Keep the sources inside the Euclid footprint
----------------------------------------------------

This step is optional, but it saves archive time for large lists. It uses a
footprint map shipped with the package and needs no login:

.. code-block:: bash

   euclidkit select-footprint -i my_sources.csv -o in_footprint.csv --survey wide

Step 2: Crossmatch with the MER catalogue
-----------------------------------------

.. code-block:: bash

   euclidkit crossmatch -i in_footprint.csv -o matches.fits --radius 1.0

``matches.fits`` contains your input columns plus, for each matched Euclid
object:

- ``object_id``, the Euclid identifier that all later steps use
- ``mer_ra``, ``mer_dec``, and ``separation_arcsec``
- VIS and NISP fluxes with errors, and fluxes in the ground-based bands
- quality flags, ``point_like_prob``, and ``extended_prob``
- ``segmentation_map_id``, which the segmentation-map workflow needs

See :ref:`crossmatch-guide` for DEEP fields, object-ID matching, and tables
with millions of rows.

Step 3: Find the spectra
------------------------

.. code-block:: bash

   euclidkit query-spectra -x matches.fits -o spectra_sources.fits

``spectra_sources.fits`` has one row per available spectrum. Each row records
which FITS file holds the spectrum and which extension to read. Objects with
no spectrum do not appear.

Step 4: Export the spectra (on Datalabs)
----------------------------------------

.. code-block:: bash

   euclidkit compile-spectra -s spectra_sources.fits -o ./spectra \
     --prefix qso -L RGS

This writes the red-grism (RGS) spectra to ``./spectra/qso_part001.parquet``
and following parts, with a manifest file that records what was exported.

If you are not on Datalabs, add ``--use-datalink`` to download the spectra
through the archive instead. That mode writes FITS files. See
:doc:`spectra_compilation`.

Step 5: Read the spectra
------------------------

Each Parquet row is one spectrum. The array columns hold the full spectrum:

.. code-block:: python

   import matplotlib.pyplot as plt
   import pandas as pd

   spectra = pd.read_parquet("spectra/qso_part001.parquet")
   row = spectra.iloc[0]

   plt.plot(row["wavelength"], row["flux"])
   plt.title(f"Euclid object {row['object_id']}")
   plt.show()

The same workflow in Python
---------------------------

.. code-block:: python

   from euclidkit.core.data_access import EuclidArchive
   from euclidkit.core.spectra_parquet import spectra_to_parquet

   archive = EuclidArchive()          # public DR1
   archive.login()                    # reads ~/.euclidkit/.cred.txt

   matches = archive.crossmatch_sources(
       user_table="in_footprint.csv",
       radius=1.0,                    # arcsec
       output_file="matches.fits",
   )
   spectra_sources = archive.query_spectra_sources(
       crossmatch_table=matches,
       output_file="spectra_sources.fits",
   )
   archive.logout()

   # On Datalabs
   stats = spectra_to_parquet(
       catalog_table="spectra_sources.fits",
       output_prefix="./spectra/qso",
       lambda_range="RGS",
   )
   print(stats.exported_rows, "spectra written")

Where to go next
----------------

- Add morphology or redshifts to your table: :doc:`catalogues`
- Work in the DEEP fields or with Q1: :doc:`environments`
- Make image cutouts: :doc:`cutouts`
