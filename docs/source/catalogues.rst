Catalogues: Photometry, Morphology, and Redshifts
=================================================

This guide covers the commands that add Euclid catalogue measurements to your
own source table. All of them run on any computer with archive credentials.
None needs ESA Datalabs.

The examples use public DR1 and the WIDE field. To change either, see
:doc:`environments`.

.. _select-footprint-guide:

Select sources inside the Euclid footprint
------------------------------------------

``select-footprint`` keeps the rows of your table that fall inside a Euclid
DR1 footprint. It uses coverage maps (MOCs) shipped with the package, so it
needs no login and no network.

.. code-block:: bash

   euclidkit select-footprint -i my_sources.fits -o in_wide.fits --survey wide

.. list-table::
   :header-rows: 1
   :widths: 20 55 25

   * - ``--survey``
     - Footprint
     - Area (deg²)
   * - ``wide``
     - DR1 WIDE catalogue footprint (default)
     - 2014.6
   * - ``deep``
     - DR1 DEEP catalogue footprint
     - 77.8
   * - ``both``
     - Union of WIDE and DEEP
     - 2058.3
   * - ``cgv``
     - CGV input footprint map
     - 2108.5

Coordinates are read from ``ra`` and ``dec`` in degrees. Use ``--ra-col`` and
``--dec-col`` for other column names, or ``--moc-path`` to supply your own MOC
file.

The maps have HEALPix order 13, which is a resolution of about 26 arcsec. A
source close to the footprint edge can pass this filter and still have no
catalogue match.

.. _crossmatch-guide:

Crossmatch with the MER catalogue
---------------------------------

``crossmatch`` finds the Euclid MER object for each of your sources and
returns its photometry and flags.

.. code-block:: bash

   euclidkit crossmatch -i my_sources.fits -o matches.fits --radius 1.0

Input
~~~~~

The input table can be FITS, CSV, VOTable, ECSV, or Parquet. It needs one of
the following:

- **Positions**: ``ra`` and ``dec`` in degrees (ICRS). Use ``--ra-col`` and
  ``--dec-col`` for other names. ``--radius`` is the match radius in arcsec
  (default 1.0).
- **Euclid object IDs**: an ``object_id`` column. ``object_id_euclid`` and
  ``source_id`` are also recognized.

With the default ``--match-mode auto``, ``euclidkit`` matches by object ID if
an ID column is present and by position otherwise. Force one with
``--match-mode object-id`` or ``--match-mode spatial``.

.. tip::

   If you already have Euclid object IDs, match by ID. It is faster than a
   positional match and cannot pick a wrong neighbour.

Output
~~~~~~

The output contains your input columns and these MER columns for each match:

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Columns
     - Content
   * - ``object_id``
     - Euclid MER identifier. Named ``mer_object_id`` if your input already
       has an ``object_id`` column.
   * - ``mer_ra``, ``mer_dec``, ``separation_arcsec``
     - MER position in degrees, and its distance from your input position
       (positional matches)
   * - ``flux_vis_psf``, ``flux_y_templfit``, ``flux_j_templfit``,
       ``flux_h_templfit``
     - VIS PSF flux and NISP template-fit fluxes, each with a ``fluxerr_``
       column
   * - ``flux_<band>_ext_<survey>_templfit``
     - Template-fit fluxes in the ground-based bands (DECam, MegaCam, JPCAM,
       Pan-STARRS, HSC), with errors
   * - ``point_like_prob``, ``extended_prob``, ``mu_max``,
       ``mumax_minus_mag``, ``kron_radius``
     - Star and galaxy separation and size
   * - ``det_quality_flag``, ``spurious_flag``, ``flag_vis``, ``flag_y``,
       ``flag_j``, ``flag_h``, ``vis_det``
     - Quality flags
   * - ``gaia_id``, ``gaia_match_quality``, ``parent_id``
     - Gaia association and deblending parent
   * - ``segmentation_map_id``, ``segmentation_area``
     - Needed for :doc:`segmentation_maps`

For the definition and units of each column, see the MER catalogue description
in the
`Euclid DR1 Data Product Description <https://euclid.esac.esa.int/dr/dr1/dpdd/>`_.

Some MER columns are empty for a given release or field. Add
``--drop-empty-columns`` to remove columns in which every value is missing.
Columns that contain zeros, ``False``, or empty strings are kept.

DEEP fields
~~~~~~~~~~~

.. code-block:: bash

   euclidkit crossmatch -i my_sources.fits -o deep_matches.fits \
     --field DEEP --deep-partition both

See :doc:`environments` for the partitions.

Large tables
~~~~~~~~~~~~

By default, ``euclidkit`` uploads your table in batches and runs short
synchronous queries. For large tables, add ``--full-async`` to run
asynchronous archive jobs instead:

.. code-block:: bash

   euclidkit crossmatch -i huge_sources.fits -o huge_matches.fits \
     --match-mode object-id --full-async

- A table with up to ``--async-chunk-size`` rows (default 500,000) runs as one
  job. ``euclidkit`` downloads the result and deletes the remote job.
- A larger table is split into chunks. Each chunk is saved as
  ``<output>_part_####.fits``, and progress is recorded in
  ``<output>.manifest.json``. When all chunks are done, they are merged into
  the output file.

``--max-sources N`` limits the run to the first ``N`` input rows, which is
useful for a test.

Reuse a table that is already in the archive
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you crossmatch the same large table several times, upload it once to your
archive user space:

.. code-block:: bash

   euclidkit upload-table -i my_sources.fits -t my_sources \
     --description "Quasar candidates"

Then crossmatch it by name, without uploading again:

.. code-block:: bash

   euclidkit crossmatch --user-table-name my_sources -o matches.fits \
     --match-mode object-id

``upload-table`` refuses to replace an existing table unless you add
``--overwrite``.

.. _morphology-guide:

Get morphology measurements
---------------------------

``query-morph`` adds MER morphology to a table that has an ``object_id``
column. You do not need to run ``crossmatch`` first if you already have the
IDs. It is available for DR1 and the internal archive, not for Q1.

.. code-block:: bash

   euclidkit query-morph -i matches.fits -o morphology.fits

The output is an inner join: it contains only the input rows that have a
morphology measurement. Each row gains these columns, and an ``_err`` column
for each:

.. list-table::
   :header-rows: 1
   :widths: 40 40 20

   * - Column
     - Measurement
     - Unit
   * - ``concentration``, ``asymmetry``, ``smoothness``
     - CAS parameters
     - None
   * - ``gini``, ``moment_20``
     - Gini coefficient and M20
     - None
   * - ``sersic_sersic_vis_radius``
     - Sérsic radius in VIS
     - arcsec
   * - ``sersic_sersic_vis_index``, ``sersic_sersic_vis_axis_ratio``
     - Sérsic index and axis ratio in VIS
     - None
   * - ``sersic_angle``
     - Sérsic position angle
     - deg
   * - ``flux_radius``
     - Flux radius
     - pixel

Things to know:

- The ``object_id`` column name is matched without regard to case. IDs must be
  integers or integer strings. Floating-point IDs are rejected, because a
  float cannot hold a Euclid ID exactly.
- If your table already has a column with one of these names, it is replaced.
- Repeated input IDs each keep their own row.
- ``--field``, ``--deep-partition``, ``--full-async``, and
  ``--async-chunk-size`` work as for ``crossmatch``.
- If a download fails in async mode, ``euclidkit`` writes a ``.job.json`` or
  ``.manifest.json`` recovery file and does not write a partial output table.

The measurements are defined in the
`DR1 MER morphology description <http://st-dm.pages.euclid-sgs.uk/data-product-doc/dmdr1/merdpd/dpcards/mer_finalcatalog.html#morphology-catalog-fits-file>`_.

.. _redshift-guide:

Get spectroscopic redshift candidates
-------------------------------------

``query-zspe`` adds SPE redshift candidates to a table that has an
``object_id`` column. It is available for DR1 and the internal archive, not
for Q1.

.. code-block:: bash

   euclidkit query-zspe -x matches.fits -o redshifts.fits --object-type qso

The SPE pipeline fits each spectrum under several object-type hypotheses.
Choose the hypothesis with ``--object-type qso`` (default) or
``--object-type galaxy``.

The output joins your input rows to the candidate solutions. The columns
include ``spe_rank``, ``spe_z``, and ``spe_z_err``. An object can have several
candidates. ``spe_rank`` orders them.

For the WIDE field, ``euclidkit`` searches the ``wide_survey`` candidate table
first and then the ``wide`` table for the objects that are still unmatched.
Use ``--field DEEP`` for the deep candidate table. ``--full-async`` and
``--async-chunk-size`` work as for ``crossmatch``.
