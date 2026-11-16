Archive Services and Data Releases
==================================

Archive service and public data release are separate choices. ``--environment``
selects the service; ``--release`` selects a dataset on the public service.

* ``PDR`` (default): public service. ``--release DR1`` is the default;
  ``--release Q1`` selects the earlier public release.
* ``IDR``: independent internal service (consortium access). Omit ``--release``.
  It is not an alias for public DR1, and this package does not yet advertise DR2.
* ``OTF`` and ``REG``: existing development/testing services. Omit ``--release``.

Python defaults to ``EuclidArchive(environment='PDR', data_release='DR1')``;
``EuclidArchive(data_release='Q1')`` selects public Q1. An explicit public release
with IDR, OTF or REG is rejected rather than redirected. Existing authentication
and async-job ownership behavior is unchanged.

Public table routing
--------------------

Public DR1 uses the ``dr1`` schema; public Q1 uses ``q1``. This applies consistently
to MER, spectra-source, segmentation-map, mosaic and raw-frame queries and the
default download/datalink schema. Explicit table or download-schema overrides
remain available. Query manifests and recovery metadata record both service and
release. IDR retains its current internal table routing separately; future IDR
schema changes will require validation against that service.

DR1 and IDR field selection
---------------------------

Use ``--field WIDE`` (default) or ``--field DEEP``. The legacy ``--idr-field``
spelling remains an alias, including for public DR1.

* WIDE queries ``dr1.mer_catalogue_wide_survey`` first, then
  ``dr1.mer_catalogue_wide_mode`` only for unmatched sources.
* DEEP uses ``--deep-partition survey`` (default), ``mode`` or ``both`` for
  ``dr1.mer_catalogue_deep_survey`` and/or ``dr1.mer_catalogue_deep_mode``.
  Survey takes precedence when an object occurs in both.
* Morphology uses the corresponding ``dr1.mer_morphology_*`` tables.
* Spectra sources use ``dr1.spectra_source_wide`` or ``dr1.spectra_source_deep``.
  SPE candidates retain the WIDE survey-first fallback.

``--idr-deep-partition`` remains an alias for ``--deep-partition``. Existing
Python ``idr_field`` and ``idr_deep_partition`` keywords remain supported.
IDR uses the corresponding internal ``catalogue.*`` catalogue tables.
Q1 has a single MER catalogue and does not support DEEP selection, morphology
queries or SPE candidate queries through these commands.

Examples
--------

.. code-block:: bash

   # Public DR1 WIDE (both choices are defaults)
   euclidkit crossmatch -i sources.fits -o dr1_matches.fits

   # Public DR1 DEEP
   euclidkit query-cutana -s sources.fits -o dr1_deep.csv \
     --environment PDR --release DR1 --field DEEP --deep-partition both

   # Earlier public Q1
   euclidkit crossmatch -i sources.fits -o q1_matches.fits --release Q1

   # Explicit internal service: no public-release flag
   euclidkit query-morph -i sources.fits -o internal_morphology.fits \
     --environment IDR --field WIDE

Migration
---------

Commands that previously used PDR to access Q1 must now pass ``--release Q1``.
Public DR1 workflows previously run against IDR should use PDR (or omit the
service option). Old field-option names remain accepted. Public crossmatch
outputs use the requested filename; the legacy automatic WIDE/DEEP filename
prefix remains restricted to explicit IDR commands.
