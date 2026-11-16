Glossary
========

Euclid and archive terms used in this documentation.

.. glossary::
   :sorted:

   ADQL
      Astronomical Data Query Language. The SQL-like language of :term:`TAP`
      services. ``euclidkit`` writes the ADQL for you.

   BGS
      The NISP blue grism. See also :term:`RGS`.

   Cutana
      An ESA tool that makes image cutouts in bulk from large collections of
      FITS tiles. ``euclidkit query-cutana`` prepares its input catalogue.

   Datalabs
      ESA Datalabs, a science platform where the Euclid data files are mounted
      as local volumes. Commands that read spectra or images in bulk run
      there.

   Datalink
      An archive service that delivers the data files linked to a catalogue
      row, such as the spectrum of a source. It works from any computer, but
      it is slower than reading the Datalabs volumes.

   DEEP
      The Euclid deep survey. In ``euclidkit``, the ``survey`` partition
      covers the Euclid Deep Fields (EDF-N, EDF-F, EDF-S) and the ``mode``
      partition covers CDFS and COSMOS. See :doc:`environments`.

   Dither
      One of the several offset exposures of a field. The combined spectrum of
      a source is built from its single-dither spectra.

   DR1
      Euclid Data Release 1. The default release in ``euclidkit``.

   IDR
      The internal archive service for Euclid Consortium members. It is
      separate from the public archive.

   MER
      The Euclid pipeline stage that merges VIS, NISP, and ground-based
      imaging into one object catalogue. The MER catalogue is the main Euclid
      source catalogue, and ``object_id`` is its identifier.

   MOC
      Multi-Order Coverage map. A HEALPix description of a sky region.
      ``euclidkit`` ships MOCs of the DR1 footprints.

   Mosaic
      A coadded image tile produced by the MER stage. Mosaic tiles overlap.

   NISP
      The Near-Infrared Spectrometer and Photometer. It takes images in the
      Y, J, and H bands and slitless spectra.

   object_id
      The MER catalogue identifier of a Euclid object. The spectra tables call
      the same value ``source_id``.

   PDR
      The public archive service. It serves the :term:`DR1` and :term:`Q1`
      releases.

   Q1
      The first Euclid Quick Release, which preceded DR1. Select it with
      ``--release Q1``.

   RGS
      The NISP red grism. Most Euclid spectra are taken with it.

   Segmentation map
      An image in which each pixel carries the label of the object that it
      belongs to.

   SIR
      The Euclid pipeline stage that extracts one-dimensional spectra from the
      NISP slitless exposures.

   SPE
      The Euclid pipeline stage that measures spectroscopic redshifts and
      spectral features from the :term:`SIR` spectra.

   TAP
      Table Access Protocol. The Virtual Observatory standard for querying
      archive tables. The Euclid archive is a TAP service.

   VIS
      The Euclid visible imager. It observes in one broad optical band.

   WIDE
      The Euclid wide survey. The default field in ``euclidkit``.
