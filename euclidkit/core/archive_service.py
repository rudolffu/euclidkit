"""Archive service and public-release table routing, independent of authentication."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ArchiveService:
    environment: str = 'PDR'
    data_release: str | None = None

    def __post_init__(self):
        environment = str(self.environment).upper()
        if environment not in {'PDR', 'IDR', 'OTF', 'REG'}:
            raise ValueError('environment must be one of: PDR, IDR, OTF, REG')
        release = str(self.data_release).upper() if self.data_release is not None else None
        if environment != 'PDR' and release is not None:
            raise ValueError('Named public releases require environment=PDR; omit data_release for internal services')
        if environment == 'PDR':
            release = release or 'DR1'
            if release not in {'DR1', 'Q1'}:
                raise ValueError('data_release must be DR1 or Q1')
        object.__setattr__(self, 'environment', environment)
        object.__setattr__(self, 'data_release', release)

    @property
    def has_fields(self):
        return self.environment == 'IDR' or self.data_release == 'DR1'

    @property
    def catalogue_schema(self):
        return self.data_release.lower() if self.environment == 'PDR' else 'catalogue'

    @property
    def product_schema(self):
        return self.data_release.lower() if self.environment == 'PDR' else 'sedm'

    def mer_tables(self, field=None, partition='survey'):
        if not self.has_fields:
            if field is not None and str(field).upper() != 'WIDE':
                raise ValueError('WIDE/DEEP selection is available for public DR1 and IDR only')
            return [{'PDR': 'q1.mer_catalogue', 'OTF': 'catalogue.mer_catalogue',
                     'REG': 'catalogue.mer_final_catalog_fits_file_regreproc1_r2'}[self.environment]]
        field = (field or 'WIDE').upper()
        if field not in {'WIDE', 'DEEP'}:
            raise ValueError('field must be WIDE or DEEP')
        if field == 'WIDE':
            suffixes = ['wide_survey', 'wide_mode']
        else:
            partition = (partition or 'survey').lower()
            if partition not in {'survey', 'mode', 'both'}:
                raise ValueError('DEEP partition must be survey, mode or both')
            suffixes = ['deep_survey', 'deep_mode'] if partition == 'both' else ['deep_' + partition]
        return [f'{self.catalogue_schema}.mer_catalogue_{suffix}' for suffix in suffixes]

    def morphology_tables(self, field=None, partition='survey'):
        self.require_fields('Morphology queries')
        return [t.replace('.mer_catalogue_', '.mer_morphology_')
                for t in self.mer_tables(field, partition)]

    def require_fields(self, operation):
        if not self.has_fields:
            raise ValueError(f'{operation} support public DR1 (PDR) and IDR only; Q1 is not supported')

    def spectra_table(self, field=None):
        if self.has_fields:
            field = (field or 'WIDE').upper()
            if field not in {'WIDE', 'DEEP'}:
                raise ValueError('field must be WIDE or DEEP')
            return f'{self.catalogue_schema}.spectra_source_{field.lower()}'
        if field is not None and str(field).upper() != 'WIDE':
            raise ValueError('WIDE/DEEP selection is available for public DR1 and IDR only')
        prefix = {'PDR': 'q1', 'OTF': 'q1', 'REG': 'dr1'}[self.environment]
        return f'{prefix}.spectra_source'

    def product_table(self, name):
        if self.environment == 'PDR':
            return f'{self.product_schema}.{name}'
        # Existing internal-service routing; these are not public release selectors.
        if self.environment == 'IDR' and name in {'raw_frame', 'mer_segmentation_map'}:
            return f'dr1.{name}'
        return f'sedm.{name}'

    def zspe_tables(self, object_type, field=None):
        self.require_fields('SPE candidate queries')
        object_type = object_type.lower()
        if object_type not in {'qso', 'galaxy'}:
            raise ValueError('object_type must be one of: qso, galaxy')
        field = (field or 'WIDE').upper()
        if field not in {'WIDE', 'DEEP'}:
            raise ValueError('field must be WIDE or DEEP')
        prefix = f'{self.catalogue_schema}.spectro_z_spe_{object_type}_candidates'
        suffixes = ['wide_survey', 'wide'] if field == 'WIDE' else ['deep']
        return [(f'{prefix}_{s}', s) for s in suffixes]


def product_schema_for_client(client, schema=None, data_release=None):
    """Resolve a datalink schema while preserving explicit schema overrides."""
    environment = getattr(client, 'environment', 'PDR')
    if not isinstance(environment, str):
        environment = 'PDR'
    service = ArchiveService(environment, data_release)
    return schema if schema is not None else service.product_schema
