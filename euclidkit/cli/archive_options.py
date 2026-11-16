"""Validate service/release choices before CLI authentication or output creation."""
from functools import wraps

import click

from euclidkit.core.archive_service import ArchiveService


def validate_archive_selection(callback):
    @wraps(callback)
    def checked(*args, **kwargs):
        try:
            service = ArchiveService(kwargs.get('environment', 'PDR'), kwargs.get('data_release'))
            if service.data_release == 'Q1':
                if kwargs.get('idr_field', 'WIDE') != 'WIDE' or kwargs.get('idr_deep_partition', 'survey') != 'survey':
                    raise ValueError('Field and DEEP partition selection are not available for Q1')
                if callback.__name__ in {'query_morph', 'query_zspe'}:
                    service.require_fields(callback.__name__)
        except ValueError as exc:
            raise click.BadParameter(str(exc), param_hint='--environment/--release') from exc
        return callback(*args, **kwargs)
    return checked
