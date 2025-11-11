#!/usr/bin/env python3
"""
Main command line interface for euclidqso package.
"""

import click
import sys
import os
from pathlib import Path
import yaml
from typing import Optional

from euclidqso.version import __version__
from euclidqso.config import generate_config_template
from euclidqso.cli.crossmatch_cli import crossmatch, query_spectra, compile_spectra, upload_table
from euclidqso.cli.cutout_cli import cutouts


@click.group()
@click.version_option(version=__version__)
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
@click.option('--config', '-c', type=click.Path(exists=True), help='Configuration file path')
@click.pass_context
def main(ctx, verbose: bool, config: Optional[str]):
    """
    euclidqso: Euclid QSO Analysis Package
    
    A comprehensive toolkit for analyzing QSO candidates from Euclid data,
    including spectroscopic and photometric analysis, template generation,
    and integration with external surveys like DESI.
    """
    # Ensure that ctx.obj exists and is a dict
    ctx.ensure_object(dict)
    ctx.obj['verbose'] = verbose
    ctx.obj['config'] = config
    
    if verbose:
        click.echo(f"euclidqso v{__version__} CLI initialized with config: {config}")


@main.command()
@click.option('--output', '-o', type=click.Path(), default='euclidqso_config.yaml',
              help='Output configuration file')
@click.option('--template', '-t', type=click.Choice(['basic', 'advanced', 'pipeline']),
              default='basic', help='Configuration template')
def init_config(output: str, template: str):
    """Generate a configuration file template."""
    try:
        config_template = generate_config_template(template)
        with open(output, 'w') as f:
            yaml.dump(config_template, f, default_flow_style=False, sort_keys=False)
        click.echo(f"Configuration template written to {output}")
    except Exception as e:
        click.echo(f"Error generating config: {e}", err=True)
        sys.exit(1)


@main.command()
@click.option('--check-deps', is_flag=True, help='Check dependencies')
@click.option('--check-data', is_flag=True, help='Check data access')
@click.option('--check-desi', is_flag=True, help='Check DESI integration')
def diagnostics(check_deps: bool, check_data: bool, check_desi: bool):
    """Run diagnostic checks for euclidqso installation."""
    if not any([check_deps, check_data, check_desi]):
        check_deps = check_data = check_desi = True
    
    results = {}
    
    # Check dependencies
    if check_deps:
        try:
            import numpy, scipy, matplotlib, astropy, astroquery, pandas
            import photutils, sep, skimage, tqdm, yaml, click, requests
            results['dependencies'] = {'passed': True, 'message': 'All core dependencies found'}
        except ImportError as e:
            results['dependencies'] = {'passed': False, 'message': f'Missing dependency: {e}'}
    
    # Check data access
    if check_data:
        q1_path = Path('/data/euclid_q1/')
        ero_path = Path('/data/euc_ero_data_01/')
        cred_path = Path('/media/user/cred.txt')
        
        if q1_path.exists() and ero_path.exists():
            results['data_access'] = {'passed': True, 'message': 'Euclid data volumes accessible'}
        else:
            results['data_access'] = {'passed': False, 'message': 'Euclid data volumes not found (expected on ESA Datalabs)'}
        
        if cred_path.exists():
            results['credentials'] = {'passed': True, 'message': 'Credentials file found'}
        else:
            results['credentials'] = {'passed': False, 'message': 'Credentials file not found at /media/user/cred.txt'}
    
    # Check DESI integration
    if check_desi:
        try:
            import sparcl.client
            results['desi'] = {'passed': True, 'message': 'DESI/SPARCL integration available'}
        except ImportError:
            results['desi'] = {'passed': False, 'message': 'DESI dependencies not installed (install with pip install euclidqso[desi])'}
    
    # Print results
    for check, result in results.items():
        status = "✓" if result['passed'] else "✗"
        click.echo(f"{status} {check}: {result['message']}")


def pipeline(config_file: str, **kwargs):
    """Pipeline entry point (placeholder for future implementation)."""
    click.echo(f"Pipeline functionality not yet implemented.")
    click.echo(f"Config file: {config_file}")
    sys.exit(1)


# Add commands to main CLI
main.add_command(crossmatch)
main.add_command(query_spectra, name='query-spectra')  
main.add_command(compile_spectra, name='compile-spectra')
main.add_command(upload_table)
main.add_command(cutouts)


if __name__ == '__main__':
    main()
