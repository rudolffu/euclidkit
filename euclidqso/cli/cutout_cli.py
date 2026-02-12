"""
Command line interface for cutout generation.
"""

import click
import os
import sys
from pathlib import Path
from typing import Optional, List
import json


@click.command()
@click.option(
    '--source', '-s',
    help='Single source: object_id or "ra,dec" (degrees)'
)
@click.option(
    '--sources-file', '-f',
    type=click.Path(exists=True),
    help='File with sources (CSV, JSON, or text file)'
)
@click.option(
    '--ra',
    type=float,
    help='Right Ascension in degrees (use with --dec)'
)
@click.option(
    '--dec', 
    type=float,
    help='Declination in degrees (use with --ra)'
)
@click.option(
    '--object-id',
    help='Euclid object ID'
)
@click.option(
    '--output', '-o',
    required=True,
    type=click.Path(),
    help='Output directory for cutouts'
)
@click.option(
    '--size',
    default=10.0,
    type=float,
    help='Cutout size in arcseconds (default: 10)'
)
@click.option(
    '--instrument',
    default='VIS',
    type=click.Choice(['VIS', 'NISP']),
    help='Instrument name (default: VIS)'
)
@click.option(
    '--product-type',
    default='mosaic',
    type=click.Choice(['mosaic', 'calibrated_frame', 'stacked_frame']),
    help='Product type (default: mosaic)'
)
@click.option(
    '--environment',
    default='PDR',
    type=click.Choice(['PDR', 'OTF', 'REG', 'IDR']),
    help='Archive environment (default: PDR)'
)
@click.option(
    '--nisp-filters',
    help='NISP filters (comma-separated): NIR_H,NIR_J,NIR_Y'
)
@click.option(
    '--verbose', '-v',
    is_flag=True,
    help='Verbose output'
)
def cutouts(
    source: Optional[str],
    sources_file: Optional[str],
    ra: Optional[float],
    dec: Optional[float], 
    object_id: Optional[str],
    output: str,
    size: float,
    instrument: str,
    product_type: str,
    environment: str,
    nisp_filters: Optional[str],
    verbose: bool
):
    """
    Generate Euclid image cutouts.
    
    Examples:
    
    \b
    # Single source by object ID
    euclidkit-cutouts --object-id 123456789 --output /tmp/cutouts/
    
    \b  
    # Single source by coordinates
    euclidkit-cutouts --ra 150.0 --dec 2.5 --output /tmp/cutouts/
    
    \b
    # Multiple sources from file
    euclidkit-cutouts --sources-file sources.csv --output /tmp/cutouts/
    
    \b
    # Custom size and instrument
    euclidkit-cutouts --source "150.0,2.5" --output /tmp/cutouts/ --size 20 --instrument NISP
    """
    
    # Set up logging
    import logging
    import astropy.units as u
    from euclidqso.core.cutouts import CutoutGenerator

    if verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)
    
    # Determine source input
    sources = None
    
    # Priority: sources_file > ra/dec > object_id > source
    if sources_file:
        sources = _load_sources_file(sources_file)
    elif ra is not None and dec is not None:
        sources = {'ra': ra, 'dec': dec}
    elif object_id:
        sources = object_id
    elif source:
        # Parse source string
        if ',' in source:
            # Assume ra,dec format
            try:
                ra_val, dec_val = map(float, source.split(','))
                sources = {'ra': ra_val, 'dec': dec_val}
            except ValueError:
                click.echo(f"Error: Could not parse coordinates from '{source}'", err=True)
                sys.exit(1)
        else:
            # Assume object_id
            sources = source
    else:
        click.echo("Error: Must provide source(s) via --source, --sources-file, --ra/--dec, or --object-id", err=True)
        sys.exit(1)
    
    # Parse NISP filters
    nisp_filter_list = None
    if nisp_filters:
        nisp_filter_list = [f.strip() for f in nisp_filters.split(',')]
    
    # Create output directory
    os.makedirs(output, exist_ok=True)
    
    try:
        # Initialize cutout generator
        generator = CutoutGenerator(environment=environment)
        
        # Generate cutouts
        cutout_size = size * u.arcsec
        
        num_cutouts = generator.make_cutouts(
            sources=sources,
            output_location=output,
            cutout_size=cutout_size,
            instrument_name=instrument,
            product_type=product_type,
            nisp_filters=nisp_filter_list
        )
        
        if num_cutouts > 0:
            click.echo(f"Successfully generated {num_cutouts} cutouts in {output}")
        else:
            click.echo("No cutouts generated. Check your sources and try again.", err=True)
            sys.exit(1)
            
    except Exception as e:
        click.echo(f"Error generating cutouts: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def _load_sources_file(filepath: str):
    """Load sources from various file formats."""
    
    path = Path(filepath)
    
    if path.suffix.lower() == '.csv':
        import pandas as pd
        return pd.read_csv(filepath)
    
    elif path.suffix.lower() == '.json':
        with open(filepath, 'r') as f:
            data = json.load(f)
        return data
    
    elif path.suffix.lower() in ['.txt', '.dat']:
        # Try to parse as simple text file
        # Format: object_id OR ra dec (one per line)
        sources = []
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split()
                if len(parts) == 1:
                    # Assume object_id
                    sources.append({'object_id': parts[0]})
                elif len(parts) == 2:
                    # Assume ra dec
                    try:
                        ra, dec = map(float, parts)
                        sources.append({'ra': ra, 'dec': dec})
                    except ValueError:
                        continue
        
        return sources
    
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")


if __name__ == '__main__':
    cutouts()
