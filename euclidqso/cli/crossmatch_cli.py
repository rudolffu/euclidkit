#!/usr/bin/env python3
"""
Command line interface for crossmatching functionality.
"""

import click
import sys
from pathlib import Path
from typing import Optional

from euclidqso.core.data_access import EuclidArchive
from euclidqso.core.spectra import SpectrumCompiler
from euclidqso.utils.io import load_table


@click.command()
@click.option('--input', '-i', required=True, type=click.Path(exists=True),
              help='Input source table (CSV, FITS, or VOTable)')
@click.option('--output', '-o', required=True, type=click.Path(),
              help='Output crossmatch results file')
@click.option('--radius', '-r', type=float, default=1.0,
              help='Search radius in arcseconds (default: 1.0)')
@click.option('--ra-col', type=str, default='ra',
              help='RA column name in input table (default: ra)')
@click.option('--dec-col', type=str, default='dec', 
              help='Dec column name in input table (default: dec)')
@click.option('--environment', '-e', type=click.Choice(['PDR', 'IDR', 'OTF', 'REG']),
              default='PDR', help='Archive environment (default: PDR)')
@click.option('--idr-field', type=click.Choice(['WIDE', 'DEEP']), default='WIDE',
              show_default=True,
              help='IDR field selection (only used when --environment=IDR)')
@click.option('--credentials', '-c', type=click.Path(exists=True),
              help='Credentials file path')
@click.option('--max-sources', type=int,
              help='Maximum number of sources to process')
@click.option('--match-mode', type=click.Choice(['auto', 'object-id', 'spatial']), default='auto',
              help='Matching mode: auto (default), object-id, or spatial')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
def crossmatch(input: str, output: str, radius: float, ra_col: str, dec_col: str,
               environment: str, idr_field: str, credentials: Optional[str],
               max_sources: Optional[int], match_mode: str, verbose: bool):
    """
    Crossmatch user source table with Euclid MER catalogue.
    
    This command crossmatches a user-provided source table with the Euclid
    MER catalogue using position-based matching within a specified radius.
    """
    import logging
    if verbose:
        logging.basicConfig(level=logging.INFO)
    
    try:
        # Initialize archive client
        archive = EuclidArchive(environment=environment)
        
        # Login
        if credentials:
            archive.login(credentials_file=credentials)
        else:
            archive.login()  # Use default credentials
        
        selected_idr_field = idr_field.upper()
        if verbose:
            click.echo(f"Connected to {environment} environment")
            click.echo(f"Input table: {input}")
            click.echo(f"Search radius: {radius} arcsec")
            click.echo(f"Match mode: {match_mode}")
            if max_sources:
                click.echo(f"Processing max {max_sources} sources")
            if environment == 'IDR':
                click.echo(f"IDR field: {selected_idr_field}")

        # Determine effective output path (IDR requires prefixed filenames)
        output_path = Path(output)
        effective_output_path = output_path
        if environment == 'IDR':
            prefix = f"{selected_idr_field.lower()}_"
            if not output_path.name.lower().startswith(prefix):
                effective_output_path = output_path.with_name(prefix + output_path.name)

        # Map match mode to use_object_id flag
        use_object_id = None
        if match_mode == 'object-id':
            use_object_id = True
        elif match_mode == 'spatial':
            use_object_id = False

        # Perform crossmatch
        crossmatch_kwargs = dict(
            user_table=input,
            radius=radius,
            output_file=effective_output_path,
            ra_col=ra_col,
            dec_col=dec_col,
            max_sources=max_sources,
            use_object_id=use_object_id
        )
        if environment == 'IDR':
            crossmatch_kwargs['idr_field'] = selected_idr_field

        results = archive.crossmatch_sources(**crossmatch_kwargs)
        
        # Report results
        click.echo(f"Crossmatch completed: {len(results)} matches found")
        click.echo(f"Results saved to: {effective_output_path}")
        
        # Show summary statistics
        if len(results) > 0 and 'separation_arcsec' in results.colnames:
            separations = results['separation_arcsec']
            click.echo(f"Separation statistics (arcsec):")
            click.echo(f"  Min: {separations.min():.3f}")
            click.echo(f"  Max: {separations.max():.3f}")
            click.echo(f"  Mean: {separations.mean():.3f}")
            click.echo(f"  Median: {separations[len(separations)//2]:.3f}")
        
    except Exception as e:
        click.echo(f"Error in crossmatch: {e}", err=True)
        sys.exit(1)
    
    finally:
        archive.logout()


@click.command()
@click.option('--crossmatch', '-x', required=True, type=click.Path(exists=True),
              help='Input crossmatch results file (must contain Euclid object_id)')
@click.option('--output', '-o', required=True, type=click.Path(),
              help='Output spectral sources file')
@click.option('--combine-output', type=click.Path(),
              help='Output FITS file for combined spectra (auto-generates after query)')
@click.option('--max-spectra', type=int,
              help='Maximum number of spectra to include in combined output')
@click.option('--environment', '-e', type=click.Choice(['PDR', 'IDR', 'OTF', 'REG']),
              default='PDR', help='Archive environment (default: PDR)')
@click.option('--credentials', '-c', type=click.Path(exists=True),
              help='Credentials file path')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
def query_spectra(crossmatch: Optional[str], output: str,
                  combine_output: Optional[str], max_spectra: Optional[int],
                  environment: str, credentials: Optional[str], verbose: bool):
    """
    Query spectral sources for objects from crossmatch or object ID list.
    
    This command queries the sedm.spectra_source table to find available
    spectra for objects identified in crossmatching or provided as a list.
    
    If --combine-output is provided, automatically combines the found spectra
    into a single FITS file after querying, similar to cell 23 in the 
    Spectra_visualization_catglobe.ipynb notebook.
    """
    import logging
    if verbose:
        logging.basicConfig(level=logging.INFO)
    
    # crossmatch is required by Click; double-check for clarity
    if not crossmatch:
        click.echo("Error: Must provide --crossmatch file", err=True)
        sys.exit(1)
    
    try:
        # Initialize archive client
        archive = EuclidArchive(environment=environment)
        
        # Login
        if credentials:
            archive.login(credentials_file=credentials)
        else:
            archive.login()  # Use default credentials
        
        if verbose:
            click.echo(f"Connected to {environment} environment")
        
        # Load crossmatch table
        crossmatch_table = load_table(crossmatch)
        if verbose:
            click.echo(f"Using crossmatch table: {crossmatch}")
        
        # Query spectra
        results = archive.query_spectra_sources(
            crossmatch_table=crossmatch_table,
            output_file=output
        )
        
        # Report results
        click.echo(f"Spectral query completed: {len(results)} spectra found")
        click.echo(f"Results saved to: {output}")
        
        if len(results) > 0:
            unique_objects = len(set(results['object_id']))
            click.echo(f"Unique objects with spectra: {unique_objects}")
            
            # Show instrument breakdown
            if 'instrument_name' in results.colnames:
                from collections import Counter
                instruments = Counter(results['instrument_name'])
                click.echo("Spectra by instrument:")
                for instrument, count in instruments.items():
                    click.echo(f"  {instrument}: {count}")
            
            # Generate combined FITS file if requested
            if combine_output:
                if verbose:
                    click.echo(f"\nGenerating combined spectra FITS file...")
                    if max_spectra:
                        click.echo(f"Limiting to {max_spectra} spectra")
                
                try:
                    combined_file = archive.combine_spectra_to_fits(
                        spectra_table=results,
                        output_file=combine_output,
                        max_spectra=max_spectra
                    )
                    
                    n_combined = min(len(results), max_spectra) if max_spectra else len(results)
                    click.echo(f"Combined FITS file created: {combined_file}")
                    click.echo(f"Contains {n_combined} spectra extensions")
                    
                except Exception as e:
                    click.echo(f"Warning: Failed to create combined FITS file: {e}", err=True)
                    # Don't exit - the query was successful
        
    except Exception as e:
        click.echo(f"Error querying spectra: {e}", err=True)
        sys.exit(1)
    
    finally:
        archive.logout()


@click.command()
@click.option('--spectra-table', '-s', required=True, type=click.Path(exists=True),
              help='Spectral sources table from query-spectra command')
@click.option('--output-dir', '-o', required=True, type=click.Path(),
              help='Output directory for compiled FITS files')
@click.option('--prefix', type=str, default='compiled_spectra',
              help='Prefix for output files (default: compiled_spectra)')
@click.option('--max-extensions', type=int, default=5000,
              help='Maximum extensions per file (default: 5000)')
@click.option('--overwrite', is_flag=True,
              help='Overwrite existing output files')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
def compile_spectra(spectra_table: str, output_dir: str, prefix: str,
                   max_extensions: int, overwrite: bool, verbose: bool):
    """
    Compile individual spectra into chunked multi-extension FITS files.
    
    This command takes the output from query-spectra and creates compiled
    FITS files with multiple spectrum extensions, following the pattern
    from your notebook example.
    """
    import logging
    if verbose:
        logging.basicConfig(level=logging.INFO)
    
    try:
        # Load spectral sources table
        sources = load_table(spectra_table)
        if verbose:
            click.echo(f"Loaded {len(sources)} spectral sources from {spectra_table}")
        
        # Initialize compiler
        compiler = SpectrumCompiler(max_extensions=max_extensions)
        
        if verbose:
            click.echo(f"Max extensions per file: {max_extensions}")
            click.echo(f"Output directory: {output_dir}")
        
        # Compile spectra
        output_files = compiler.compile_spectra(
            spectra_table=sources,
            output_dir=output_dir,
            output_prefix=prefix,
            overwrite=overwrite
        )
        
        # Create metadata table
        metadata_file = compiler.create_metadata_table(
            spectra_table=sources,
            output_files=output_files,
            output_dir=output_dir,
            output_name=f"{prefix}_metadata.fits"
        )
        
        # Report results
        click.echo(f"Compilation completed successfully!")
        click.echo(f"Created {len(output_files)} FITS files:")
        for i, file_path in enumerate(output_files, 1):
            file_name = Path(file_path).name
            click.echo(f"  {i:2d}. {file_name}")
        
        click.echo(f"Metadata saved to: {Path(metadata_file).name}")
        click.echo(f"Total spectra processed: {len(sources)}")
        
    except Exception as e:
        click.echo(f"Error compiling spectra: {e}", err=True)
        sys.exit(1)


# Main command group for crossmatching functionality
@click.group()
def crossmatch_commands():
    """Crossmatching and spectral compilation commands."""
    pass


# Add commands to group
crossmatch_commands.add_command(crossmatch, name='crossmatch')
crossmatch_commands.add_command(query_spectra, name='query-spectra') 
crossmatch_commands.add_command(compile_spectra, name='compile-spectra')


if __name__ == '__main__':
    crossmatch_commands()
