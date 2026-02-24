#!/usr/bin/env python3
"""
Command line interface for crossmatching functionality.
"""

import click
import sys
from pathlib import Path
from typing import Optional


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
@click.option('--full-async', is_flag=True,
              help='Upload entire table in a single asynchronous job (no batching)')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
def crossmatch(input: str, output: str, radius: float, ra_col: str, dec_col: str,
               environment: str, idr_field: str, credentials: Optional[str],
               max_sources: Optional[int], match_mode: str, full_async: bool,
               verbose: bool):
    """
    Crossmatch user source table with Euclid MER catalogue.
    
    This command crossmatches a user-provided source table with the Euclid
    MER catalogue using position-based matching within a specified radius.
    """
    import logging
    from euclidkit.core.data_access import EuclidArchive

    if verbose:
        logging.basicConfig(level=logging.INFO)
    
    archive = None

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
            click.echo(f"Match mode: {match_mode}")
            if match_mode == 'object-id':
                click.echo("Search radius: n/a (object-id mode)")
            else:
                click.echo(f"Search radius: {radius} arcsec")
            if max_sources:
                click.echo(f"Processing max {max_sources} sources")
            if environment == 'IDR':
                click.echo(f"IDR field: {selected_idr_field}")
            if full_async:
                click.echo("Full-table async mode enabled (no batching)")

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
            use_object_id=use_object_id,
            full_async=full_async,
        )
        if environment == 'IDR':
            crossmatch_kwargs['idr_field'] = selected_idr_field

        results = archive.crossmatch_sources(**crossmatch_kwargs)
        
        if full_async:
            if results.get('results_downloaded'):
                click.echo("Crossmatch async query completed and results were downloaded.")
                click.echo(f"Results saved to: {effective_output_path}")
                if results.get('result_row_count') is not None:
                    click.echo(f"Rows downloaded: {results['result_row_count']}")
                job_id = results.get('job_id')
                if job_id:
                    click.echo(f"Job ID: {job_id}")
            else:
                job_id = results.get('job_id')
                click.echo("Crossmatch job submitted asynchronously.")
                if job_id:
                    click.echo(f"Job ID: {job_id}")
                job_info_file = results.get('job_info_file', str(effective_output_path))
                click.echo(f"Job info saved to: {job_info_file}")
        else:
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
        if archive is not None:
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
    
    This command queries the environment-specific spectra_source table to find available
    spectra for objects identified in crossmatching or provided as a list.

    Credits: Kristin Anett Remmelgas and Héctor Cánovas Cabrera.
    
    If --combine-output is provided, automatically combines the found spectra
    into a single FITS file after querying, similar to cell 23 in the 
    Spectra_visualization_catglobe.ipynb notebook.
    """
    import logging
    from euclidkit.core.data_access import EuclidArchive
    from euclidkit.utils.io import load_table

    archive = None

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
        if archive is not None:
            archive.logout()


@click.command(name='query-cutana')
@click.option('--sources', '-s', required=True, type=click.Path(exists=True),
              help='Input table with object_id or RA/Dec columns')
@click.option('--output', '-o', required=True, type=click.Path(),
              help='Output CSV file path for Cutana input')
@click.option('--instrument', type=click.Choice(['VIS', 'NISP']), default='VIS',
              show_default=True, help='Instrument used for mosaic selection')
@click.option('--nisp-filters', type=str,
              help='Comma-separated NISP filters (e.g. NIR_Y,NIR_H)')
@click.option('--cutout-size', type=click.Choice(['pixel', 'arcsec']), default='arcsec',
              show_default=True, help='Diameter unit for Cutana file')
@click.option('--cutout-size-value', type=float, default=15.0, show_default=True,
              help='Constant diameter value for all sources')
@click.option('--drop-noncutana-cols/--keep-noncutana-cols', default=True, show_default=True,
              help='Drop or keep non-Cutana columns from the input table')
@click.option('--environment', '-e', type=click.Choice(['PDR', 'IDR', 'OTF', 'REG']),
              default='PDR', help='Archive environment (default: PDR)')
@click.option('--idr-field', type=click.Choice(['WIDE', 'DEEP']), default='WIDE',
              show_default=True,
              help='IDR field selection (only used when --environment=IDR)')
@click.option('--credentials', '-c', type=click.Path(exists=True),
              help='Credentials file path')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
def query_cutana(sources: str, output: str, instrument: str, nisp_filters: Optional[str],
                 cutout_size: str, cutout_size_value: float, drop_noncutana_cols: bool,
                 environment: str, idr_field: str, credentials: Optional[str], verbose: bool):
    """
    Generate a Cutana input catalogue from a source table.

    The input table can contain either Euclid ``object_id`` values or source
    coordinates (``ra``/``dec`` or ``right_ascension``/``declination``).
    """
    import logging
    from euclidkit.core.cutouts import CutoutGenerator
    from euclidkit.utils.io import load_table

    if verbose:
        logging.basicConfig(level=logging.INFO)

    generator = CutoutGenerator(environment=environment)
    archive = generator.archive

    try:
        if credentials:
            archive.login(credentials_file=credentials)
        else:
            archive.login()

        input_table = load_table(sources)
        filters = None
        if nisp_filters:
            filters = [f.strip() for f in nisp_filters.split(',') if f.strip()]
            if instrument != 'NISP':
                raise ValueError("--nisp-filters requires --instrument NISP")

        result_df = generator.generate_cutana_input(
            sources=input_table,
            output_file=output,
            instrument_name=instrument,
            nisp_filters=filters,
            cutout_size=cutout_size,
            cutout_size_value=cutout_size_value,
            drop_noncutana_cols=drop_noncutana_cols,
            idr_field=idr_field.upper() if environment == 'IDR' else None,
        )

        click.echo(f"Cutana query completed: {len(result_df)} sources with mosaic matches")
        click.echo(f"Cutana input file saved to: {output}")
    except Exception as e:
        click.echo(f"Error generating Cutana input: {e}", err=True)
        sys.exit(1)
    finally:
        archive.logout()


@click.command(name='upload-table')
@click.option('--input', '-i', required=True, type=click.Path(exists=True),
              help='Local table file to upload (FITS, CSV, VOTable, etc.)')
@click.option('--table-name', '-t', required=True,
              help='Name for the table in the Euclid user workspace')
@click.option('--description', '-d', type=str,
              help='Optional table description stored in TAP metadata')
@click.option('--format', '-f', 'fmt', type=str,
              help='Explicit input format (e.g., votable, fits, csv)')
@click.option('--overwrite/--no-overwrite', default=False, show_default=True,
              help='Delete any existing table with the same name before uploading')
@click.option('--environment', '-e', type=click.Choice(['PDR', 'IDR', 'OTF', 'REG']),
              default='PDR', help='Archive environment (default: PDR)')
@click.option('--credentials', '-c', type=click.Path(exists=True),
              help='Credentials file path')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
def upload_table(input: str, table_name: str, description: Optional[str], fmt: Optional[str],
                 overwrite: bool, environment: str, credentials: Optional[str], verbose: bool):
    """
    Upload a local table to the Euclid TAP user workspace.
    """
    import logging
    from euclidkit.core.data_access import EuclidArchive

    archive = None

    if verbose:
        logging.basicConfig(level=logging.INFO)
    
    try:
        archive = EuclidArchive(environment=environment)
        if credentials:
            archive.login(credentials_file=credentials)
        else:
            archive.login()

        result = archive.upload_user_table(
            table=input,
            table_name=table_name,
            description=description,
            fmt=fmt,
            overwrite=overwrite,
            verbose=verbose,
        )

        job_id = result.get('job_id')
        if job_id:
            click.echo(f"Upload job submitted (ID: {job_id}). Use astroquery TAP tools to monitor completion.")
        else:
            click.echo("Table uploaded successfully.")
        click.echo(f"Table name: {table_name}")
        click.echo(f"Format: {result.get('format')}")
        if description:
            click.echo(f"Description: {description}")

    except Exception as e:
        click.echo(f"Error uploading table: {e}", err=True)
        sys.exit(1)
    finally:
        if archive is not None:
            archive.logout()


@click.command()
@click.option('--spectra-table', '-s', required=True, type=click.Path(exists=True),
              help='Spectral sources table from query-spectra command')
@click.option('--output-dir', '-o', required=True, type=click.Path(),
              help='Output directory for compiled FITS files')
@click.option('--prefix', type=str, default='compiled_spectra',
              help='Prefix for output files (default: compiled_spectra)')
@click.option('--max-extensions', type=int, default=1000,
              help='Maximum extensions per file (default: 1000)')
@click.option('--overwrite', is_flag=True,
              help='Overwrite existing output files')
@click.option('--use-datalink', is_flag=True,
              help='Retrieve spectra via Euclid datalink instead of local datalabs_path/file_name')
@click.option('--environment', '-e', type=click.Choice(['PDR', 'IDR', 'OTF', 'REG']),
              default='PDR', help='Archive environment for datalink mode (default: PDR)')
@click.option('--credentials', '-c', type=click.Path(exists=True),
              help='Credentials file path for datalink mode')
@click.option('--retrieval-type', type=click.Choice(['ALL', 'SPECTRA_BGS', 'SPECTRA_RGS']),
              default='SPECTRA_RGS', show_default=True,
              help='Datalink retrieval type (used with --use-datalink)')
@click.option('--schema', type=str, default='sedm', show_default=True,
              help='Datalink schema value (used with --use-datalink)')
@click.option('--limit', type=int, default=None,
              help='Process only the first N rows from spectra table (useful for quick tests)')
@click.option('--workers', type=int, default=1, show_default=True,
              help='Number of chunk workers (canonical compile only)')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
def compile_spectra(spectra_table: str, output_dir: str, prefix: str,
                   max_extensions: int, overwrite: bool, use_datalink: bool,
                   environment: str, credentials: Optional[str], retrieval_type: str,
                   schema: str, limit: Optional[int], workers: int, verbose: bool):
    """
    Compile individual spectra into chunked multi-extension FITS files.
    
    This command takes the output from query-spectra and creates compiled
    FITS files with multiple spectrum extensions, following the pattern
    from your notebook example.
    """
    import logging
    from euclidkit.core.spectra import SpectrumCompiler
    from euclidkit.utils.io import load_table
    from euclidkit.core.data_access import EuclidArchive

    if verbose:
        logging.basicConfig(level=logging.INFO)
    
    archive = None
    try:
        target_output_dir = Path(output_dir)
        if target_output_dir.exists() and target_output_dir.is_dir():
            has_existing_content = any(target_output_dir.iterdir())
            if has_existing_content and overwrite:
                confirmed = click.confirm(
                    f"Output directory '{target_output_dir}' is not empty. Overwrite existing files?",
                    default=False
                )
                if not confirmed:
                    overwrite = False
            elif has_existing_content and not overwrite:
                click.echo(
                    f"Output directory '{target_output_dir}' is not empty. "
                    "Resume mode: existing chunk files will be kept."
                )
        elif target_output_dir.exists() and target_output_dir.is_file():
            raise ValueError(f"Output path '{target_output_dir}' is a file; please provide a directory path.")

        output_dir = str(target_output_dir)

        # Load spectral sources table
        sources = load_table(spectra_table)
        if limit is not None:
            if limit <= 0:
                raise ValueError("--limit must be a positive integer")
            sources = sources[:limit]
        if verbose:
            click.echo(f"Loaded {len(sources)} spectral sources from {spectra_table}")
        
        # Initialize compiler
        compiler = SpectrumCompiler(max_extensions=max_extensions)
        
        if verbose:
            click.echo(f"Max extensions per file: {max_extensions}")
            click.echo(f"Output directory: {output_dir}")
            if not use_datalink:
                click.echo(f"Chunk workers: {workers}")
        
        # Compile spectra
        if use_datalink:
            if workers != 1:
                click.echo("Note: --workers is currently applied to canonical compile only; datalink runs with one worker.")
            archive = EuclidArchive(environment=environment)
            if credentials:
                archive.login(credentials_file=credentials)
            else:
                archive.login()
            output_files = compiler.compile_spectra_datalink(
                spectra_table=sources,
                euclid_client=archive.euclid,
                output_dir=output_dir,
                output_prefix=prefix,
                retrieval_type=retrieval_type,
                schema=schema,
                overwrite=overwrite,
            )
        else:
            output_files = compiler.compile_spectra(
                spectra_table=sources,
                output_dir=output_dir,
                output_prefix=prefix,
                overwrite=overwrite,
                workers=workers,
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
    finally:
        if archive is not None:
            archive.logout()


# Main command group for crossmatching functionality
@click.group()
def crossmatch_commands():
    """Crossmatching and spectral compilation commands."""
    pass


# Add commands to group
crossmatch_commands.add_command(crossmatch, name='crossmatch')
crossmatch_commands.add_command(query_spectra, name='query-spectra') 
crossmatch_commands.add_command(query_cutana, name='query-cutana')
crossmatch_commands.add_command(compile_spectra, name='compile-spectra')


if __name__ == '__main__':
    crossmatch_commands()
