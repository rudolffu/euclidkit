#!/usr/bin/env python3
"""
Command line interface for crossmatching functionality.
"""

import click
import sys
from pathlib import Path
from typing import Optional


@click.command()
@click.option('--input', '-i', required=False, type=click.Path(exists=True),
              help='Input source table (CSV, FITS, or VOTable)')
@click.option('--user-table-name', type=str,
              help='Archive user table name (e.g. my_table for user_<username>.my_table)')
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
@click.option('--idr-deep-partition',
              type=click.Choice(['survey', 'mode', 'both']),
              default='survey',
              show_default=True,
              help='IDR DEEP MER partition to query')
@click.option('--credentials', '-c', type=click.Path(exists=True),
              help='Credentials file path')
@click.option('--max-sources', type=int,
              help='Maximum number of sources to process')
@click.option('--match-mode', type=click.Choice(['auto', 'object-id', 'spatial']), default='auto',
              help='Matching mode: auto (default), object-id, or spatial')
@click.option('--full-async', is_flag=True,
              help='Use asynchronous TAP mode; very large tables are split into async chunks')
@click.option('--async-chunk-size', type=int, default=500000, show_default=True,
              help='Rows per async chunk when --full-async is used on large tables')
@click.option('--drop-empty-columns', is_flag=True,
              help='Drop entirely null/missing columns from the final crossmatch output')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
def crossmatch(input: Optional[str], user_table_name: Optional[str], output: str, radius: float, ra_col: str, dec_col: str,
               environment: str, idr_field: str, idr_deep_partition: str, credentials: Optional[str],
               max_sources: Optional[int], match_mode: str, full_async: bool,
               async_chunk_size: int, drop_empty_columns: bool, verbose: bool):
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
        # Exactly one input mode is required.
        if bool(input) == bool(user_table_name):
            raise click.ClickException("Provide exactly one of --input or --user-table-name")

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
            if input:
                click.echo(f"Input table: {input}")
            else:
                click.echo(f"Archive user table: {user_table_name}")
            click.echo(f"Match mode: {match_mode}")
            if match_mode == 'object-id':
                click.echo("Search radius: n/a (object-id mode)")
            else:
                click.echo(f"Search radius: {radius} arcsec")
            if max_sources:
                click.echo(f"Processing max {max_sources} sources")
            if environment == 'IDR':
                click.echo(f"IDR field: {selected_idr_field}")
                if selected_idr_field == 'DEEP':
                    click.echo(f"IDR DEEP partition: {idr_deep_partition}")
            if full_async:
                click.echo("Full-table async mode enabled (no batching)")
                click.echo(f"Async chunk size: {async_chunk_size}")
            if drop_empty_columns:
                click.echo("Dropping entirely empty columns from final output")

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
        if user_table_name:
            crossmatch_kwargs = dict(
                user_table_name=user_table_name,
                radius=radius,
                output_file=effective_output_path,
                ra_col=ra_col,
                dec_col=dec_col,
                max_sources=max_sources,
                use_object_id=use_object_id,
                full_async=full_async,
                drop_empty_columns=drop_empty_columns,
            )
            if environment == 'IDR':
                crossmatch_kwargs['idr_field'] = selected_idr_field
                crossmatch_kwargs['idr_deep_partition'] = idr_deep_partition
            results = archive.crossmatch_user_table(**crossmatch_kwargs)
        else:
            crossmatch_kwargs = dict(
                user_table=input,
                radius=radius,
                output_file=effective_output_path,
                ra_col=ra_col,
                dec_col=dec_col,
                max_sources=max_sources,
                use_object_id=use_object_id,
                full_async=full_async,
                async_chunk_size=async_chunk_size,
                drop_empty_columns=drop_empty_columns,
            )
            if environment == 'IDR':
                crossmatch_kwargs['idr_field'] = selected_idr_field
                crossmatch_kwargs['idr_deep_partition'] = idr_deep_partition
            results = archive.crossmatch_sources(**crossmatch_kwargs)
        
        if full_async:
            if results.get('results_downloaded'):
                click.echo("Crossmatch async query completed and results were downloaded.")
                click.echo(f"Results saved to: {effective_output_path}")
                if results.get('result_row_count') is not None:
                    click.echo(f"Rows downloaded: {results['result_row_count']}")
                if results.get('chunk_count') is not None:
                    click.echo(f"Chunks processed: {results['chunk_count']} (chunk size: {results.get('chunk_size')})")
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
@click.option('--idr-field', type=click.Choice(['WIDE', 'DEEP']), default='WIDE',
              show_default=True,
              help='IDR field selection (only used when --environment=IDR)')
@click.option('--credentials', '-c', type=click.Path(exists=True),
              help='Credentials file path')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
def query_spectra(crossmatch: Optional[str], output: str,
                  combine_output: Optional[str], max_spectra: Optional[int],
                  environment: str, idr_field: str, credentials: Optional[str], verbose: bool):
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
            if environment == 'IDR':
                click.echo(f"IDR field: {idr_field.upper()}")
        
        # Load crossmatch table
        crossmatch_table = load_table(crossmatch)
        if verbose:
            click.echo(f"Using crossmatch table: {crossmatch}")
        
        # Query spectra
        results = archive.query_spectra_sources(
            crossmatch_table=crossmatch_table,
            output_file=output,
            idr_field=idr_field.upper() if environment == 'IDR' else None,
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


@click.command(name='query-zspe')
@click.option('--crossmatch', '-x', required=True, type=click.Path(exists=True),
              help='Input table containing Euclid object_id values')
@click.option('--output', '-o', required=True, type=click.Path(),
              help='Output SPE redshift candidate table')
@click.option('--object-type', type=click.Choice(['qso', 'galaxy']), default='qso',
              show_default=True, help='SPE candidate object type')
@click.option('--idr-field', type=click.Choice(['WIDE', 'DEEP']), default='WIDE',
              show_default=True, help='IDR field selection')
@click.option('--credentials', '-c', type=click.Path(exists=True),
              help='Credentials file path')
@click.option('--full-async', is_flag=True,
              help='Use asynchronous TAP mode with server-side WIDE fallback')
@click.option('--async-chunk-size', type=int, default=500000, show_default=True,
              help='Rows per async chunk when --full-async is used on large tables')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
def query_zspe(crossmatch: str, output: str, object_type: str, idr_field: str,
               credentials: Optional[str], full_async: bool, async_chunk_size: int,
               verbose: bool):
    """
    Query IDR SPE redshift candidates by object_id.

    WIDE queries search wide_survey first, then query wide only for objects
    without a wide_survey candidate.
    """
    import logging
    from euclidkit.core.data_access import EuclidArchive
    from euclidkit.utils.io import load_table

    archive = None

    if verbose:
        logging.basicConfig(level=logging.INFO)

    try:
        archive = EuclidArchive(environment='IDR')

        if credentials:
            archive.login(credentials_file=credentials)
        else:
            archive.login()

        if verbose:
            click.echo("Connected to IDR environment")
            click.echo(f"Input table: {crossmatch}")
            click.echo(f"Object type: {object_type}")
            click.echo(f"IDR field: {idr_field.upper()}")
            if full_async:
                click.echo("Full-table async mode enabled")
                click.echo(f"Async chunk size: {async_chunk_size}")

        crossmatch_table = load_table(crossmatch)
        results = archive.query_zspe_candidates(
            crossmatch_table=crossmatch_table,
            output_file=output,
            object_type=object_type,
            idr_field=idr_field.upper(),
            full_async=full_async,
            async_chunk_size=async_chunk_size,
        )

        if full_async:
            if results.get('results_downloaded'):
                click.echo("SPE redshift async query completed and results were downloaded.")
                click.echo(f"Results saved to: {output}")
                if results.get('result_row_count') is not None:
                    click.echo(f"Rows downloaded: {results['result_row_count']}")
                if results.get('chunk_count') is not None:
                    click.echo(f"Chunks processed: {results['chunk_count']} (chunk size: {results.get('chunk_size')})")
                job_id = results.get('job_id')
                if job_id:
                    click.echo(f"Job ID: {job_id}")
            else:
                click.echo("SPE redshift job submitted asynchronously.")
                job_id = results.get('job_id')
                if job_id:
                    click.echo(f"Job ID: {job_id}")
                job_info_file = results.get('job_info_file', str(output))
                click.echo(f"Job info saved to: {job_info_file}")
            return

        click.echo(f"SPE redshift query completed: {len(results)} candidates found")
        click.echo(f"Results saved to: {output}")

        if len(results) > 0 and 'object_id' in results.colnames:
            unique_objects = len(set(results['object_id']))
            click.echo(f"Unique objects with SPE redshifts: {unique_objects}")
        if len(results) > 0 and 'source_table' in results.colnames:
            click.echo("Candidates by source table:")
            for source in sorted(set(results['source_table'])):
                count = sum(results['source_table'] == source)
                click.echo(f"  {source}: {count}")

    except Exception as e:
        click.echo(f"Error querying SPE redshifts: {e}", err=True)
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
@click.option('--idr-deep-partition',
              type=click.Choice(['survey', 'mode', 'both']),
              default='survey',
              show_default=True,
              help='IDR DEEP MER partition to query')
@click.option('--credentials', '-c', type=click.Path(exists=True),
              help='Credentials file path')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
def query_cutana(sources: str, output: str, instrument: str, nisp_filters: Optional[str],
                 cutout_size: str, cutout_size_value: float, drop_noncutana_cols: bool,
                 environment: str, idr_field: str, idr_deep_partition: str,
                 credentials: Optional[str], verbose: bool):
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
            idr_deep_partition=idr_deep_partition,
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
              help='Output directory for compiled spectra')
@click.option('--prefix', type=str, default='compiled_spectra',
              help='Prefix for output files (default: compiled_spectra)')
@click.option('--output-format', type=click.Choice(['parquet', 'fits']), default='parquet',
              show_default=True,
              help='Output format for local Datalabs mode; ignored with --use-datalink')
@click.option('--max-extensions', type=int, default=1000,
              help='Maximum extensions per FITS file (FITS and Datalink modes only)')
@click.option('--chunk-size', type=int, default=2000, show_default=True,
              help='Rows per parquet part (parquet mode only)')
@click.option('--overwrite', is_flag=True,
              help='Overwrite existing output files for the selected prefix')
@click.option('--use-datalink', is_flag=True,
              help='Retrieve spectra via Euclid datalink instead of local datalabs_path/file_name')
@click.option('--environment', '-e', type=click.Choice(['PDR', 'IDR', 'OTF', 'REG']),
              default='PDR', help='Archive environment for datalink/FITS compatibility mode')
@click.option('--idr-field', type=click.Choice(['WIDE', 'DEEP']), default='WIDE',
              show_default=True,
              help='IDR field selection (used by FITS mode to resolve DEEP BGS/RGS XMLs)')
@click.option('--credentials', '-c', type=click.Path(exists=True),
              help='Credentials file path for datalink mode')
@click.option('--retrieval-type', type=click.Choice(['ALL', 'SPECTRA_BGS', 'SPECTRA_RGS']),
              default='SPECTRA_RGS', show_default=True,
              help='Datalink retrieval type (legacy; unified arm selection is controlled by --lambda-range/-L)')
@click.option('--schema', type=str, default='sedm', show_default=True,
              help='Datalink schema value (used with --use-datalink)')
@click.option('--limit', type=int, default=None,
              help='Process only the first N rows from spectra table (useful for quick tests)')
@click.option('--workers', type=int, default=None,
              help='Worker count: parquet defaults to min(os.cpu_count(), 8); FITS/Datalink default to 1')
@click.option('--lambda-range', '-L', type=click.Choice(['RGS', 'BGS', 'BOTH']),
              default='RGS', show_default=True,
              help='Spectral arm selection (RGS/BGS/BOTH)')
@click.option('--on-error', type=click.Choice(['fail', 'skip']), default='fail', show_default=True,
              help='Fail on bad parquet rows or skip them and write a failures JSONL')
@click.option('--progress', 'progress', is_flag=True, default=None,
              help='Show parquet export progress')
@click.option('--no-progress', 'progress', flag_value=False,
              help='Disable parquet export progress')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
@click.pass_context
def compile_spectra(ctx: click.Context, spectra_table: str, output_dir: str, prefix: str,
                   output_format: str, max_extensions: int, chunk_size: int, overwrite: bool,
                   use_datalink: bool, environment: str, idr_field: str,
                   credentials: Optional[str], retrieval_type: str, schema: str,
                   limit: Optional[int], workers: Optional[int], lambda_range: str,
                   on_error: str, progress: Optional[bool], verbose: bool):
    """
    Compile local Datalabs spectra to parquet by default, or FITS via compatibility modes.

    Non-Datalink mode defaults to raw parquet parts read directly from
    datalabs_path/file_name/hdu_index catalog rows. Use --output-format fits for
    the legacy multi-extension FITS compiler. --use-datalink remains FITS-only.
    """
    import logging
    from euclidkit.core.spectra import SpectrumCompiler
    from euclidkit.utils.io import load_table
    from euclidkit.core.data_access import EuclidArchive
    from euclidkit.core.spectra_parquet import default_raw_parquet_workers, spectra_to_parquet

    if verbose:
        logging.basicConfig(level=logging.INFO)

    archive = None
    try:
        target_output_dir = Path(output_dir)
        if target_output_dir.exists() and target_output_dir.is_file():
            raise ValueError(f"Output path '{target_output_dir}' is a file; please provide a directory path.")
        target_output_dir.mkdir(parents=True, exist_ok=True)
        output_dir = str(target_output_dir)

        if limit is not None and limit <= 0:
            raise ValueError("--limit must be a positive integer")

        selected_lambda = lambda_range.upper()
        fits_workers = workers if workers is not None else 1
        parquet_workers = workers if workers is not None else default_raw_parquet_workers()
        show_progress = bool(progress) if progress is not None else False

        # Datalink is intentionally kept on the existing FITS implementation.
        if use_datalink:
            if output_format != 'fits':
                click.echo("Note: --use-datalink is FITS-only; ignoring --output-format parquet.")
            sources = load_table(spectra_table)
            if limit is not None:
                sources = sources[:limit]
            compiler = SpectrumCompiler(max_extensions=max_extensions)
            if verbose:
                click.echo(f"Loaded {len(sources)} spectral sources from {spectra_table}")
                click.echo(f"Max extensions per file: {max_extensions}")
                click.echo(f"Output directory: {output_dir}")

            if idr_field.upper() != 'WIDE':
                click.echo(
                    "Note: --idr-field is ignored in --use-datalink mode "
                    "(datalink selection depends on --retrieval-type and --schema)."
                )
            if workers is not None and workers != 1:
                click.echo("Note: --workers is currently applied to local modes only; datalink runs with one worker.")

            lambda_to_retrieval = {
                'RGS': 'SPECTRA_RGS',
                'BGS': 'SPECTRA_BGS',
                'BOTH': 'ALL',
            }
            retrieval_to_lambda = {v: k for k, v in lambda_to_retrieval.items()}
            lambda_selected = selected_lambda
            retrieval_selected = retrieval_type.upper()

            lambda_src = ctx.get_parameter_source('lambda_range')
            retrieval_src = ctx.get_parameter_source('retrieval_type')
            lambda_explicit = lambda_src != click.core.ParameterSource.DEFAULT
            retrieval_explicit = retrieval_src != click.core.ParameterSource.DEFAULT

            if lambda_explicit:
                mapped = lambda_to_retrieval[lambda_selected]
                if retrieval_explicit and retrieval_selected != mapped:
                    click.echo(
                        "Warning: --lambda-range/-L and --retrieval-type disagree in datalink mode; "
                        f"using --lambda-range ({lambda_selected} -> {mapped})."
                    )
                retrieval_selected = mapped
            elif retrieval_explicit:
                lambda_selected = retrieval_to_lambda.get(retrieval_selected, 'RGS')
            else:
                retrieval_selected = lambda_to_retrieval[lambda_selected]

            click.echo(
                f"Datalink arm selection: lambda_range={lambda_selected}, "
                f"retrieval_type={retrieval_selected}"
            )
            source_id_col = compiler._resolve_source_id_column(sources, preferred='source_id')
            dedup_sources, dedup_stats = compiler.deduplicate_by_source_id(
                spectra_table=sources,
                source_id_col=source_id_col,
            )
            if dedup_stats['duplicate_rows_removed'] > 0:
                click.echo(
                    "Datalink deduplication: "
                    f"{dedup_stats['input_rows']} rows -> {dedup_stats['unique_sources']} unique {source_id_col} "
                    f"(removed {dedup_stats['duplicate_rows_removed']} duplicates)"
                )
            sources_for_datalink = dedup_sources

            archive = EuclidArchive(environment=environment)
            if credentials:
                archive.login(credentials_file=credentials)
            else:
                archive.login()
            if lambda_selected == 'BOTH':
                click.echo("Datalink BOTH mode: writing separate RGS and BGS compiled files.")
                output_files = []
                metadata_files = []

                rgs_prefix = f"{prefix}_rgs"
                rgs_files = compiler.compile_spectra_datalink(
                    spectra_table=sources_for_datalink,
                    euclid_client=archive.euclid,
                    output_dir=output_dir,
                    output_prefix=rgs_prefix,
                    retrieval_type='SPECTRA_RGS',
                    schema=schema,
                    overwrite=overwrite,
                )
                output_files.extend(rgs_files)
                metadata_files.append(
                    compiler.create_metadata_table(
                        spectra_table=sources_for_datalink,
                        output_files=rgs_files,
                        output_dir=output_dir,
                        output_name=f"{rgs_prefix}_metadata.fits",
                    )
                )

                bgs_prefix = f"{prefix}_bgs"
                bgs_files = compiler.compile_spectra_datalink(
                    spectra_table=sources_for_datalink,
                    euclid_client=archive.euclid,
                    output_dir=output_dir,
                    output_prefix=bgs_prefix,
                    retrieval_type='SPECTRA_BGS',
                    schema=schema,
                    overwrite=overwrite,
                )
                output_files.extend(bgs_files)
                metadata_files.append(
                    compiler.create_metadata_table(
                        spectra_table=sources_for_datalink,
                        output_files=bgs_files,
                        output_dir=output_dir,
                        output_name=f"{bgs_prefix}_metadata.fits",
                    )
                )
                total_processed = len(sources_for_datalink) * 2
            else:
                output_files = compiler.compile_spectra_datalink(
                    spectra_table=sources_for_datalink,
                    euclid_client=archive.euclid,
                    output_dir=output_dir,
                    output_prefix=prefix,
                    retrieval_type=retrieval_selected,
                    schema=schema,
                    overwrite=overwrite,
                )
                metadata_file = compiler.create_metadata_table(
                    spectra_table=sources_for_datalink,
                    output_files=output_files,
                    output_dir=output_dir,
                    output_name=f"{prefix}_metadata.fits"
                )
                metadata_files = [metadata_file]
                total_processed = len(sources_for_datalink)

            click.echo("Compilation completed successfully!")
            click.echo(f"Created {len(output_files)} FITS files:")
            for i, file_path in enumerate(output_files, 1):
                click.echo(f"  {i:2d}. {Path(file_path).name}")
            if len(metadata_files) == 1:
                click.echo(f"Metadata saved to: {Path(metadata_files[0]).name}")
            else:
                click.echo("Metadata files:")
                for i, meta_path in enumerate(metadata_files, 1):
                    click.echo(f"  {i:2d}. {Path(meta_path).name}")
            click.echo(f"Total spectra processed: {total_processed}")
            return

        if output_format == 'parquet':
            if verbose:
                click.echo(f"Output directory: {output_dir}")
                click.echo(f"Output format: parquet")
                click.echo(f"Parquet workers: {parquet_workers}")
                click.echo(f"Parquet chunk size: {chunk_size}")
                click.echo(f"LambdaRange selection: {selected_lambda}")

            if selected_lambda == 'BOTH':
                click.echo("Parquet BOTH mode: writing separate RGS and BGS parquet outputs.")
                stats_by_arm = []
                for arm in ('RGS', 'BGS'):
                    arm_prefix = str(target_output_dir / f"{prefix}_{arm.lower()}")
                    stats_by_arm.append((arm, spectra_to_parquet(
                        catalog_table=spectra_table,
                        output_prefix=arm_prefix,
                        chunk_size=chunk_size,
                        workers=parquet_workers,
                        lambda_range=arm,
                        overwrite=overwrite,
                        on_error=on_error,
                        show_progress=show_progress,
                        limit=limit,
                    )))
                total_exported = sum(stats.exported_rows for _, stats in stats_by_arm)
                total_requested = sum(stats.requested_rows for _, stats in stats_by_arm)
                click.echo("Parquet export completed successfully!")
                for arm, stats in stats_by_arm:
                    _echo_raw_parquet_stats(stats, label=f"{arm} parquet")
                click.echo(f"Total spectra exported: {total_exported}")
                click.echo(f"Total rows scanned: {total_requested}")
            else:
                stats = spectra_to_parquet(
                    catalog_table=spectra_table,
                    output_prefix=str(target_output_dir / prefix),
                    chunk_size=chunk_size,
                    workers=parquet_workers,
                    lambda_range=selected_lambda,
                    overwrite=overwrite,
                    on_error=on_error,
                    show_progress=show_progress,
                    limit=limit,
                )
                click.echo("Parquet export completed successfully!")
                _echo_raw_parquet_stats(stats, label="parquet")
                click.echo(f"Total spectra exported: {stats.exported_rows}")
            return

        # Legacy local FITS compatibility mode.
        sources = load_table(spectra_table)
        if limit is not None:
            sources = sources[:limit]
        compiler = SpectrumCompiler(max_extensions=max_extensions)
        if verbose:
            click.echo(f"Loaded {len(sources)} spectral sources from {spectra_table}")
            click.echo(f"Max extensions per file: {max_extensions}")
            click.echo(f"Output directory: {output_dir}")
            click.echo("Output format: fits")
            click.echo(f"Chunk workers: {fits_workers}")
            click.echo(f"LambdaRange selection: {selected_lambda}")
            if environment == 'IDR':
                click.echo(f"IDR field: {idr_field.upper()}")

        if target_output_dir.exists() and any(target_output_dir.iterdir()) and not overwrite:
            click.echo(
                f"Output directory '{target_output_dir}' is not empty. "
                "Resume mode: existing chunk files will be kept."
            )
        elif target_output_dir.exists() and any(target_output_dir.iterdir()) and overwrite:
            confirmed = click.confirm(
                f"Output directory '{target_output_dir}' is not empty. Overwrite existing files?",
                default=False,
            )
            if not confirmed:
                overwrite = False

        is_idr_deep = environment == 'IDR' and idr_field.upper() == 'DEEP'
        if not is_idr_deep and selected_lambda in {'BGS', 'BOTH'}:
            click.echo(
                f"Warning: {environment} {idr_field.upper()} only provides RGS in FITS mode; "
                f"falling back to RGS (requested {selected_lambda})."
            )
            selected_lambda = 'RGS'

        metadata_files = []
        if is_idr_deep:
            sources, xml_stats = compiler.annotate_lambda_range_from_xml(
                spectra_table=sources,
                datalabs_path_col='datalabs_path',
                file_name_col='file_name',
                lambda_col='lambda_range',
            )
            click.echo(
                "LambdaRange XML summary: "
                f"total={xml_stats['total']}, resolved={xml_stats['resolved']}, "
                f"RGS={xml_stats['rgs']}, BGS={xml_stats['bgs']}, "
                f"unresolved={xml_stats['unresolved']}, ambiguous={xml_stats['ambiguous']}"
            )

            if selected_lambda == 'BOTH':
                rgs_table, rgs_stats = compiler.filter_table_by_lambda_range(sources, lambda_range='RGS')
                bgs_table, bgs_stats = compiler.filter_table_by_lambda_range(sources, lambda_range='BGS')
                click.echo(
                    "LambdaRange selection summary (BOTH): "
                    f"RGS selected={rgs_stats['selected']}, BGS selected={bgs_stats['selected']}, "
                    f"unresolved={rgs_stats['unresolved']}, ambiguous={rgs_stats['ambiguous']}"
                )

                output_files = []
                total_processed = 0
                if len(rgs_table) > 0:
                    rgs_prefix = f"{prefix}_rgs"
                    rgs_files = compiler.compile_spectra(
                        spectra_table=rgs_table,
                        output_dir=output_dir,
                        output_prefix=rgs_prefix,
                        overwrite=overwrite,
                        workers=fits_workers,
                    )
                    output_files.extend(rgs_files)
                    metadata_files.append(
                        compiler.create_metadata_table(
                            spectra_table=rgs_table,
                            output_files=rgs_files,
                            output_dir=output_dir,
                            output_name=f"{rgs_prefix}_metadata.fits",
                        )
                    )
                    total_processed += len(rgs_table)
                if len(bgs_table) > 0:
                    bgs_prefix = f"{prefix}_bgs"
                    bgs_files = compiler.compile_spectra(
                        spectra_table=bgs_table,
                        output_dir=output_dir,
                        output_prefix=bgs_prefix,
                        overwrite=overwrite,
                        workers=fits_workers,
                    )
                    output_files.extend(bgs_files)
                    metadata_files.append(
                        compiler.create_metadata_table(
                            spectra_table=bgs_table,
                            output_files=bgs_files,
                            output_dir=output_dir,
                            output_name=f"{bgs_prefix}_metadata.fits",
                        )
                    )
                    total_processed += len(bgs_table)
            else:
                selected_sources, sel_stats = compiler.filter_table_by_lambda_range(
                    sources, lambda_range=selected_lambda
                )
                click.echo(
                    "LambdaRange selection summary: "
                    f"selected={sel_stats['selected']}, skipped_other_band={sel_stats['skipped_other_band']}, "
                    f"unresolved={sel_stats['unresolved']}, ambiguous={sel_stats['ambiguous']}"
                )
                output_files = compiler.compile_spectra(
                    spectra_table=selected_sources,
                    output_dir=output_dir,
                    output_prefix=prefix,
                    overwrite=overwrite,
                    workers=fits_workers,
                )
                metadata_files.append(
                    compiler.create_metadata_table(
                        spectra_table=selected_sources,
                        output_files=output_files,
                        output_dir=output_dir,
                        output_name=f"{prefix}_metadata.fits"
                    )
                )
                total_processed = len(selected_sources)
        else:
            output_files = compiler.compile_spectra(
                spectra_table=sources,
                output_dir=output_dir,
                output_prefix=prefix,
                overwrite=overwrite,
                workers=fits_workers,
            )
            metadata_files.append(
                compiler.create_metadata_table(
                    spectra_table=sources,
                    output_files=output_files,
                    output_dir=output_dir,
                    output_name=f"{prefix}_metadata.fits"
                )
            )
            total_processed = len(sources)

        click.echo("Compilation completed successfully!")
        click.echo(f"Created {len(output_files)} FITS files:")
        for i, file_path in enumerate(output_files, 1):
            click.echo(f"  {i:2d}. {Path(file_path).name}")
        if len(metadata_files) == 1:
            click.echo(f"Metadata saved to: {Path(metadata_files[0]).name}")
        else:
            click.echo("Metadata files:")
            for i, meta_path in enumerate(metadata_files, 1):
                click.echo(f"  {i:2d}. {Path(meta_path).name}")
        click.echo(f"Total spectra processed: {total_processed}")

    except Exception as e:
        click.echo(f"Error compiling spectra: {e}", err=True)
        sys.exit(1)
    finally:
        if archive is not None:
            archive.logout()


def _echo_raw_parquet_stats(stats, *, label: str) -> None:
    click.echo(
        f"{label}: requested_rows={stats.requested_rows} exported_rows={stats.exported_rows} "
        f"skipped_rows={stats.skipped_rows} failed_rows={stats.failed_rows} "
        f"file_missing_rows={getattr(stats, 'file_missing_rows', 0)} "
        f"parts={len(stats.output_files)} manifest={stats.manifest_path}"
    )
    if stats.failures_path:
        click.echo(f"failures={stats.failures_path}")


@click.command(name='dithers-to-parquet')
@click.option('--catalog-table', required=True, type=click.Path(exists=True),
              help='Catalog table with datalabs_path, file_name, and hdu_index')
@click.option('--output-prefix', required=True, type=click.Path(),
              help='Output prefix for combined/dither parquet parts')
@click.option('--chunk-size', type=int, default=2000, show_default=True,
              help='Catalog rows per processing chunk')
@click.option('--workers', type=int, default=None,
              help='Parallel FITS read workers (default: min(os.cpu_count(), 8))')
@click.option('--lambda-range', type=click.Choice(['BGS', 'RGS']), default=None,
              help='Optional LRANGE filter')
@click.option('--dithers-only', is_flag=True,
              help='Write only per-dither parquet parts, skipping combined rows')
@click.option('--overwrite', is_flag=True,
              help='Replace existing parquet parts for this prefix')
@click.option('--on-error', type=click.Choice(['fail', 'skip']), default='fail', show_default=True,
              help='Fail on bad rows or skip them and write a failures JSONL')
@click.option('--progress', 'progress', is_flag=True, default=None,
              help='Show export progress')
@click.option('--no-progress', 'progress', flag_value=False,
              help='Disable export progress')
def dithers_to_parquet_cli(catalog_table: str, output_prefix: str, chunk_size: int,
                           workers: Optional[int], lambda_range: Optional[str],
                           dithers_only: bool, overwrite: bool, on_error: str,
                           progress: Optional[bool]):
    """Export local Datalabs combined and per-dither SIR spectra to parquet."""
    from euclidkit.core.spectra_parquet import dithers_to_parquet

    try:
        show_progress = bool(progress) if progress is not None else False
        stats = dithers_to_parquet(
            catalog_table=catalog_table,
            output_prefix=output_prefix,
            chunk_size=chunk_size,
            workers=workers,
            lambda_range=lambda_range,
            include_combined=not dithers_only,
            overwrite=overwrite,
            on_error=on_error,
            show_progress=show_progress,
        )
        click.echo(
            "Dither parquet export completed successfully!\n"
            "objects_requested={objects_requested} objects_exported={objects_exported} "
            "combined_rows={combined_rows} dither_rows={dither_rows} skipped_rows={skipped_rows} "
            "failed_rows={failed_rows} file_missing_rows={file_missing_rows} "
            "combined_parts={combined_parts} dither_parts={dither_parts} manifest={manifest}".format(
                objects_requested=stats.objects_requested,
                objects_exported=stats.objects_exported,
                combined_rows=stats.combined_rows,
                dither_rows=stats.dither_rows,
                skipped_rows=stats.skipped_rows,
                failed_rows=stats.failed_rows,
                file_missing_rows=stats.file_missing_rows,
                combined_parts=len(stats.combined_output_files),
                dither_parts=len(stats.dither_output_files),
                manifest=stats.manifest_path,
            )
        )
        if stats.failures_path:
            click.echo(f"failures={stats.failures_path}")
    except Exception as e:
        click.echo(f"Error exporting dither spectra: {e}", err=True)
        sys.exit(1)


# Main command group for crossmatching functionality
@click.group()
def crossmatch_commands():
    """Crossmatching and spectral compilation commands."""
    pass


# Add commands to group
crossmatch_commands.add_command(crossmatch, name='crossmatch')
crossmatch_commands.add_command(query_spectra, name='query-spectra') 
crossmatch_commands.add_command(query_zspe, name='query-zspe')
crossmatch_commands.add_command(query_cutana, name='query-cutana')
crossmatch_commands.add_command(compile_spectra, name='compile-spectra')
crossmatch_commands.add_command(dithers_to_parquet_cli, name='dithers-to-parquet')


if __name__ == '__main__':
    crossmatch_commands()
