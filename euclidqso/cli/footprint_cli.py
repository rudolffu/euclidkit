"""
Command line interface for selecting catalog entries within Euclid footprints.
"""

import sys
from typing import Optional

import click


@click.command(name="select-footprint")
@click.option(
    "--input",
    "-i",
    required=True,
    type=click.Path(exists=True),
    help="Input catalog (FITS, CSV, VOTable, etc.).",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(),
    help="Output path for the filtered catalog.",
)
@click.option(
    "--survey",
    "-s",
    type=click.Choice(["WIDE", "DEEP"], case_sensitive=False),
    default="WIDE",
    show_default=True,
    help="Survey footprint to use.",
)
@click.option(
    "--data-release",
    "-d",
    type=click.Choice(["DR1"], case_sensitive=False),
    default="DR1",
    show_default=True,
    help="Data release version (default and currently available: DR1).",
)
@click.option(
    "--ra-col",
    default="ra",
    show_default=True,
    help="Right ascension column name (degrees).",
)
@click.option(
    "--dec-col",
    default="dec",
    show_default=True,
    help="Declination column name (degrees).",
)
@click.option(
    "--moc-path",
    type=click.Path(exists=True),
    help="Optional custom MOC file path (overrides packaged MOCs).",
)
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose output.")
def select_footprint(
    input: str,
    output: str,
    survey: str,
    data_release: str,
    ra_col: str,
    dec_col: str,
    moc_path: Optional[str],
    verbose: bool,
):
    """
    Extract catalog entries that fall within the Euclid MOC footprint.
    """
    import logging
    from euclidqso.core.footprint import filter_catalog_by_moc

    logging.basicConfig(level=logging.INFO if verbose else logging.WARNING)

    try:
        subset = filter_catalog_by_moc(
            input_catalog=input,
            output_catalog=output,
            survey=survey,
            data_release=data_release,
            ra_col=ra_col,
            dec_col=dec_col,
            moc_path=moc_path,
        )
        click.echo(
            f"Selected {len(subset)} sources within {survey.upper()} {data_release.upper()} footprint."
        )
        click.echo(f"Saved filtered catalog to: {output}")
    except Exception as e:
        click.echo(f"Error selecting by MOC footprint: {e}", err=True)
        if verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)
