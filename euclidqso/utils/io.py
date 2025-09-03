"""
Input/output utilities for euclidqso package.
"""

import os
import logging
from pathlib import Path
from typing import Union, Optional, List, Any

import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.io import fits

logger = logging.getLogger(__name__)


def ensure_dir(directory: Union[str, Path]) -> Path:
    """
    Ensure directory exists, create if necessary.
    
    Parameters
    ----------
    directory : str or Path
        Directory path
        
    Returns
    -------
    Path
        Directory path as Path object
    """
    dir_path = Path(directory)
    dir_path.mkdir(parents=True, exist_ok=True)
    return dir_path


def load_table(file_path: Union[str, Path], **kwargs) -> Table:
    """
    Load table from various formats.
    
    Parameters
    ----------
    file_path : str or Path
        Path to table file
    **kwargs
        Additional arguments passed to format-specific readers
        
    Returns
    -------
    Table
        Loaded table
    """
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    suffix = file_path.suffix.lower()
    
    try:
        if suffix == '.csv':
            # Load CSV via pandas for better handling, then convert
            df = pd.read_csv(file_path, **kwargs)
            table = Table.from_pandas(df)
            logger.info(f"Loaded CSV table with {len(table)} rows from {file_path}")
            
        elif suffix in ['.fits', '.fit']:
            table = Table.read(file_path, format='fits', **kwargs)
            logger.info(f"Loaded FITS table with {len(table)} rows from {file_path}")
            
        elif suffix in ['.vot', '.xml']:
            table = Table.read(file_path, format='votable', **kwargs)
            logger.info(f"Loaded VOTable with {len(table)} rows from {file_path}")
            
        elif suffix in ['.ecsv']:
            table = Table.read(file_path, format='ascii.ecsv', **kwargs)
            logger.info(f"Loaded ECSV table with {len(table)} rows from {file_path}")
            
        else:
            # Try astropy's automatic format detection
            table = Table.read(file_path, **kwargs)
            logger.info(f"Loaded table with {len(table)} rows from {file_path}")
            
    except Exception as e:
        logger.error(f"Failed to load table from {file_path}: {e}")
        raise
    
    return table


def save_table(
    table: Table, 
    file_path: Union[str, Path], 
    format: Optional[str] = None,
    overwrite: bool = True,
    **kwargs
) -> None:
    """
    Save table to file.
    
    Parameters
    ----------
    table : Table
        Table to save
    file_path : str or Path
        Output file path
    format : str, optional
        Output format (auto-detected from extension if not provided)
    overwrite : bool, default True
        Whether to overwrite existing files
    **kwargs
        Additional arguments passed to format-specific writers
    """
    file_path = Path(file_path)
    ensure_dir(file_path.parent)
    
    if format is None:
        suffix = file_path.suffix.lower()
        format_map = {
            '.csv': 'ascii.csv',
            '.fits': 'fits',
            '.fit': 'fits',
            '.vot': 'votable',
            '.xml': 'votable',
            '.ecsv': 'ascii.ecsv'
        }
        format = format_map.get(suffix, 'fits')  # Default to FITS
    
    try:
        table.write(file_path, format=format, overwrite=overwrite, **kwargs)
        logger.info(f"Saved table with {len(table)} rows to {file_path}")
        
    except Exception as e:
        logger.error(f"Failed to save table to {file_path}: {e}")
        raise


def load_spectra(input_path: Union[str, Path, List[str]]) -> List[str]:
    """
    Load list of spectrum files from various inputs.
    
    Parameters
    ----------
    input_path : str, Path, or list of str
        Input spectra (file, directory, or list of files)
        
    Returns
    -------
    List[str]
        List of spectrum file paths
    """
    if isinstance(input_path, list):
        return [str(p) for p in input_path]
    
    input_path = Path(input_path)
    
    if input_path.is_file():
        # Single file or file list
        if input_path.suffix.lower() in ['.fits', '.fit']:
            return [str(input_path)]
        else:
            # Text file with list of files
            with open(input_path, 'r') as f:
                files = [line.strip() for line in f if line.strip()]
            return files
    
    elif input_path.is_dir():
        # Directory with FITS files
        fits_files = list(input_path.glob('*.fits')) + list(input_path.glob('*.fit'))
        return [str(f) for f in sorted(fits_files)]
    
    else:
        raise ValueError(f"Input path does not exist: {input_path}")


def load_sources(sources_input: Union[str, Path, Table, pd.DataFrame]) -> Table:
    """
    Load source catalog from various inputs.
    
    Parameters
    ----------
    sources_input : str, Path, Table, or DataFrame
        Source catalog input
        
    Returns
    -------
    Table
        Source catalog as astropy Table
    """
    if isinstance(sources_input, Table):
        return sources_input.copy()
    elif isinstance(sources_input, pd.DataFrame):
        return Table.from_pandas(sources_input)
    else:
        return load_table(sources_input)


class DataLoader:
    """General data loader with caching and validation."""
    
    def __init__(self, cache_dir: Optional[Union[str, Path]] = None):
        """
        Initialize data loader.
        
        Parameters
        ----------
        cache_dir : str or Path, optional
            Cache directory for temporary files
        """
        if cache_dir is None:
            cache_dir = Path.home() / '.euclidqso' / 'cache'
        
        self.cache_dir = ensure_dir(cache_dir)
        logger.info(f"Initialized DataLoader with cache: {self.cache_dir}")
    
    def clear_cache(self):
        """Clear cache directory."""
        import shutil
        if self.cache_dir.exists():
            shutil.rmtree(self.cache_dir)
            self.cache_dir.mkdir()
            logger.info("Cleared data cache")


class FileManager:
    """File management utilities."""
    
    @staticmethod
    def get_file_info(file_path: Union[str, Path]) -> dict:
        """
        Get file information.
        
        Parameters
        ----------
        file_path : str or Path
            File path
            
        Returns
        -------
        dict
            File information
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            return {'exists': False}
        
        stat = file_path.stat()
        
        info = {
            'exists': True,
            'size_bytes': stat.st_size,
            'size_mb': stat.st_size / (1024 * 1024),
            'modified_time': stat.st_mtime,
            'is_file': file_path.is_file(),
            'is_dir': file_path.is_dir(),
            'suffix': file_path.suffix,
            'name': file_path.name
        }
        
        # Add FITS-specific information
        if file_path.suffix.lower() in ['.fits', '.fit'] and file_path.is_file():
            try:
                with fits.open(file_path) as hdul:
                    info['n_hdus'] = len(hdul)
                    info['primary_header_cards'] = len(hdul[0].header)
            except Exception as e:
                info['fits_error'] = str(e)
        
        return info
    
    @staticmethod
    def validate_paths(paths: List[Union[str, Path]]) -> dict:
        """
        Validate multiple file paths.
        
        Parameters
        ----------
        paths : list of str or Path
            File paths to validate
            
        Returns
        -------
        dict
            Validation results
        """
        results = {
            'valid': [],
            'invalid': [],
            'total_size_mb': 0
        }
        
        for path in paths:
            info = FileManager.get_file_info(path)
            if info['exists']:
                results['valid'].append(str(path))
                results['total_size_mb'] += info.get('size_mb', 0)
            else:
                results['invalid'].append(str(path))
        
        return results


class FormatConverter:
    """Format conversion utilities."""
    
    @staticmethod
    def table_to_csv(table: Table, output_path: Union[str, Path], **kwargs):
        """Convert astropy Table to CSV."""
        save_table(table, output_path, format='ascii.csv', **kwargs)
    
    @staticmethod
    def csv_to_fits(csv_path: Union[str, Path], fits_path: Union[str, Path], **kwargs):
        """Convert CSV to FITS table."""
        table = load_table(csv_path)
        save_table(table, fits_path, format='fits', **kwargs)
    
    @staticmethod
    def pandas_to_table(df: pd.DataFrame) -> Table:
        """Convert pandas DataFrame to astropy Table."""
        return Table.from_pandas(df)
    
    @staticmethod
    def table_to_pandas(table: Table) -> pd.DataFrame:
        """Convert astropy Table to pandas DataFrame."""
        return table.to_pandas()