"""
Spectral data handling module for euclidkit package.

Provides functionality for loading, processing, and combining spectral data.
"""

import os
import logging
from pathlib import Path
from typing import Union, Optional, List, Dict, Any, TYPE_CHECKING

if TYPE_CHECKING:
    from astropy.io.fits.hdu.base import ExtensionHDU
import warnings

import numpy as np
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm

from euclidkit.utils.io import ensure_dir
from euclidkit.config import config

logger = logging.getLogger(__name__)


class SpectrumLoader:
    """Loader for individual Euclid spectra."""
    
    def __init__(self):
        self.data_volumes = config.get('data.euclid_volumes', {})
    
    def load_spectrum(self, file_path: str, hdu_index: int = 1) -> 'ExtensionHDU':
        """
        Load a single spectrum from FITS file.
        
        Parameters
        ----------
        file_path : str
            Path to FITS file
        hdu_index : int
            HDU index containing the spectrum
            
        Returns
        -------
        ExtensionHDU
            Spectrum HDU
        """
        try:
            with fits.open(file_path) as hdul:
                if hdu_index >= len(hdul):
                    raise ValueError(f"HDU index {hdu_index} out of range for {file_path}")
                
                # Create a copy to avoid issues with closed file
                spectrum_hdu = hdul[hdu_index].copy()
                
            return spectrum_hdu
            
        except Exception as e:
            logger.error(f"Failed to load spectrum from {file_path}[{hdu_index}]: {e}")
            raise
    
    def verify_spectrum_path(self, datalabs_path: str, file_name: str) -> str:
        """
        Verify and construct full spectrum file path.
        
        Parameters
        ----------
        datalabs_path : str
            Datalabs path from archive query
        file_name : str
            FITS file name
            
        Returns
        -------
        str
            Full file path
        """
        full_path = os.path.join(datalabs_path, file_name)
        
        if not os.path.exists(full_path):
            logger.warning(f"Spectrum file not found: {full_path}")
            return None
        
        return full_path


class SpectrumCompiler:
    """Compiler for creating multi-extension FITS files from individual spectra."""
    
    def __init__(self, max_extensions: int = 5000):
        """
        Initialize spectrum compiler.
        
        Parameters
        ----------
        max_extensions : int, default 5000
            Maximum number of extensions per output file
        """
        self.max_extensions = max_extensions
        self.loader = SpectrumLoader()
    
    def compile_spectra(
        self,
        spectra_table: Table,
        output_dir: Union[str, Path],
        output_prefix: str = "compiled_spectra",
        source_id_col: str = 'source_id',
        datalabs_path_col: str = 'datalabs_path',
        file_name_col: str = 'file_name', 
        hdu_index_col: str = 'hdu_index',
        overwrite: bool = True
    ) -> List[str]:
        """
        Compile spectra into chunked FITS files.
        
        Parameters
        ----------
        spectra_table : Table
            Table with spectral source information
        output_dir : str or Path
            Output directory for compiled FITS files
        output_prefix : str, default "compiled_spectra"
            Prefix for output file names
        source_id_col : str, default 'source_id'
            Column name for source IDs
        datalabs_path_col : str, default 'datalabs_path'
            Column name for datalabs paths
        file_name_col : str, default 'file_name'
            Column name for file names
        hdu_index_col : str, default 'hdu_index'
            Column name for HDU indices
        overwrite : bool, default True
            Whether to overwrite existing files
            
        Returns
        -------
        List[str]
            List of output file paths
        """
        output_dir = Path(output_dir)
        ensure_dir(output_dir)
        
        # Validate input table
        required_cols = [source_id_col, datalabs_path_col, file_name_col, hdu_index_col]
        missing_cols = [col for col in required_cols if col not in spectra_table.colnames]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        logger.info(f"Compiling {len(spectra_table)} spectra into FITS files")
        logger.info(f"Max extensions per file: {self.max_extensions}")
        
        output_files = []
        chunk_number = 1
        failed_count = 0
        
        # Process spectra in chunks
        for i in range(0, len(spectra_table), self.max_extensions):
            chunk = spectra_table[i:i + self.max_extensions]
            
            # Create output file for this chunk
            output_file = output_dir / f"{output_prefix}_chunk_{chunk_number:03d}.fits"
            output_files.append(str(output_file))
            
            logger.info(f"Creating chunk {chunk_number}: {len(chunk)} spectra -> {output_file.name}")
            
            # Create primary HDU
            primary_hdu = fits.PrimaryHDU()
            primary_hdu.header['COMMENT'] = f'Compiled Euclid spectra - Chunk {chunk_number}'
            primary_hdu.header['NSPECTRA'] = len(chunk)
            primary_hdu.header['CREATOR'] = 'euclidkit.core.spectra.SpectrumCompiler'
            
            hdul_new = fits.HDUList([primary_hdu])
            
            # Add each spectrum as an extension
            chunk_failed = 0
            for j, source in enumerate(tqdm(chunk, desc=f"Chunk {chunk_number}")):
                try:
                    # Construct file path
                    file_path = self.loader.verify_spectrum_path(
                        source[datalabs_path_col], 
                        source[file_name_col]
                    )
                    
                    if file_path is None:
                        chunk_failed += 1
                        continue
                    
                    # Load spectrum
                    spectrum_hdu = self.loader.load_spectrum(
                        file_path, 
                        source[hdu_index_col]
                    )
                    
                    # Set extension name to source ID
                    spectrum_hdu.name = str(source[source_id_col])
                    
                    # Add metadata to header
                    spectrum_hdu.header['SOURCE_ID'] = source[source_id_col]
                    
                    if 'ra_obj' in source.colnames:
                        spectrum_hdu.header['RA'] = source['ra_obj']
                    if 'dec_obj' in source.colnames:
                        spectrum_hdu.header['DEC'] = source['dec_obj']

                    hdul_new.append(spectrum_hdu)
                    
                except Exception as e:
                    logger.warning(f"Failed to add spectrum {source[source_id_col]}: {e}")
                    chunk_failed += 1
                    continue
            
            # Write chunk file
            try:
                hdul_new.writeto(output_file, overwrite=overwrite)
                logger.info(f"Wrote {len(hdul_new) - 1} spectra to {output_file.name}")
                
                if chunk_failed > 0:
                    logger.warning(f"Failed to include {chunk_failed} spectra in chunk {chunk_number}")
                    failed_count += chunk_failed
                
            except Exception as e:
                logger.error(f"Failed to write chunk {chunk_number}: {e}")
                if output_file.exists():
                    output_file.unlink()
                output_files.pop()  # Remove failed file from list
            
            finally:
                hdul_new.close()
            
            chunk_number += 1
        
        logger.info(f"Compilation complete: {len(output_files)} files created")
        if failed_count > 0:
            logger.warning(f"Total failed spectra: {failed_count}")
        
        return output_files
    
    def create_metadata_table(
        self,
        spectra_table: Table,
        output_files: List[str],
        output_dir: Union[str, Path],
        output_name: str = "spectra_metadata.fits"
    ) -> str:
        """
        Create metadata table for compiled spectra.
        
        Parameters
        ----------
        spectra_table : Table
            Original spectral source table
        output_files : List[str]
            List of output FITS files
        output_dir : str or Path
            Output directory
        output_name : str, default "spectra_metadata.fits"
            Metadata file name
            
        Returns
        -------
        str
            Path to metadata file
        """
        output_path = Path(output_dir) / output_name
        
        # Create extended metadata table
        metadata = spectra_table.copy()
        
        # Add chunk information
        chunk_info = []
        for i, source in enumerate(spectra_table):
            chunk_num = (i // self.max_extensions) + 1
            extension_num = (i % self.max_extensions) + 1  # +1 for primary HDU
            chunk_file = output_files[chunk_num - 1] if chunk_num <= len(output_files) else None
            
            chunk_info.append({
                'chunk_number': chunk_num,
                'chunk_file': Path(chunk_file).name if chunk_file else None,
                'extension_number': extension_num
            })
        
        # Add new columns
        metadata['chunk_number'] = [info['chunk_number'] for info in chunk_info]
        metadata['chunk_file'] = [info['chunk_file'] for info in chunk_info]
        metadata['extension_number'] = [info['extension_number'] for info in chunk_info]
        
        # Save metadata
        metadata.write(output_path, format='fits', overwrite=True)
        logger.info(f"Saved metadata table to {output_path}")
        
        return str(output_path)
    
    def compile_single_fits(
        self,
        spectra_table: Table,
        output_file: Union[str, Path],
        source_id_col: str = 'source_id',
        datalabs_path_col: str = 'datalabs_path',
        file_name_col: str = 'file_name',
        hdu_index_col: str = 'hdu_index',
        overwrite: bool = True
    ) -> str:
        """
        Compile spectra into a single FITS file.
        
        This is a convenience method for creating a single FITS file
        containing all spectra, similar to the notebook cell 23 pattern.
        
        Parameters
        ----------
        spectra_table : Table
            Table with spectral source information
        output_file : str or Path
            Output FITS file path
        source_id_col : str, default 'source_id'
            Column name for source IDs
        datalabs_path_col : str, default 'datalabs_path'
            Column name for datalabs paths
        file_name_col : str, default 'file_name'
            Column name for file names
        hdu_index_col : str, default 'hdu_index'
            Column name for HDU indices
        overwrite : bool, default True
            Whether to overwrite existing files
            
        Returns
        -------
        str
            Path to output file
        """
        output_path = Path(output_file)
        output_dir = output_path.parent
        output_stem = output_path.stem
        
        # Temporarily set max_extensions to accommodate all spectra
        original_max = self.max_extensions
        self.max_extensions = len(spectra_table) + 1000
        
        try:
            output_files = self.compile_spectra(
                spectra_table,
                output_dir=output_dir,
                output_prefix=output_stem,
                source_id_col=source_id_col,
                datalabs_path_col=datalabs_path_col,
                file_name_col=file_name_col,
                hdu_index_col=hdu_index_col,
                overwrite=overwrite
            )
            
            if len(output_files) == 1:
                # Rename the chunked file to the requested name
                chunked_file = Path(output_files[0])
                if chunked_file != output_path:
                    chunked_file.rename(output_path)
                return str(output_path)
            else:
                logger.warning(f"Expected single file but got {len(output_files)} files")
                return output_files[0] if output_files else None
                
        finally:
            # Restore original max_extensions
            self.max_extensions = original_max


class SpectrumProcessor:
    """Processor for spectral analysis and manipulation."""
    
    def __init__(self):
        pass
    
    def validate_spectrum(self, spectrum_hdu: 'ExtensionHDU') -> Dict[str, Any]:
        """
        Validate a spectrum HDU.
        
        Parameters
        ----------
        spectrum_hdu : ExtensionHDU
            Spectrum HDU to validate
            
        Returns
        -------
        Dict[str, Any]
            Validation results
        """
        validation = {
            'valid': True,
            'issues': [],
            'n_pixels': 0,
            'wavelength_range': None,
            'flux_range': None
        }
        
        try:
            data = spectrum_hdu.data
            if data is None:
                validation['valid'] = False
                validation['issues'].append('No data in HDU')
                return validation
            
            # Check for required columns
            required_cols = ['WAVELENGTH', 'SIGNAL']  # Standard Euclid spectrum columns
            missing_cols = []
            
            if hasattr(data, 'columns'):  # Binary table
                existing_cols = [col.name for col in data.columns]
                missing_cols = [col for col in required_cols if col not in existing_cols]
                if missing_cols:
                    validation['issues'].append(f'Missing columns: {missing_cols}')
                else:
                    validation['n_pixels'] = len(data)
                    validation['wavelength_range'] = (
                        float(np.min(data['WAVELENGTH'])), 
                        float(np.max(data['WAVELENGTH']))
                    )
                    validation['flux_range'] = (
                        float(np.min(data['SIGNAL'])),
                        float(np.max(data['SIGNAL']))
                    )
            else:
                validation['issues'].append('Unexpected data format')
            
            if validation['issues']:
                validation['valid'] = False
                
        except Exception as e:
            validation['valid'] = False
            validation['issues'].append(f'Validation error: {str(e)}')
        
        return validation