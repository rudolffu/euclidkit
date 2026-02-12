"""
Configuration management for euclidkit.

Handles default settings, user configuration files, and environment variables.
"""

import os
import yaml
from pathlib import Path
from typing import Dict, Any, Optional, Union


class EuclidKitConfig:
    """Configuration manager for euclidkit package."""
    
    def __init__(self):
        self.config_dir = Path.home() / '.euclidkit'
        self.config_file = self.config_dir / 'config.yaml'
        self._config = self._load_default_config()
        self._load_user_config()
        self._load_env_overrides()
    
    def _load_default_config(self) -> Dict[str, Any]:
        """Load default configuration."""
        return {
            # Data access settings
            'data': {
                'euclid_volumes': {
                    'q1': '/data/euclid_q1/',
                    'ero': '/data/euc_ero_data_01/',
                },
                'credentials_file': '/media/user/cred.txt',
                'cache_dir': str(Path.home() / '.euclidkit' / 'cache'),
            },
            # Analysis settings
            'analysis': {
                'redshift': {
                    'method': 'xcorr',
                    'z_range': [0.0, 6.0],
                    'template_dir': None,
                },
                'spectroscopy': {
                    'snr_min': 3.0,
                    'wavelength_range': [3500, 25000],  # Angstroms
                },
                'photometry': {
                    'mag_systems': ['AB', 'Vega'],
                    'error_floor': 0.01,
                },
            },
            # Processing settings
            'processing': {
                'parallel': True,
                'n_cores': None,  # Use all available
                'chunk_size': 100,
                'memory_limit': '8GB',
            },
            # Output settings
            'output': {
                'format': 'fits',
                'compression': True,
                'plot_format': 'png',
                'plot_dpi': 150,
            },
            # External survey integration
            'external': {
                'desi': {
                    'enabled': False,
                    'data_release': 'DESI-DR1',
                    'search_radius': 2.0,  # arcseconds
                },
                'wise': {
                    'enabled': True,
                    'catalog': 'AllWISE',
                },
                'gaia': {
                    'enabled': True,
                    'data_release': 'DR3',
                },
            },
            # Logging
            'logging': {
                'level': 'INFO',
                'file': None,
                'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            },
        }
    
    def _load_user_config(self):
        """Load user configuration file if it exists."""
        if self.config_file.exists():
            try:
                with open(self.config_file, 'r') as f:
                    user_config = yaml.safe_load(f)
                    self._merge_config(self._config, user_config)
            except Exception as e:
                print(f"Warning: Could not load user config: {e}")
    
    def _load_env_overrides(self):
        """Load configuration overrides from environment variables."""
        env_mappings = {
            'EUCLIDKIT_CREDENTIALS': 'data.credentials_file',
            'EUCLIDKIT_CACHE_DIR': 'data.cache_dir',
            'EUCLIDKIT_LOG_LEVEL': 'logging.level',
            'EUCLIDKIT_PARALLEL': 'processing.parallel',
            'EUCLIDKIT_N_CORES': 'processing.n_cores',
        }
        
        for env_var, config_path in env_mappings.items():
            if env_var in os.environ:
                self._set_nested_config(config_path, os.environ[env_var])
    
    def _merge_config(self, base: Dict, override: Dict):
        """Recursively merge configuration dictionaries."""
        for key, value in override.items():
            if key in base and isinstance(base[key], dict) and isinstance(value, dict):
                self._merge_config(base[key], value)
            else:
                base[key] = value
    
    def _set_nested_config(self, path: str, value: Any):
        """Set a nested configuration value using dot notation."""
        keys = path.split('.')
        config = self._config
        for key in keys[:-1]:
            config = config.setdefault(key, {})
        config[keys[-1]] = value
    
    def get(self, key: str, default: Any = None) -> Any:
        """Get a configuration value using dot notation."""
        keys = key.split('.')
        config = self._config
        for k in keys:
            if isinstance(config, dict) and k in config:
                config = config[k]
            else:
                return default
        return config
    
    def set(self, key: str, value: Any):
        """Set a configuration value using dot notation."""
        self._set_nested_config(key, value)
    
    def save_user_config(self):
        """Save current configuration to user config file."""
        self.config_dir.mkdir(parents=True, exist_ok=True)
        with open(self.config_file, 'w') as f:
            yaml.dump(self._config, f, default_flow_style=False, sort_keys=False)


def generate_config_template(template_type: str = 'basic') -> Dict[str, Any]:
    """Generate configuration template."""
    config = EuclidKitConfig()
    
    if template_type == 'basic':
        return {
            'data': config.get('data'),
            'analysis': {
                'redshift': config.get('analysis.redshift'),
                'spectroscopy': config.get('analysis.spectroscopy'),
            },
            'output': config.get('output'),
        }
    elif template_type == 'advanced':
        return config._config
    elif template_type == 'pipeline':
        return {
            'input': 'sources.csv',
            'output': 'results/',
            'steps': ['cutouts', 'spectra', 'photometry', 'validation', 'composites'],
            'cutouts': {
                'instrument': 'both',
                'size': '1arcmin',
                'format': 'fits',
            },
            'spectra': {
                'method': 'xcorr',
                'z_range': [0.0, 6.0],
                'snr_min': 3.0,
            },
            'composites': {
                'z_bins': [[0.0, 2.0], [2.0, 4.0], [4.0, 6.0]],
                'weighting': 'snr',
                'bootstrap': 100,
            },
        }
    else:
        raise ValueError(f"Unknown template type: {template_type}")


# Global configuration instance
config = EuclidKitConfig()
