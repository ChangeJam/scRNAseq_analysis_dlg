"""
Centralized configuration loader.

All scripts should import get_config() from this module
instead of hardcoding paths or parameters.
"""

import yaml
from pathlib import Path
from typing import Any, Dict

def get_config(config_path:str = "config/config.yaml") -> Dict[str, Any]:
    """
    Load and return the analysis configuration from YAML file.

    Args:
        config_path: Relative or absolute path to config.yaml.

    Returns:
        Nested dictionary containing all pipeline parameters.

    Raises:
        FileNotFoundError: If config file does not exist.
    """

    config_file = Path(config_path)
    if not config_file.exists():
        raise FileNotFoundError(
            f"Configuration file not found: {config_file.resolve()}"
        )
    
    with open(config_file, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    return config