"""Reference file manager for downloading and caching remote reference files.

Uses the pooch library for reliable file fetching with caching support.
"""

import logging
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional

import pooch
import yaml

logger = logging.getLogger(__name__)


class ReferenceFileManager:
    """Manager for downloading, caching, and clearing reference files using pooch."""

    def __init__(self, cache_dir: str, config_path: str):
        self.cache_dir = Path(cache_dir)
        with open(config_path) as f:
            self.config: Dict[str, Any] = yaml.safe_load(f)

    def fetch(self, category: str, key: str) -> List[str]:
        """Fetch and extract a reference file, returning paths to extracted files.

        Uses pooch's built-in caching - only downloads if not already cached.
        """
        file_config = self._get_file_config(category, key)
        url: str = file_config["url"]

        # Choose processor based on file type
        if url.endswith(".tar.gz"):
            processor = pooch.Untar()
        elif url.endswith(".zip"):
            processor = pooch.Unzip()
        else:
            processor = None

        # pooch handles caching automatically
        result = pooch.retrieve(
            url=url,
            known_hash=None,
            path=str(self.cache_dir / category),
            processor=processor,
        )
        # Normalize to always return a list
        return result if isinstance(result, list) else [result]

    def _get_file_config(self, category: str, key: str) -> Dict[str, Any]:
        """Get configuration for a specific reference file."""
        if category not in self.config:
            raise KeyError(f"Category '{category}' not found")
        if key not in self.config[category]:
            raise KeyError(f"Key '{key}' not found in '{category}'")
        result: Dict[str, Any] = self.config[category][key]
        return result

    def get_key_by_organism_id(self, category: str, organism_id: str) -> Optional[str]:
        """Look up a file key by organism ID."""
        for key, cfg in self.config.get(category, {}).items():
            if cfg.get("organism_id") == organism_id:
                return str(key)
        return None

    def clear_cache(self, category: Optional[str] = None) -> List[Path]:
        """Clear cached reference files."""
        cleared: List[Path] = []
        target = self.cache_dir / category if category else self.cache_dir

        if target.exists():
            for item in target.iterdir():
                if item.name.startswith("."):
                    continue
                if item.is_dir():
                    shutil.rmtree(item)
                else:
                    item.unlink()
                cleared.append(item)

        return cleared

    def list_available_files(self) -> Dict[str, Any]:
        """List all available reference files."""
        return self.config
