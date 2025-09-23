#!/usr/bin/env python3
"""
Data Manager Utility
Handles data loading, processing, and management for all data types
"""

import streamlit as st
import pandas as pd
import numpy as np
import anndata
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

class DataManager:
    """Manages data loading and processing for all three data types"""
    
    def __init__(self, data_dir: str = None):
        # Set up different data directories for different data types
        base_dir = "/Users/mayongzhi/Desktop/FS_lab/dashboard/test_data"
        self.scrna_data_dir = Path(base_dir) / "scRNA-seq"
        self.spatial_data_dir = Path(base_dir)  # Spatial data is in the main test_data directory
        self.tangram_data_dir = Path(base_dir)  # Tangram data is also in the main test_data directory
        
        # Use provided data_dir if specified (for backward compatibility)
        if data_dir:
            self.data_dir = Path(data_dir)
        else:
            self.data_dir = self.spatial_data_dir  # Default to spatial data directory
        
        self.scrna_data = {}  # Dictionary to store scRNA-seq data by sample
        self.registered_data: Dict[str, anndata.AnnData] = {}
        self.finalized_data: Dict[str, anndata.AnnData] = {}
        self._registered_catalog: Optional[Dict[str, Dict]] = None
        self._registered_alias_map: Dict[str, str] = {}
        self._dapi_cache: Dict[str, Optional[Path]] = {}

    # --- Shared Helpers ---

    def _infer_spatial_metadata(self, file_path: Path) -> Tuple[str, str, List[str], str]:
        """Infer spatial sample identifiers and metadata from a file path."""
        name = file_path.stem
        lower_name = name.lower()

        sample_pattern = re.compile(r"e\d+(?:\.\d+)?(?:_\d+)?", re.IGNORECASE)
        match = sample_pattern.search(lower_name)

        if match:
            primary_sample = match.group(0).upper()
            decimal_match = re.match(r"(E\d+\.\d+)", primary_sample)
            base_match = re.match(r"(E\d+)", primary_sample)

            aliases: List[str] = []
            if decimal_match and decimal_match.group(0) != primary_sample:
                aliases.append(decimal_match.group(0))
            if base_match and base_match.group(0) not in aliases and base_match.group(0) != primary_sample:
                aliases.append(base_match.group(0))

            display_suffix: List[str] = []
            if base_match and base_match.group(0) != primary_sample:
                display_suffix.append(f"stage {base_match.group(0)}")
        else:
            primary_sample = name
            aliases = []
            display_suffix = []

        file_name_lower = file_path.name.lower()
        if 'tangram' in file_name_lower:
            data_type = 'tangram'
        elif 'predicted' in file_name_lower:
            data_type = 'predicted'
        elif 'spatial' in file_name_lower:
            data_type = 'spatial'
        else:
            data_type = 'unknown'

        if data_type not in ('spatial', 'unknown'):
            display_suffix.append(data_type.replace('_', ' ').title())

        if display_suffix and match:
            display_name = f"{primary_sample} ({', '.join(display_suffix)})"
        elif display_suffix:
            display_name = f"{name} ({', '.join(display_suffix)})"
        else:
            display_name = primary_sample if match else name

        aliases = sorted({alias for alias in aliases if alias})
        return primary_sample, display_name, aliases, data_type

    def _finalize_spatial_dataset(self, adata: anndata.AnnData) -> anndata.AnnData:
        """Ensure spatial coordinates and basic metrics are present."""
        if 'spatial' in adata.obsm:
            spatial_coords = adata.obsm['spatial']
            adata.obs['spatial_x'] = spatial_coords[:, 0]
            adata.obs['spatial_y'] = spatial_coords[:, 1]
            adata.obsm['spatial'] = spatial_coords
        elif 'X_spatial' in adata.obsm:
            spatial_coords = adata.obsm['X_spatial']
            adata.obs['spatial_x'] = spatial_coords[:, 0]
            adata.obs['spatial_y'] = spatial_coords[:, 1]
            adata.obsm['spatial'] = spatial_coords
        elif 'spatial_x' in adata.obs.columns and 'spatial_y' in adata.obs.columns:
            spatial_coords = np.column_stack([adata.obs['spatial_x'], adata.obs['spatial_y']])
            adata.obsm['spatial'] = spatial_coords

        # Basic QC metrics for convenience in downstream plots
        adata.obs['total_counts'] = np.sum(adata.X, axis=1)
        adata.obs['n_genes_by_counts'] = np.sum(adata.X > 0, axis=1)
        return adata

    def _register_catalog_aliases(self, canonical_key: str, entry: Dict, extra_alias: Optional[str] = None):
        """Record aliases for a dataset across spatial and Tangram contexts."""
        if not canonical_key or entry is None:
            return

        alias_values = set(entry.get('aliases', []))
        alias_values.update({
            canonical_key,
            entry.get('primary_sample'),
            entry.get('file_stem'),
            entry.get('file_name'),
            extra_alias,
        })
        alias_values.discard(None)

        sorted_aliases = sorted(alias_values)
        entry['aliases'] = sorted_aliases

        for alias in sorted_aliases:
            lower_alias = alias.lower()
            self._registered_alias_map[lower_alias] = canonical_key

    # --- scRNA-seq (Page 2) ---

    def load_scrna_data(self, sample: str) -> anndata.AnnData:
        """Load scRNA-seq data for specific sample"""
        try:
            # Specific path for Cartana.h5ad file
            cartana_path = self.scrna_data_dir / "Cartana_simplified.h5ad"
            
            # Check if Cartana.h5ad exists and load it
            if cartana_path.exists():
                file_path = str(cartana_path)
                adata = anndata.read_h5ad(file_path)
                
                # Filter data based on selected sample (orig.ident)
                if 'orig.ident' in adata.obs.columns:
                    # Filter to only include cells from the selected sample
                    sample_mask = adata.obs['orig.ident'] == sample
                    adata = adata[sample_mask].copy()
                    
                    if adata.n_obs == 0:
                        st.error(f"No cells found for sample '{sample}' in orig.ident")
                        return None
                
                # Calculate quality metrics
                adata.obs['total_counts'] = np.sum(adata.X, axis=1)
                adata.obs['n_genes_by_counts'] = np.sum(adata.X > 0, axis=1)
                adata.var['total_counts'] = np.sum(adata.X, axis=0)
                adata.var['n_cells_by_counts'] = np.sum(adata.X > 0, axis=0)
                
                self.scrna_data[sample] = adata
                return adata
            else:
                # Fallback to original logic for other samples
                sample_files = list(self.scrna_data_dir.glob(f"*{sample}*scrna*.h5ad"))
                if not sample_files:
                    sample_files = list(self.scrna_data_dir.glob(f"*{sample}*.h5ad"))
                
                if sample_files:
                    file_path = str(sample_files[0])
                    adata = anndata.read_h5ad(file_path)
                    
                    # Calculate quality metrics
                    adata.obs['total_counts'] = np.sum(adata.X, axis=1)
                    adata.obs['n_genes_by_counts'] = np.sum(adata.X > 0, axis=1)
                    adata.var['total_counts'] = np.sum(adata.X, axis=0)
                    adata.var['n_cells_by_counts'] = np.sum(adata.X > 0, axis=0)
                    
                    self.scrna_data[sample] = adata
                    return adata
                else:
                    st.error(f"No scRNA-seq data found for sample {sample}")
                    return None
        except Exception as e:
            st.error(f"Error loading scRNA-seq data for {sample}: {e}")
            return None

    def get_scrna_available_samples(self) -> List[str]:
        """Get list of available samples (E12, E14, E17) for scRNA-seq data"""
        samples = []
        for file_path in self.scrna_data_dir.glob("*.h5ad"):
            file_name = file_path.name.lower()
            if 'e12' in file_name and 'e12' not in samples:
                samples.append('E12')
            elif 'e14' in file_name and 'e14' not in samples:
                samples.append('E14')
            elif 'e17' in file_name and 'e17' not in samples:
                samples.append('E17')
        
        # If no samples found, return default list
        if not samples:
            samples = ['E12', 'E14', 'E17']  # Default options
        
        return sorted(samples)

    def get_scrna_sample_options(self) -> List[str]:
        """Get unique orig.ident values from Cartana.h5ad file for scRNA-seq sample selection"""
        try:
            cartana_path = self.scrna_data_dir / "Cartana.h5ad"
            
            if cartana_path.exists():
                # Load the Cartana.h5ad file temporarily to get orig.ident values
                adata = anndata.read_h5ad(cartana_path)
                
                # Check if orig.ident column exists
                if 'orig.ident' in adata.obs.columns:
                    unique_identities = sorted(adata.obs['orig.ident'].unique().tolist())
                    return unique_identities
                else:
                    # If orig.ident doesn't exist, return default options
                    return ['E12', 'E14', 'E17']
            else:
                # Fallback to default options if file doesn't exist
                return ['E12', 'E14', 'E17']
        except Exception as e:
            print(f"Error reading Cartana.h5ad for sample options: {e}")
            return ['E12', 'E14', 'E17']

    # --- Registered Data (Spatial & Tangram Pages) ---

    def _build_registered_catalog(self) -> Dict[str, Dict]:
        """Discover available registered datasets."""
        catalog: Dict[str, Dict] = {}

        candidate_paths = sorted(self.tangram_data_dir.rglob("*.h5ad"))

        for file_path in candidate_paths:
            if not file_path.is_file():
                continue

            if any(part.lower() == "scrna-seq" for part in file_path.parts):
                continue

            if "tangram" not in file_path.name.lower():
                continue

            primary_sample, display_name, aliases, _ = self._infer_spatial_metadata(file_path)
            data_type = "registered"

            key = primary_sample or file_path.stem
            base_key = key
            counter = 2
            while key in catalog:
                key = f"{base_key}_{counter}"
                counter += 1

            alias_values = set(aliases)
            alias_values.update(filter(None, [primary_sample, file_path.stem, file_path.name, key]))

            entry = {
                "key": key,
                "primary_sample": primary_sample,
                "display_name": display_name,
                "aliases": sorted(alias_values),
                "data_type": data_type,
                "path": file_path,
                "file_name": file_path.name,
                "file_stem": file_path.stem,
            }
            catalog[key] = entry

        return catalog

    def refresh_registered_catalog(self):
        """Clear the cached registered catalog."""
        self._registered_catalog = None
        self._registered_alias_map.clear()
        self._dapi_cache.clear()

    def get_registered_catalog(self) -> Dict[str, Dict]:
        """Return the cached registered catalog, rebuilding if needed."""
        if self._registered_catalog is None:
            self._registered_catalog = self._build_registered_catalog()
        return self._registered_catalog

    def _resolve_registered_entry(self, sample: str) -> Optional[Dict]:
        """Resolve a sample identifier to a registered catalog entry."""
        if not sample:
            return None

        catalog = self.get_registered_catalog()

        if sample in catalog:
            return catalog[sample]

        sample_lower = sample.lower()

        alias_key = self._registered_alias_map.get(sample_lower)
        if alias_key and alias_key in catalog:
            return catalog[alias_key]

        for entry in catalog.values():
            candidate_aliases = [
                entry.get("key"),
                entry.get("primary_sample"),
                entry.get("file_stem"),
                entry.get("file_name"),
            ] + entry.get("aliases", [])

            for alias in candidate_aliases:
                if alias and alias.lower() == sample_lower:
                    return entry

        return None

    def _resolve_registered_key(self, sample: str, *, finalized: bool = False) -> Optional[str]:
        """Resolve a sample identifier to an in-memory dataset key."""
        if not sample:
            return None

        data_map = self.finalized_data if finalized else self.registered_data
        if sample in data_map:
            return sample

        alias = self._registered_alias_map.get(sample.lower())
        if alias and alias in data_map:
            return alias

        return alias

    def get_registered_dataset(self, sample: str, *, finalized: bool = False) -> Optional[anndata.AnnData]:
        """Retrieve a loaded registered dataset by sample or alias."""
        key = self._resolve_registered_key(sample, finalized=finalized)
        if not key:
            return None
        data_map = self.finalized_data if finalized else self.registered_data
        return data_map.get(key)

    def load_registered_data(self, sample: str, *, finalize: bool = False) -> Optional[anndata.AnnData]:
        """Load a registered dataset from disk, optionally finalising spatial coordinates."""
        try:
            existing = self.get_registered_dataset(sample, finalized=finalize)
            if existing is not None:
                return existing

            entry = self._resolve_registered_entry(sample)
            if entry is None:
                self.refresh_registered_catalog()
                entry = self._resolve_registered_entry(sample)

            file_path: Optional[Path] = None
            if entry:
                file_path = Path(entry["path"])
            else:
                sample_files = list(self.tangram_data_dir.glob(f"*{sample}*tangram*.h5ad"))
                if not sample_files:
                    sample_files = list(self.tangram_data_dir.glob(f"*{sample}*Tangram*.h5ad"))
                if sample_files:
                    file_path = Path(sample_files[0])
                    entry = {
                        "key": sample,
                        "primary_sample": sample,
                        "display_name": sample,
                        "aliases": [sample],
                        "data_type": "registered",
                        "path": file_path,
                        "file_name": file_path.name,
                        "file_stem": file_path.stem,
                    }

            if file_path is None or not file_path.exists():
                st.error(f"No registered data found for sample {sample}")
                return None

            adata = anndata.read_h5ad(str(file_path))

            adata.obs["total_counts"] = np.asarray(np.sum(adata.X, axis=1)).ravel()
            adata.obs["n_genes_by_counts"] = np.asarray(np.sum(adata.X > 0, axis=1)).ravel()

            store_key = entry.get("key", sample) if entry else sample
            self.registered_data[store_key] = adata

            if entry:
                self._register_catalog_aliases(store_key, entry, extra_alias=sample)
                if self._registered_catalog is not None and store_key not in self._registered_catalog:
                    self._registered_catalog[store_key] = entry

            if finalize:
                finalized = self._finalize_spatial_dataset(adata.copy())
                self.finalized_data[store_key] = finalized
                return finalized

            return adata
        except Exception as exc:
            st.error(f"Error loading registered data for {sample}: {exc}")
            return None

    def get_available_registered_samples(self) -> List[str]:
        """List available dataset keys for UI selectors."""
        catalog = self.get_registered_catalog()
        if catalog:
            sorted_entries = sorted(
                catalog.values(),
                key=lambda entry: entry.get("display_name", entry.get("key", "")).lower()
            )
            return [entry["key"] for entry in sorted_entries]

        return ['E14']

    def _get_spatial_sample_tokens(self, entry: Dict, sample: str) -> List[str]:
        """Collect token variants that help match auxiliary files for a sample."""
        tokens: set = set()
        candidates = [
            sample,
            entry.get('key'),
            entry.get('primary_sample'),
            entry.get('file_stem'),
        ]

        for candidate in candidates:
            if not candidate:
                continue
            upper_value = str(candidate).upper()
            tokens.add(upper_value)
            tokens.add(upper_value.replace('.', ''))
            tokens.update(part for part in re.split(r'[_\-]', upper_value) if part)

        primary_sample = entry.get('primary_sample', '')
        if primary_sample:
            base_match = re.match(r'(E\d+)', primary_sample, re.IGNORECASE)
            decimal_match = re.match(r'(E\d+\.\d+)', primary_sample, re.IGNORECASE)
            if base_match:
                tokens.add(base_match.group(1).upper())
            if decimal_match:
                tokens.add(decimal_match.group(1).upper())

        return [token for token in tokens if token]

    def _candidate_dapi_dirs(self, entry: Dict) -> List[Path]:
        """Determine directories likely to contain DAPI images for a sample."""
        directories: List[Path] = []
        path = Path(entry['path'])

        base_dir = self.spatial_data_dir.resolve()
        current = path.parent.resolve()
        while True:
            if current.exists() and current not in directories:
                directories.append(current)
            if current == base_dir:
                break
            if base_dir not in current.parents:
                break
            current = current.parent

        spatial_dir = base_dir / 'spatial'
        if spatial_dir.exists() and spatial_dir not in directories:
            directories.append(spatial_dir)

        if base_dir not in directories:
            directories.append(base_dir)

        return directories

    def get_dapi_image_path(self, sample: str) -> Optional[Path]:
        """Return the DAPI image path for a spatial sample if available.

        Masked DAPI images (``*_masked.tif``) are preferred when both masked and
        unmasked variants are present within a directory.
        """
        if not sample:
            return None

        sample_key = sample.lower()
        if sample_key in self._dapi_cache:
            return self._dapi_cache[sample_key]

        entry = self._resolve_registered_entry(sample)
        if entry is None:
            self.refresh_registered_catalog()
            entry = self._resolve_registered_entry(sample)

        if entry is None:
            return None

        catalog_key = (entry.get('key') or sample).lower()
        if catalog_key in self._dapi_cache:
            return self._dapi_cache[catalog_key]

        tokens = self._get_spatial_sample_tokens(entry, sample)
        best_file: Optional[Path] = None
        best_score = -1.0

        for directory in self._candidate_dapi_dirs(entry):
            dapi_files = sorted(directory.glob("DAPI*.tif*"))
            for file_path in dapi_files:
                if file_path.suffix.lower() not in ('.tif', '.tiff'):
                    continue
                name_upper = file_path.name.upper()
                if 'DAPI' not in name_upper:
                    continue

                score = 0.0
                if '_MASKED' in name_upper:
                    score += 50.0  # Prefer masked imagery when available
                for token in tokens:
                    if token and token in name_upper:
                        score += len(token)

                if '(' in file_path.name:
                    score -= 0.5

                if score > best_score:
                    best_score = score
                    best_file = file_path

            if best_score > 0:
                break

        # Cache the result for all aliases even if not found (None)
        cached_value = best_file
        for alias in set(entry.get('aliases', []) + [entry.get('key'), sample]):
            if alias:
                self._dapi_cache[alias.lower()] = cached_value

        entry['dapi_path'] = best_file
        return best_file

    # --- Dashboard Summary & Maintenance ---

    def get_data_summary(self) -> Dict:
        """Get summary of all loaded data"""
        summary = {
            'scrna_samples': list(self.scrna_data.keys()),
            'finalized_samples': list(self.finalized_data.keys()),
            'registered_samples': list(self.registered_data.keys()),
            'total_samples': len(set(list(self.scrna_data.keys()) +
                                   list(self.finalized_data.keys()) +
                                   list(self.registered_data.keys())))
        }
        return summary

    def clear_data(self, data_type: str = None, sample: str = None):
        """Clear loaded data"""
        if data_type == 'scrna' or data_type is None:
            if sample:
                self.scrna_data.pop(sample, None)
            else:
                self.scrna_data.clear()
        
        if data_type == 'spatial' or data_type is None:
            if sample:
                spatial_key = self._resolve_registered_key(sample, finalized=True)
                if spatial_key:
                    self.finalized_data.pop(spatial_key, None)
                    aliases_to_purge = [alias for alias, key in list(self._registered_alias_map.items()) if key == spatial_key]
                    for alias in aliases_to_purge:
                        self._dapi_cache.pop(alias.lower(), None)
                    self._dapi_cache.pop(spatial_key.lower(), None)
            else:
                self.finalized_data.clear()
                self._dapi_cache.clear()

        if data_type == 'tangram' or data_type is None:
            if sample:
                canonical = self._registered_alias_map.get(sample.lower(), sample)
                self.registered_data.pop(canonical, None)
                self.finalized_data.pop(canonical, None)
                aliases_to_remove = [alias for alias, key in list(self._registered_alias_map.items()) if key == canonical]
                for alias in aliases_to_remove:
                    self._registered_alias_map.pop(alias, None)
                    self._dapi_cache.pop(alias, None)
                self._dapi_cache.pop(canonical.lower(), None)
            else:
                self.registered_data.clear()
                self.finalized_data.clear()
                self._registered_alias_map.clear()
                self._dapi_cache.clear()
