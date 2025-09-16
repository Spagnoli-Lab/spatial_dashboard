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
        
        self.scrna_data = {}  # Dictionary to store data by sample
        self.spatial_data = {}  # Dictionary to store data by sample
        self.tangram_data = {}  # Dictionary to store data by sample
        self._spatial_catalog: Optional[Dict[str, Dict]] = None
        self._spatial_alias_map: Dict[str, str] = {}
        self._dapi_cache: Dict[str, Optional[Path]] = {}

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

    def _build_spatial_catalog(self) -> Dict[str, Dict]:
        """Build a catalog of available spatial datasets."""
        catalog: Dict[str, Dict] = {}

        candidate_paths: List[Path] = []
        candidate_paths.extend(sorted(self.spatial_data_dir.glob("*.h5ad")))

        spatial_dir = self.spatial_data_dir / "spatial"
        if spatial_dir.exists():
            candidate_paths.extend(sorted(spatial_dir.rglob("*.h5ad")))

        for file_path in candidate_paths:
            if not file_path.is_file():
                continue

            if any(part.lower() == 'scrna-seq' for part in file_path.parts):
                continue

            primary_sample, display_name, aliases, data_type = self._infer_spatial_metadata(file_path)

            key = primary_sample or file_path.stem
            base_key = key
            counter = 2
            while key in catalog:
                key = f"{base_key}_{counter}"
                counter += 1

            alias_values = set(aliases)
            alias_values.update(filter(None, [primary_sample, file_path.stem, file_path.name, key]))

            entry = {
                'key': key,
                'primary_sample': primary_sample,
                'display_name': display_name,
                'aliases': sorted(alias_values),
                'data_type': data_type,
                'path': file_path,
                'file_name': file_path.name,
                'file_stem': file_path.stem,
            }
            catalog[key] = entry

        return catalog

    def refresh_spatial_catalog(self):
        """Force rebuilding of the spatial catalog on next access."""
        self._spatial_catalog = None

    def get_spatial_catalog(self) -> Dict[str, Dict]:
        """Return the cached spatial catalog, rebuilding if needed."""
        if self._spatial_catalog is None:
            self._spatial_catalog = self._build_spatial_catalog()
        return self._spatial_catalog

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

    def _resolve_spatial_entry(self, sample: str) -> Optional[Dict]:
        """Resolve a sample identifier to a catalog entry."""
        if not sample:
            return None

        catalog = self.get_spatial_catalog()

        if sample in catalog:
            return catalog[sample]

        sample_lower = sample.lower()
        for entry in catalog.values():
            candidate_aliases = [
                entry.get('key'),
                entry.get('primary_sample'),
                entry.get('file_stem'),
                entry.get('file_name'),
            ] + entry.get('aliases', [])

            for alias in candidate_aliases:
                if alias and alias.lower() == sample_lower:
                    return entry

        return None

    def _register_spatial_aliases(self, canonical_key: str, entry: Dict, extra_alias: Optional[str] = None):
        """Update alias bookkeeping for loaded spatial datasets."""
        if not canonical_key:
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

        entry['aliases'] = sorted(alias_values)

        for alias in alias_values:
            self._spatial_alias_map[alias.lower()] = canonical_key

    def _resolve_spatial_key(self, sample: str) -> Optional[str]:
        """Resolve a sample identifier to a loaded spatial dataset key."""
        if not sample:
            return None

        if sample in self.spatial_data:
            return sample

        return self._spatial_alias_map.get(sample.lower())

    def get_spatial_dataset(self, sample: str) -> Optional[anndata.AnnData]:
        """Retrieve a loaded spatial dataset by sample or alias."""
        spatial_key = self._resolve_spatial_key(sample)
        if spatial_key:
            return self.spatial_data.get(spatial_key)
        return None

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
        """Return the DAPI image path for a spatial sample if available."""
        if not sample:
            return None

        sample_key = sample.lower()
        if sample_key in self._dapi_cache:
            return self._dapi_cache[sample_key]

        entry = self._resolve_spatial_entry(sample)
        if entry is None:
            self.refresh_spatial_catalog()
            entry = self._resolve_spatial_entry(sample)

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

    def load_spatial_data(self, sample: str) -> anndata.AnnData:
        """Load spatial data for specific sample"""
        try:
            catalog_entry = self._resolve_spatial_entry(sample)
            if catalog_entry is None:
                # Catalog may be stale if new files were added; refresh once
                self.refresh_spatial_catalog()
                catalog_entry = self._resolve_spatial_entry(sample)

            file_path: Optional[Path] = None
            if catalog_entry:
                file_path = Path(catalog_entry['path'])
            else:
                # Fallback to original glob patterns for backward compatibility
                if sample == "E14":
                    tangram_file = self.spatial_data_dir / "E14.5_2_Tangram.h5ad"
                    if tangram_file.exists():
                        file_path = tangram_file

                if file_path is None:
                    sample_files = list(self.spatial_data_dir.glob(f"*{sample}*spatial*.h5ad"))
                    if not sample_files:
                        sample_files = list(self.spatial_data_dir.glob(f"*{sample}*predicted*.h5ad"))
                    if sample_files:
                        file_path = sample_files[0]

                if file_path is None:
                    st.error(f"No spatial data found for sample {sample}")
                    return None

                file_path = Path(file_path)
                catalog_entry = {
                    'key': sample,
                    'primary_sample': sample,
                    'display_name': sample,
                    'aliases': [sample],
                    'data_type': 'unknown',
                    'path': file_path,
                    'file_name': file_path.name,
                    'file_stem': file_path.stem,
                }

            if not file_path.exists():
                st.error(f"Spatial data file not found: {file_path}")
                return None

            adata = anndata.read_h5ad(str(file_path))
            adata = self._finalize_spatial_dataset(adata)

            store_key = catalog_entry.get('key', sample)
            self.spatial_data[store_key] = adata

            # Remember aliases for subsequent lookups and updates
            self._register_spatial_aliases(store_key, catalog_entry, extra_alias=sample)

            # Update cached catalog with fallback entries so selectors stay in sync
            if self._spatial_catalog is not None and store_key not in self._spatial_catalog:
                self._spatial_catalog[store_key] = catalog_entry

            return adata
        except Exception as e:
            st.error(f"Error loading spatial data for {sample}: {e}")
            return None
    
    def load_tangram_data(self, sample: str) -> anndata.AnnData:
        """Load Tangram data for specific sample"""
        try:
            # Look for Tangram files for the specific sample
            sample_files = list(self.tangram_data_dir.glob(f"*{sample}*tangram*.h5ad"))
            if not sample_files:
                sample_files = list(self.tangram_data_dir.glob(f"*{sample}*Tangram*.h5ad"))
            
            if sample_files:
                file_path = str(sample_files[0])
                adata = anndata.read_h5ad(file_path)
                
                # Calculate quality metrics
                adata.obs['total_counts'] = np.sum(adata.X, axis=1)
                adata.obs['n_genes_by_counts'] = np.sum(adata.X > 0, axis=1)
                
                self.tangram_data[sample] = adata
                return adata
            else:
                st.error(f"No Tangram data found for sample {sample}")
                return None
        except Exception as e:
            st.error(f"Error loading Tangram data for {sample}: {e}")
            return None
    
    def get_available_samples(self) -> List[str]:
        """Get list of available samples (E12, E14, E17) for spatial data"""
        catalog = self.get_spatial_catalog()
        if catalog:
            sorted_entries = sorted(
                catalog.values(),
                key=lambda entry: entry.get('display_name', entry.get('key', '')).lower()
            )
            return [entry['key'] for entry in sorted_entries]

        # Fallback to original glob-based detection if the catalog is empty
        samples = []
        for file_path in self.spatial_data_dir.glob("*.h5ad"):
            file_name = file_path.name.lower()
            if 'e14.5' in file_name and 'e14' not in samples:
                samples.append('E14')
            elif 'e12' in file_name and 'e12' not in samples:
                samples.append('E12')
            elif 'e14' in file_name and 'e14' not in samples:
                samples.append('E14')
            elif 'e17' in file_name and 'e17' not in samples:
                samples.append('E17')

        if not samples:
            samples = ['E14']

        return sorted(samples)
    
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
    
    def get_data_summary(self) -> Dict:
        """Get summary of all loaded data"""
        summary = {
            'scrna_samples': list(self.scrna_data.keys()),
            'spatial_samples': list(self.spatial_data.keys()),
            'tangram_samples': list(self.tangram_data.keys()),
            'total_samples': len(set(list(self.scrna_data.keys()) + 
                                   list(self.spatial_data.keys()) + 
                                   list(self.tangram_data.keys())))
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
                spatial_key = self._resolve_spatial_key(sample)
                if spatial_key:
                    self.spatial_data.pop(spatial_key, None)
                    aliases_to_remove = [alias for alias, key in self._spatial_alias_map.items() if key == spatial_key]
                    for alias in aliases_to_remove:
                        self._spatial_alias_map.pop(alias, None)
            else:
                self.spatial_data.clear()
                self._spatial_alias_map.clear()

        if data_type == 'tangram' or data_type is None:
            if sample:
                self.tangram_data.pop(sample, None)
            else:
                self.tangram_data.clear()
