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
from pathlib import Path
from typing import Dict, List, Optional, Tuple

class DataManager:
    """Manages data loading and processing for all three data types"""
    
    def __init__(self, data_dir: str = "/Users/mayongzhi/Desktop/FS_lab/dashboard/test_data/scRNA-seq"):
        self.data_dir = Path(data_dir)
        self.scrna_data = {}  # Dictionary to store data by sample
        self.spatial_data = {}  # Dictionary to store data by sample
        self.tangram_data = {}  # Dictionary to store data by sample
        
    def load_scrna_data(self, sample: str) -> anndata.AnnData:
        """Load scRNA-seq data for specific sample"""
        try:
            # Specific path for Cartana.h5ad file
            cartana_path = "/Users/mayongzhi/Desktop/FS_lab/dashboard/test_data/scRNA-seq/Cartana.h5ad"
            
            # Check if Cartana.h5ad exists and load it
            if Path(cartana_path).exists():
                file_path = cartana_path
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
                sample_files = list(self.data_dir.glob(f"*{sample}*scrna*.h5ad"))
                if not sample_files:
                    sample_files = list(self.data_dir.glob(f"*{sample}*.h5ad"))
                
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
            # Look for spatial files for the specific sample
            sample_files = list(self.data_dir.glob(f"*{sample}*spatial*.h5ad"))
            if not sample_files:
                sample_files = list(self.data_dir.glob(f"*{sample}*predicted*.h5ad"))
            
            if sample_files:
                file_path = str(sample_files[0])
                adata = anndata.read_h5ad(file_path)
                
                # Extract spatial coordinates
                if 'spatial' in adata.obsm:
                    adata.obs['spatial_x'] = adata.obsm['spatial'][:, 0]
                    adata.obs['spatial_y'] = adata.obsm['spatial'][:, 1]
                elif 'X_spatial' in adata.obsm:
                    adata.obs['spatial_x'] = adata.obsm['X_spatial'][:, 0]
                    adata.obs['spatial_y'] = adata.obsm['X_spatial'][:, 1]
                
                # Calculate quality metrics
                adata.obs['total_counts'] = np.sum(adata.X, axis=1)
                adata.obs['n_genes_by_counts'] = np.sum(adata.X > 0, axis=1)
                
                self.spatial_data[sample] = adata
                return adata
            else:
                st.error(f"No spatial data found for sample {sample}")
                return None
        except Exception as e:
            st.error(f"Error loading spatial data for {sample}: {e}")
            return None
    
    def load_tangram_data(self, sample: str) -> anndata.AnnData:
        """Load Tangram data for specific sample"""
        try:
            # Look for Tangram files for the specific sample
            sample_files = list(self.data_dir.glob(f"*{sample}*tangram*.h5ad"))
            if not sample_files:
                sample_files = list(self.data_dir.glob(f"*{sample}*Tangram*.h5ad"))
            
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
        """Get list of available samples (E12, E14, E17)"""
        samples = []
        for file_path in self.data_dir.glob("*.h5ad"):
            file_name = file_path.name.lower()
            if 'e12' in file_name and 'e12' not in samples:
                samples.append('E12')
            elif 'e14' in file_name and 'e14' not in samples:
                samples.append('E14')
            elif 'e17' in file_name and 'e17' not in samples:
                samples.append('E17')
        return sorted(samples)
    
    def get_scrna_sample_options(self) -> List[str]:
        """Get unique orig.ident values from Cartana.h5ad file for scRNA-seq sample selection"""
        try:
            cartana_path = "/Users/mayongzhi/Desktop/FS_lab/dashboard/test_data/scRNA-seq/Cartana.h5ad"
            
            if Path(cartana_path).exists():
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
                self.spatial_data.pop(sample, None)
            else:
                self.spatial_data.clear()
        
        if data_type == 'tangram' or data_type is None:
            if sample:
                self.tangram_data.pop(sample, None)
            else:
                self.tangram_data.clear()
