#!/usr/bin/env python3
"""
Spatial Visualization Manager Utility
Handles squidpy-based spatial visualizations for the dashboard
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import anndata
from pathlib import Path
import sys
import copy
from typing import Optional, Tuple

# Try to import squidpy
try:
    import squidpy as sq
    import scanpy as sc
    SQUIDPY_AVAILABLE = True
except ImportError:
    SQUIDPY_AVAILABLE = False

class SpatialVisualizationManager:
    """Manages spatial visualizations using squidpy"""
    
    def __init__(self):
        """Initialize the spatial visualization manager"""
        if not SQUIDPY_AVAILABLE:
            st.warning("⚠️ Squidpy not available. Please install with: pip install squidpy")
    
    @staticmethod
    def plot_spatial_scatter(
        adata,
        color=None,
        size=20,
        shape=None,
        title="Spatial Scatter Plot",
        figsize: Optional[Tuple[float, float]] = None,
    ):
        """
        Create spatial scatter plot using squidpy
        
        Parameters:
        -----------
        adata : AnnData
            AnnData object with spatial coordinates
        color : str or list
            Column name(s) to color by
        size : int
            Point size
        shape : str
            Point shape (None, 'circle', 'square', 'triangle')
        title : str
            Plot title
        figsize : tuple(float, float), optional
            Figure size in inches (width, height)
        """
        if not SQUIDPY_AVAILABLE:
            st.error("Squidpy is required for spatial scatter plots")
            return None
        
        try:
            # Check if spatial coordinates exist
            if 'spatial' not in adata.obsm:
                st.error("No spatial coordinates found in adata.obsm['spatial']")
                return None
            
            # Create the plot
            resolved_size = figsize if figsize else (6.0, 6.0)
            fig, ax = plt.subplots(figsize=resolved_size)
            
            # Add spatial key to uns if not present
            if 'spatial' not in adata.uns:
                adata.uns['spatial'] = {'spatial': True}
            
            # Create spatial scatter plot
            sq.pl.spatial_scatter(
                adata,
                color=color,
                size=size,
                shape=shape,
                ax=ax
            )
            
            ax.set_title(title)
            plt.tight_layout()
            
            return fig
            
        except Exception as e:
            st.error(f"Error creating spatial scatter plot: {e}")
            return None
    
    @staticmethod
    def plot_neighborhood_enrichment(
        adata,
        cluster_key: str = "celltype",
        figsize: Optional[Tuple[float, float]] = None,
        dpi: int = 110,
    ):
        """
        Create neighborhood enrichment plot using squidpy
        
        Parameters:
        -----------
        adata : AnnData
            AnnData object with spatial coordinates
        cluster_key : str
            Column name for clustering
        method : str
            Clustering method
        n_perms : int
            Number of permutations
        title : str
            Plot title
        """
        if not SQUIDPY_AVAILABLE:
            st.error("Squidpy is required for neighborhood enrichment analysis")
            return None
        
        try:
            # Check if spatial coordinates exist
            if 'spatial' not in adata.obsm:
                st.error("No spatial coordinates found in adata.obsm['spatial']")
                return None
            
            # Check if cluster key exists
            if cluster_key not in adata.obs.columns:
                st.error(f"Cluster key '{cluster_key}' not found in adata.obs")
                return None
            
            # Reuse previously computed enrichment results when available.
            required_graph_keys = {'spatial_connectivities', 'spatial_distances'}
            if not required_graph_keys.issubset(adata.obsp.keys()):
                st.error(
                    "Spatial neighbor graph not found. Run the preprocessing script "
                    "to populate spatial neighbors and enrichment matrices."
                )
                return None

            cache_bucket = adata.uns.get("cached_nhood_enrichment", {})
            cached_result = cache_bucket.get(cluster_key)

            if cached_result is not None:
                adata.uns["nhood_enrichment"] = copy.deepcopy(cached_result)
            else:
                existing = adata.uns.get("nhood_enrichment")
                if not existing or existing.get("cluster_key") != cluster_key:
                    st.error(
                        "Precomputed neighborhood enrichment not found for "
                        f"cluster key '{cluster_key}'. Run the preprocessing "
                        "script and reload this dataset."
                    )
                    return None

            # Resolve desired figure size
            if figsize is None:
                fig_width, fig_height = (9.0, 9.0)
            else:
                fig_width, fig_height = figsize

            # Let squidpy create the axes but control overall scale via figsize/dpi
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            sq.pl.nhood_enrichment(
                adata,
                cluster_key=cluster_key,
                #cmap="coolwarm",
                title="Neighborhood enrichment",
                annotate=True,
                cbar_kwargs={"shrink": 0.7},
                ax=ax
            )

            return fig
            
        except Exception as e:
            st.error(f"Error creating neighborhood enrichment plot: {e}")
            return None
    
    @staticmethod
    def plot_gene_expression_spatial(
        adata,
        genes,
        size=20,
        title="Gene Expression Spatial Plot",
        figsize: Optional[Tuple[float, float]] = None,
    ):
        """
        Create spatial scatter plot colored by gene expression
        
        Parameters:
        -----------
        adata : AnnData
            AnnData object with spatial coordinates
        genes : list
            List of gene names to plot
        size : int
            Point size
        title : str
            Plot title
        figsize : tuple(float, float), optional
            Figure size in inches (width, height)
        """
        if not SQUIDPY_AVAILABLE:
            st.error("Squidpy is required for spatial gene expression plots")
            return None
        
        try:
            # Check if spatial coordinates exist
            if 'spatial' not in adata.obsm:
                st.error("No spatial coordinates found in adata.obsm['spatial']")
                return None
            
            # Filter genes to those that exist in the dataset
            available_genes = [gene for gene in genes if gene in adata.var_names]
            if not available_genes:
                st.error("None of the specified genes found in the dataset")
                return None
            
            # Add these debug lines before spatial_scatter:
            print("Spatial coordinates shape:", adata.obsm['spatial'].shape)
            print("Spatial coordinates type:", type(adata.obsm['spatial']))
            print("Spatial coordinates sample:", adata.obsm['spatial'][:5])
            print("Size parameter:", size, "Type:", type(size))
            print("Available genes for color:", available_genes)

            # Create the plot
            resolved_size = figsize if figsize else (6.0, 6.0)
            fig, ax = plt.subplots(figsize=resolved_size)
            
            # Add spatial key to uns if not present
            # Add spatial key to uns if not present
            if 'spatial' not in adata.uns:
                adata.uns['spatial'] = {
                    'sample1': {
                        'images': {},
                        'scalefactors': {
                            'hires': 1.0,    
                            'lowres': 1.0      
                        }
                    }
                }

            # Also set library_id in obs if not present
            #if 'library_id' not in adata.obs.columns:
            #    adata.obs['library_id'] = "sample1"

            # Create spatial scatter plot for gene expression
            sq.pl.spatial_scatter(
                adata,
                color=available_genes,
                size=size,
                ax=ax
            )
            
            ax.set_title(title)
            plt.tight_layout()
            
            return fig
            
        except Exception as e:
            st.error(f"Error creating gene expression spatial plot: {e}")
            return None
    
    @staticmethod
    def compute_spatial_statistics(adata, cluster_key="active.ident"):
        """
        Compute basic spatial statistics
        
        Parameters:
        -----------
        adata : AnnData
            AnnData object with spatial coordinates
        cluster_key : str
            Column name for clustering
        """
        if not SQUIDPY_AVAILABLE:
            st.error("Squidpy is required for spatial statistics")
            return None
        
        try:
            # Check if spatial coordinates exist
            if 'spatial' not in adata.obsm:
                st.error("No spatial coordinates found in adata.obsm['spatial']")
                return None
            
            # Compute spatial neighbors
            if 'spatial_neighbors' not in adata.obsp:
                sq.gr.spatial_neighbors(adata, coord_type="generic", radius=3)
            
            # Get spatial statistics
            stats = {}
            
            # Basic spatial info
            spatial_coords = adata.obsm['spatial']
            stats['spatial_range_x'] = (spatial_coords[:, 0].min(), spatial_coords[:, 0].max())
            stats['spatial_range_y'] = (spatial_coords[:, 1].min(), spatial_coords[:, 1].max())
            stats['tissue_area'] = (stats['spatial_range_x'][1] - stats['spatial_range_x'][0]) * \
                                  (stats['spatial_range_y'][1] - stats['spatial_range_y'][0])
            stats['cell_density'] = len(adata) / stats['tissue_area']
            
            # Spatial neighbors info
            if 'spatial_neighbors' in adata.obsp:
                n_neighbors = adata.obsp['spatial_neighbors'].sum(axis=1).A1
                stats['mean_neighbors'] = np.mean(n_neighbors)
                stats['median_neighbors'] = np.median(n_neighbors)
                stats['max_neighbors'] = np.max(n_neighbors)
            
            # Cluster-specific statistics
            if cluster_key in adata.obs.columns:
                cluster_stats = {}
                for cluster in adata.obs[cluster_key].unique():
                    cluster_mask = adata.obs[cluster_key] == cluster
                    cluster_coords = spatial_coords[cluster_mask]
                    
                    if len(cluster_coords) > 0:
                        cluster_stats[cluster] = {
                            'count': len(cluster_coords),
                            'spatial_center': (np.mean(cluster_coords[:, 0]), np.mean(cluster_coords[:, 1])),
                            'spatial_spread': np.std(cluster_coords, axis=0)
                        }
                
                stats['cluster_statistics'] = cluster_stats
            
            return stats
            
        except Exception as e:
            st.error(f"Error computing spatial statistics: {e}")
            return None
    
    @staticmethod
    def plot_spatial_statistics(adata, cluster_key="active.ident", title="Spatial Statistics"):
        """
        Create visualization of spatial statistics
        
        Parameters:
        -----------
        adata : AnnData
            AnnData object with spatial coordinates
        cluster_key : str
            Column name for clustering
        title : str
            Plot title
        """
        if not SQUIDPY_AVAILABLE:
            st.error("Squidpy is required for spatial statistics plots")
            return None
        
        try:
            # Compute statistics
            stats = SpatialVisualizationManager.compute_spatial_statistics(adata, cluster_key)
            if stats is None:
                return None
            
            # Create subplots
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            fig.suptitle(title, fontsize=16)
            
            # Plot 1: Cell density heatmap
            spatial_coords = adata.obsm['spatial']
            axes[0, 0].hexbin(spatial_coords[:, 0], spatial_coords[:, 1], gridsize=20, cmap='viridis')
            axes[0, 0].set_title('Cell Density Heatmap')
            axes[0, 0].set_xlabel('X Coordinate')
            axes[0, 0].set_ylabel('Y Coordinate')
            
            # Plot 2: Cluster distribution
            if cluster_key in adata.obs.columns:
                cluster_counts = adata.obs[cluster_key].value_counts()
                axes[0, 1].bar(range(len(cluster_counts)), cluster_counts.values)
                axes[0, 1].set_title('Cell Count by Cluster')
                axes[0, 1].set_xlabel('Cluster')
                axes[0, 1].set_ylabel('Cell Count')
                axes[0, 1].set_xticks(range(len(cluster_counts)))
                axes[0, 1].set_xticklabels(cluster_counts.index, rotation=45)
            
            # Plot 3: Spatial range
            axes[1, 0].scatter(spatial_coords[:, 0], spatial_coords[:, 1], alpha=0.6, s=1)
            axes[1, 0].set_title('Spatial Distribution')
            axes[1, 0].set_xlabel('X Coordinate')
            axes[1, 0].set_ylabel('Y Coordinate')
            
            # Plot 4: Statistics summary
            axes[1, 1].axis('off')
            stats_text = f"""
            Tissue Area: {stats['tissue_area']:.2f} units²
            Cell Density: {stats['cell_density']:.2f} cells/unit²
            Mean Neighbors: {stats.get('mean_neighbors', 'N/A'):.2f}
            Total Cells: {len(adata):,}
            """
            axes[1, 1].text(0.1, 0.5, stats_text, transform=axes[1, 1].transAxes, 
                           fontsize=12, verticalalignment='center')
            
            plt.tight_layout()
            return fig
            
        except Exception as e:
            st.error(f"Error creating spatial statistics plot: {e}")
            return None
    
    @staticmethod
    def plot_multiple_spatial_scatters(adata, color_options, size=20, n_cols=2, title="Multiple Spatial Plots"):
        """
        Create multiple spatial scatter plots
        
        Parameters:
        -----------
        adata : AnnData
            AnnData object with spatial coordinates
        color_options : list
            List of column names to color by
        size : int
            Point size
        n_cols : int
            Number of columns in the subplot grid
        title : str
            Plot title
        """
        if not SQUIDPY_AVAILABLE:
            st.error("Squidpy is required for multiple spatial plots")
            return None
        
        try:
            # Check if spatial coordinates exist
            if 'spatial' not in adata.obsm:
                st.error("No spatial coordinates found in adata.obsm['spatial']")
                return None
            
            # Filter to available columns
            available_colors = [col for col in color_options if col in adata.obs.columns]
            if not available_colors:
                st.error("None of the specified color columns found in the dataset")
                return None
            
            # Calculate subplot layout
            n_plots = len(available_colors)
            n_rows = (n_plots + n_cols - 1) // n_cols
            
            # Create subplots
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
            if n_plots == 1:
                axes = [axes]
            elif n_rows == 1:
                axes = axes.reshape(1, -1)
            else:
                axes = axes.flatten()
            
            # Add spatial key to uns if not present
            if 'spatial' not in adata.uns:
                adata.uns['spatial'] = {'spatial': True}
            
            # Create plots
            for i, color_col in enumerate(available_colors):
                if i < len(axes):
                    sq.pl.spatial_scatter(
                        adata,
                        color=color_col,
                        size=size,
                        ax=axes[i],
                        show=False
                    )
                    axes[i].set_title(f'Color by: {color_col}')
            
            # Hide empty subplots
            for i in range(n_plots, len(axes)):
                axes[i].set_visible(False)
            
            fig.suptitle(title, fontsize=16)
            plt.tight_layout()
            
            return fig
            
        except Exception as e:
            st.error(f"Error creating multiple spatial plots: {e}")
            return None
