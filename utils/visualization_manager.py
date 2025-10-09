#!/usr/bin/env python3
"""
Visualization Manager Utility
Handles all plotting and visualization functions for the dashboard
"""

import math

from typing import Optional

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import anndata
import matplotlib.pyplot as plt
import scanpy as sc

class VisualizationManager:
    """Manages visualizations for all data types"""
    
    @staticmethod
    def plot_umap(adata, groupby="active.ident", title: Optional[str] = None):
        """Create UMAP plot"""
        try:
            # Check if UMAP coordinates exist
            if 'X_umap' in adata.obsm:
                umap_coords = adata.obsm['X_umap']
            elif 'umap' in adata.obsm:
                umap_coords = adata.obsm['umap']
            else:
                st.warning("No UMAP coordinates found. Please run UMAP first.")
                return None
            
            # Get grouping variable
            if groupby in adata.obs.columns:
                group_values = adata.obs[groupby]
            else:
                st.warning(f"Grouping variable '{groupby}' not found. Using default.")
                group_values = adata.obs.get('active.ident', pd.Series(['Unknown'] * len(adata)))
            
            # Create plot
            scatter_kwargs = dict(
                x=umap_coords[:, 0],
                y=umap_coords[:, 1],
                color=group_values,
                labels={'x': 'UMAP1', 'y': 'UMAP2'}
            )
            if title:
                scatter_kwargs['title'] = title

            fig = px.scatter(**scatter_kwargs)
            fig.update_layout(height=600)
            return fig
        except Exception as e:
            st.error(f"Error creating UMAP plot: {e}")
            return None

    @staticmethod
    def plot_violinplot(adata, groupby="active.ident", feature="nCount_SCT",title="Violin Plot"):
        """Create violin plot"""
        try:
            if groupby in adata.obs.columns and feature in adata.obs.columns:
                keys = [feature]
                fig = sc.pl.violin(
                    adata,
                    keys=keys,
                    groupby=groupby,
                    rotation=90,
                    show=False
                )

                # Get the current figure and return it
                fig = plt.gcf()  # Get current figure
                
                st.pyplot(fig)
            else:
                st.warning("Check that 'groupby' is in adata.obs and 'objects' is a valid obs column or gene.")

        except Exception as e:
            st.error(f"Error creating violin plot: {e}")
            return None

    # NOTE: Unused across dashboard pages; keeping implementation for potential future use.
    # @staticmethod
    # def plot_violinplot_plotly(adata, groupby="active.ident", feature="nCount_SCT", title="Violin Plot"):
    #     """Create violin plot using Plotly (interactive)"""
    #     try:
    #         if groupby not in adata.obs.columns or feature not in adata.obs.columns:
    #             st.warning("Check that 'groupby' and 'feature' exist in adata.obs.")
    #             return None
    #         df = pd.DataFrame({
    #             groupby: adata.obs[groupby].astype(str).values,
    #             feature: adata.obs[feature].values
    #         })
    #         fig = px.violin(
    #             df,
    #             x=groupby,
    #             y=feature,
    #             color=groupby,
    #             title=title
    #         )

    #         fig.update_layout(height=600, showlegend=False)
    #         return fig
    #     except Exception as e:
    #         st.error(f"Error creating Plotly violin plot: {e}")
    #         return None

    @staticmethod
    def plot_dotplot(adata, groupby="active.ident", genes=None):
        """Create dot plot"""
        try:
            if genes is None:
                # Get top genes by default
                top_genes = adata.var['total_counts'].nlargest(10).index.tolist()
                genes = top_genes[:5]  # Use top 5 genes
            
            # Get grouping variable
            if groupby in adata.obs.columns:
                group_values = adata.obs[groupby]
            else:
                st.warning(f"Grouping variable '{groupby}' not found.")
                return None
            
            # Prepare data for dot plot
            plot_data = []
            for gene in genes:
                if gene in adata.var_names:
                    gene_idx = list(adata.var_names).index(gene)
                    expression = adata.X[:, gene_idx]
                    
                    for group in group_values.unique():
                        group_mask = group_values == group
                        group_expression = expression[group_mask]
                        
                        # Calculate mean expression and percentage
                        mean_expr = np.mean(group_expression)
                        pct_expr = np.sum(group_expression > 0) / len(group_expression) * 100
                        
                        plot_data.append({
                            'Gene': gene,
                            'Group': group,
                            'Mean_Expression': mean_expr,
                            'Percent_Expressed': pct_expr
                        })
            
            if plot_data:
                df = pd.DataFrame(plot_data)
                
                # Create dot plot using scatter
                fig = px.scatter(
                    df,
                    x='Group',
                    y='Gene',
                    size='Percent_Expressed',
                    color='Mean_Expression',
                    color_continuous_scale='viridis',
                    size_max=20
                )
                fig.update_layout(height=500)
                return fig
            else:
                st.warning("No valid genes found for dot plot")
                return None
        except Exception as e:
            st.error(f"Error creating dot plot: {e}")
            return None
    
    @staticmethod
    def plot_feature_plot(adata, genes=None, title="Feature Plot"):
        """Create feature plot (UMAP with gene expression)"""
        try:
            if genes is None:
                # Get top gene by default
                top_gene = adata.var['total_counts'].nlargest(1).index[0]
                genes = [top_gene]

            # Check if UMAP coordinates exist
            if 'X_umap' in adata.obsm:
                umap_coords = adata.obsm['X_umap']
            elif 'umap' in adata.obsm:
                umap_coords = adata.obsm['umap']
            else:
                st.warning("No UMAP coordinates found for feature plot.")
                return None

            # Create subplots for each gene
            n_genes = len(genes)
            fig = make_subplots(
                rows=1, cols=n_genes,
                subplot_titles=genes,
                specs=[[{"secondary_y": False}] * n_genes]
            )

            for i, gene in enumerate(genes):
                if gene in adata.var_names:
                    gene_idx = list(adata.var_names).index(gene)
                    expression = adata.X[:, gene_idx]

                    # Convert to array if sparse
                    if hasattr(expression, 'toarray'):
                        expression = expression.toarray().flatten()

                    fig.add_trace(
                        go.Scatter(
                            x=umap_coords[:, 0],
                            y=umap_coords[:, 1],
                            mode='markers',
                            marker=dict(
                                color=expression,
                                colorscale='viridis',
                                showscale=True if i == 0 else False
                            ),
                            name=gene
                        ),
                        row=1, col=i+1
                    )

            fig.update_layout(height=500, title=title)
            return fig
        except Exception as e:
            st.error(f"Error creating feature plot: {e}")
            return None

    @staticmethod
    def plot_feature_plot_scanpy(adata, genes=None):
        """Create feature plots using Scanpy's matplotlib interface."""
        try:
            if genes is None or len(genes) == 0:
                top_gene = adata.var['total_counts'].nlargest(1).index[0]
                genes = [top_gene]

            n_cols = 3
            n_rows = math.ceil(len(genes) / n_cols)

            fig, axs = plt.subplots(n_rows, n_cols, figsize=(3 * n_cols, 3 * n_rows))
            axs = np.atleast_2d(axs)

            for idx, gene in enumerate(genes):
                row = idx // n_cols
                col = idx % n_cols
                ax = axs[row, col]

                sc.pl.umap(
                    adata,
                    color=gene,
                    ax=ax,
                    show=False,
                    frameon=False,
                    s=20,
                )

                ax.axis("off")
                ax.set_title(gene)

            total_plots = n_rows * n_cols
            for extra_idx in range(len(genes), total_plots):
                row = extra_idx // n_cols
                col = extra_idx % n_cols
                fig.delaxes(axs[row, col])

            fig.tight_layout()
            return fig
        except Exception as exc:
            st.error(f"Error creating Scanpy feature plots: {exc}")
            return None

    @staticmethod
    def plot_celltype_proportion(adata, sample, title="Cell Type Proportion"):
        """Create cell type proportion plot"""
        try:
            # Look for cell type column
            celltype_cols = [col for col in adata.obs.columns if 'cell' in col.lower() or 'type' in col.lower() or 'ident' in col.lower()]
            
            if celltype_cols:
                celltype_col = celltype_cols[0]
                cell_counts = adata.obs[celltype_col].value_counts()
                
                fig = px.pie(
                    values=cell_counts.values,
                    names=cell_counts.index,
                    title=f"{title} - {sample}"
                )
                fig.update_layout(height=500)
                return fig
            else:
                st.warning("No cell type column found")
                return None
        except Exception as e:
            st.error(f"Error creating cell type proportion plot: {e}")
            return None
    
    # NOTE: Unused across dashboard pages; keeping implementation for potential future use.
    # @staticmethod
    # def plot_spatial_scatter(adata, color_by="active.ident", title="Spatial Scatter Plot"):
    #     """Create spatial scatter plot"""
    #     try:
    #         if 'spatial_x' not in adata.obs.columns or 'spatial_y' not in adata.obs.columns:
    #             st.warning("No spatial coordinates found")
    #             return None
    #         
    #         # Get color variable
    #         if color_by in adata.obs.columns:
    #             color_values = adata.obs[color_by]
    #         else:
    #             color_values = None
    #         
    #         fig = px.scatter(
    #             x=adata.obs['spatial_x'],
    #             y=adata.obs['spatial_y'],
    #             color=color_values,
    #             title=title,
    #             labels={'x': 'X Coordinate', 'y': 'Y Coordinate'}
    #         )
    #         fig.update_layout(height=600)
    #         return fig
    #     except Exception as e:
    #         st.error(f"Error creating spatial scatter plot: {e}")
    #         return None
    
    # NOTE: Unused across dashboard pages; keeping implementation for potential future use.
    # @staticmethod
    # def plot_neighborhood_enrichment(adata, groupby="active.ident", title="Neighborhood Enrichment"):
    #     """Create neighborhood enrichment plot (placeholder)"""
    #     try:
    #         # This is a placeholder - you would need squidpy for actual neighborhood analysis
    #         st.info("Neighborhood enrichment analysis requires squidpy. This is a placeholder plot.")
    #         
    #         # Create a simple correlation matrix as placeholder
    #         if groupby in adata.obs.columns:
    #             group_values = adata.obs[groupby]
    #             group_counts = group_values.value_counts()
    #             
    #             fig = px.bar(
    #                 x=group_counts.index,
    #                 y=group_counts.values,
    #                 title=f"{title} - {groupby}",
    #                 labels={'x': 'Cell Type', 'y': 'Count'}
    #             )
    #             fig.update_layout(height=500)
    #             return fig
    #         else:
    #             st.warning(f"Grouping variable '{groupby}' not found")
    #             return None
    #     except Exception as e:
    #         st.error(f"Error creating neighborhood enrichment plot: {e}")
    #         return None
    
    # NOTE: Unused across dashboard pages; keeping implementation for potential future use.
    # @staticmethod
    # def plot_quality_metrics(adata, title="Quality Metrics"):
    #     """Plot quality control metrics"""
    #     try:
    #         fig = make_subplots(
    #             rows=2, cols=2,
    #             subplot_titles=('Total Counts Distribution', 'Genes per Cell Distribution', 
    #                            'Total Counts vs Genes', 'Top Expressed Genes'),
    #             specs=[[{"secondary_y": False}, {"secondary_y": False}],
    #                    [{"secondary_y": False}, {"secondary_y": False}]]
    #         )
    #         
    #         # Total counts histogram
    #         fig.add_trace(
    #             go.Histogram(x=adata.obs['total_counts'], name='Total Counts', opacity=0.7),
    #             row=1, col=1
    #         )
    #         
    #         # Genes per cell histogram
    #         fig.add_trace(
    #             go.Histogram(x=adata.obs['n_genes_by_counts'], name='Genes per Cell', opacity=0.7),
    #             row=1, col=2
    #         )
    #         
    #         # Scatter plot
    #         fig.add_trace(
    #             go.Scatter(x=adata.obs['total_counts'], y=adata.obs['n_genes_by_counts'], 
    #                       mode='markers', name='Cells', opacity=0.6),
    #             row=2, col=1
    #         )
    #         
    #         # Top genes bar plot
    #         top_10_genes = adata.var['total_counts'].nlargest(10)
    #         fig.add_trace(
    #             go.Bar(x=top_10_genes.index, y=top_10_genes.values, name='Top Genes'),
    #             row=2, col=2
    #         )
    #         
    #         fig.update_layout(height=800, showlegend=False, title=title)
    #         return fig
    #     except Exception as e:
    #         st.error(f"Error creating quality metrics plot: {e}")
    #         return None
    
    # NOTE: Unused across dashboard pages; keeping implementation for potential future use.
    # @staticmethod
    # def plot_gene_expression(adata, gene_names, plot_type='histogram'):
    #     """Plot gene expression"""
    #     try:
    #         available_genes = [gene for gene in gene_names if gene in adata.var_names]
    #         if not available_genes:
    #             return None
    #         
    #         gene_indices = [list(adata.var_names).index(gene) for gene in available_genes]
    #         expression_data = adata.X[:, gene_indices]
    #         
    #         if plot_type == 'histogram':
    #             fig = make_subplots(rows=1, cols=len(available_genes), 
    #                                subplot_titles=available_genes)
    #             
    #             for i, gene in enumerate(available_genes):
    #                 fig.add_trace(
    #                     go.Histogram(x=expression_data[:, i], name=gene, opacity=0.7),
    #                     row=1, col=i+1
    #                 )
    #             
    #             fig.update_layout(height=500, showlegend=False)
    #             return fig
    #         
    #         elif plot_type == 'violin':
    #             fig = go.Figure()
    #             
    #             for i, gene in enumerate(available_genes):
    #                 fig.add_trace(go.Violin(
    #                     y=expression_data[:, i],
    #                     name=gene,
    #                     box_visible=True,
    #                     meanline_visible=True
    #                 ))
    #             
    #             fig.update_layout(height=500, title="Gene Expression Distribution")
    #             return fig
    #     except Exception as e:
    #         st.error(f"Error creating gene expression plot: {e}")
    #         return None
