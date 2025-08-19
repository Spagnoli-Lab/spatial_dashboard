#!/usr/bin/env python3
"""
Spatial Data Analysis Page
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import anndata
from pathlib import Path
import sys

# Add the project root to the Python path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

# Import utilities
from utils.data_manager import DataManager
from utils.visualization_manager import VisualizationManager

st.set_page_config(
    page_title="Spatial Data - Spatial Transcriptomics Dashboard",
    page_icon="🗺️",
    layout="wide"
)

# Initialize managers
@st.cache_resource
def get_managers():
    data_manager = DataManager()
    viz_manager = VisualizationManager()
    return data_manager, viz_manager

data_manager, viz_manager = get_managers()

# Sidebar
st.sidebar.title("🗺️ Spatial Data Analysis")

# Sample selection
st.sidebar.header("📁 Sample Selection")
available_samples = data_manager.get_available_samples()

if not available_samples:
    st.sidebar.error("No sample data found in data directory")
    st.info("Please ensure your data files contain sample identifiers (E12, E14, E17)")
    st.stop()

selected_sample = st.sidebar.selectbox(
    "Select Sample:",
    available_samples,
    help="Choose the sample to analyze (E12, E14, E17)"
)

# Data status and loading
st.sidebar.header("📊 Data Status")

# Show current data status
if selected_sample in data_manager.spatial_data:
    st.sidebar.success(f"✅ Spatial\n{data_manager.spatial_data[selected_sample].n_obs} cells")
else:
    st.sidebar.info("⏳ Spatial\nNot loaded")

# Auto-load data if not loaded
if selected_sample not in data_manager.spatial_data:
    with st.spinner(f"Loading spatial data for {selected_sample}..."):
        data_manager.load_spatial_data(selected_sample)

# Manual reload button
if st.sidebar.button("🔄 Reload Data"):
    with st.spinner(f"Reloading spatial data for {selected_sample}..."):
        data_manager.load_spatial_data(selected_sample)
        if selected_sample in data_manager.spatial_data:
            st.sidebar.success(f"Reloaded {data_manager.spatial_data[selected_sample].n_obs} cells")

# Check if data is available
if selected_sample not in data_manager.spatial_data:
    st.error(f"Failed to load spatial data for {selected_sample}")
    st.info("Please check that spatial data files exist for this sample")
    st.stop()

# Main content
adata = data_manager.spatial_data[selected_sample]

st.title(f"🗺️ Spatial Data - {selected_sample}")

# Summary statistics
col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric("Cells", f"{adata.n_obs:,}")
with col2:
    st.metric("Genes", f"{adata.n_vars:,}")
with col3:
    st.metric("Avg Counts/Cell", f"{np.mean(adata.obs['total_counts']):.1f}")
with col4:
    if 'spatial_x' in adata.obs.columns:
        st.metric("Spatial Coords", "Available")
    else:
        st.metric("Spatial Coords", "Not found")

# Spatial Scatter Plot
st.subheader("🗺️ Spatial Scatter Plot")

col1, col2 = st.columns([3, 1])

with col1:
    # Get available color variables
    color_options = ['active.ident'] + [col for col in adata.obs.columns if col not in ['spatial_x', 'spatial_y', 'active.ident']]
    selected_color = st.selectbox("Color by:", color_options, key="spatial_color")
    
    spatial_fig = viz_manager.plot_spatial_scatter(adata, color_by=selected_color, title=f"Spatial Scatter - {selected_sample}")
    if spatial_fig:
        st.plotly_chart(spatial_fig, use_container_width=True)

# Neighborhood Enrichment Plot
st.subheader("🔗 Neighborhood Enrichment")

col1, col2 = st.columns([3, 1])

with col1:
    # Get available grouping variables
    groupby_options = ['active.ident'] + [col for col in adata.obs.columns if col != 'active.ident']
    selected_groupby = st.selectbox("Group by:", groupby_options, key="neighborhood_groupby")
    
    neighborhood_fig = viz_manager.plot_neighborhood_enrichment(adata, groupby=selected_groupby, title=f"Neighborhood Enrichment - {selected_sample}")
    if neighborhood_fig:
        st.plotly_chart(neighborhood_fig, use_container_width=True)

# Cell Type Proportion Plot
st.subheader("🥧 Cell Type Proportion")

col1, col2 = st.columns([3, 1])

with col1:
    sample_options = ['E12', 'E14', 'E17']
    selected_sample_day = st.selectbox("Select sample day:", sample_options, index=sample_options.index(selected_sample), key="spatial_proportion")
    
    proportion_fig = viz_manager.plot_celltype_proportion(adata, selected_sample_day, title="Cell Type Proportion")
    if proportion_fig:
        st.plotly_chart(proportion_fig, use_container_width=True)

# UMAP Plot
st.subheader("📊 UMAP Plot")

col1, col2 = st.columns([3, 1])

with col1:
    # Get available grouping variables
    groupby_options = ['active.ident'] + [col for col in adata.obs.columns if col != 'active.ident']
    selected_groupby = st.selectbox("Group by:", groupby_options, key="spatial_umap_groupby")
    
    umap_fig = viz_manager.plot_umap(adata, groupby=selected_groupby, title=f"UMAP - {selected_sample}")
    if umap_fig:
        st.plotly_chart(umap_fig, use_container_width=True)

# Data information sidebar
with col2:
    st.subheader("📋 Spatial Data Info")
    st.write(f"**Sample**: {selected_sample}")
    st.write(f"**Cells**: {adata.n_obs:,}")
    st.write(f"**Genes**: {adata.n_vars:,}")
    
    # Spatial coordinates info
    if 'spatial_x' in adata.obs.columns and 'spatial_y' in adata.obs.columns:
        st.write("**Spatial Coordinates**: ✅ Available")
        st.write(f"**X Range**: {adata.obs['spatial_x'].min():.1f} - {adata.obs['spatial_x'].max():.1f}")
        st.write(f"**Y Range**: {adata.obs['spatial_y'].min():.1f} - {adata.obs['spatial_y'].max():.1f}")
    else:
        st.write("**Spatial Coordinates**: ❌ Not found")
    
    # Cell type information
    if 'active.ident' in adata.obs.columns:
        cell_types = adata.obs['active.ident'].value_counts()
        st.write("**Cell Types**:")
        for ct, count in cell_types.head(5).items():
            st.write(f"  - {ct}: {count}")
        if len(cell_types) > 5:
            st.write(f"  - ... and {len(cell_types) - 5} more")
    
    # Quality metrics
    st.subheader("📊 Quality Metrics")
    st.write(f"**Mean Counts/Cell**: {np.mean(adata.obs['total_counts']):.1f}")
    st.write(f"**Median Counts/Cell**: {np.median(adata.obs['total_counts']):.1f}")
    st.write(f"**Mean Genes/Cell**: {np.mean(adata.obs['n_genes_by_counts']):.1f}")
    st.write(f"**Median Genes/Cell**: {np.median(adata.obs['n_genes_by_counts']):.1f}")
    
    # Spatial statistics
    if 'spatial_x' in adata.obs.columns:
        st.subheader("🗺️ Spatial Statistics")
        st.write(f"**Tissue Area**: {np.ptp(adata.obs['spatial_x']) * np.ptp(adata.obs['spatial_y']):.1f} units²")
        st.write(f"**Cell Density**: {len(adata) / (np.ptp(adata.obs['spatial_x']) * np.ptp(adata.obs['spatial_y'])):.2f} cells/unit²")
    
    # Top genes
    st.subheader("🔝 Top Expressed Genes")
    top_10_genes = adata.var['total_counts'].nlargest(10)
    for i, (gene, counts) in enumerate(top_10_genes.items(), 1):
        st.write(f"{i}. {gene}: {counts:,.0f}")

# Additional spatial analysis options
st.subheader("🔬 Additional Spatial Analysis")

col1, col2 = st.columns(2)

with col1:
    st.markdown("""
    ### Available Analyses
    
    - **Spatial Clustering**: Identify spatial domains
    - **Spatial Autocorrelation**: Measure spatial correlation
    - **Cell-Cell Interactions**: Analyze neighbor relationships
    - **Spatial Gene Expression**: Map gene expression in space
    
    ### Export Options
    
    - **Spatial Coordinates**: Export cell positions
    - **Expression Matrix**: Export gene expression data
    - **Cell Metadata**: Export cell annotations
    - **Plots**: Save visualizations as images
    """)

with col2:
    st.markdown("""
    ### Analysis Parameters
    
    - **Neighborhood Radius**: Adjust spatial neighborhood size
    - **Clustering Method**: Choose clustering algorithm
    - **Gene Selection**: Select genes for spatial analysis
    - **Cell Type Filter**: Filter by specific cell types
    
    ### Visualization Settings
    
    - **Color Schemes**: Customize plot colors
    - **Point Sizes**: Adjust cell point sizes
    - **Transparency**: Set plot transparency levels
    - **Legend Position**: Customize legend placement
    """)
