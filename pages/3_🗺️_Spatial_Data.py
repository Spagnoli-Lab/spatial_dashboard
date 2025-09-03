#!/usr/bin/env python3
"""
Spatial Data Analysis Page
"""

import streamlit as st
from streamlit import session_state as ss
import pandas as pd
import numpy as np
import anndata
from pathlib import Path
import sys
import matplotlib.pyplot as plt

# Add the project root to the Python path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

# Import utilities
from utils.data_manager import DataManager
from utils.visualization_manager import VisualizationManager
from utils.spatial_visualization_manager import SpatialVisualizationManager

# Try to import squidpy
try:
    import squidpy as sq
    SQUIDPY_AVAILABLE = True
    # Show Squidpy version for debugging
    import squidpy
    st.sidebar.info(f"Squidpy version: {squidpy.__version__}")
    st.sidebar.success("✅ Spatial Visualization Manager available")
except ImportError:
    SQUIDPY_AVAILABLE = False
    st.sidebar.warning("⚠️ Squidpy not available. Spatial analysis features will be limited.")

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
    spatial_viz_manager = SpatialVisualizationManager()
    return data_manager, viz_manager, spatial_viz_manager

data_manager, viz_manager, spatial_viz_manager = get_managers()

# Sidebar
st.sidebar.title("🗺️ Spatial Data Analysis")

# Sample selection
st.sidebar.header("📁 Sample Selection")

# Debug: Show what files are found
# Safety check for new attributes
if hasattr(data_manager, 'spatial_data_dir'):
    st.sidebar.info(f"Spatial data directory: {data_manager.spatial_data_dir}")
    h5ad_files = list(data_manager.spatial_data_dir.glob("*.h5ad"))
    st.sidebar.info(f"Found H5AD files: {[f.name for f in h5ad_files]}")
else:
    # Fallback for old cached version
    st.sidebar.info(f"Data directory: {data_manager.data_dir}")
    h5ad_files = list(data_manager.data_dir.glob("*.h5ad"))
    st.sidebar.info(f"Found H5AD files: {[f.name for f in h5ad_files]}")
    st.sidebar.warning("⚠️ Using cached DataManager. Please restart the app for full functionality.")

# Get available samples with safety check
if hasattr(data_manager, 'spatial_data_dir'):
    available_samples = data_manager.get_available_samples()
else:
    # Fallback: manually check for E14.5_2_Tangram.h5ad
    available_samples = ['E14']  # Default since we know this file exists
    st.sidebar.warning("⚠️ Using fallback sample detection")

st.sidebar.info(f"Available samples: {available_samples}")

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
        try:
            data_manager.load_spatial_data(selected_sample)
        except Exception as e:
            st.error(f"Error loading data: {e}")
            st.info("Please restart the app to use the updated DataManager")
            st.stop()

# Manual reload button
if st.sidebar.button("🔄 Reload Data"):
    with st.spinner(f"Reloading spatial data for {selected_sample}..."):
        data_manager.load_spatial_data(selected_sample)
        if selected_sample in data_manager.spatial_data:
            st.sidebar.success(f"Reloaded {data_manager.spatial_data[selected_sample].n_obs} cells")

# Available Spatial Visualizations
if SQUIDPY_AVAILABLE:
    st.sidebar.header("🔬 Available Spatial Visualizations")
    st.sidebar.markdown("""
    - **Spatial Scatter Plots**: Color by annotations
    - **Neighborhood Enrichment**: Spatial relationships
    - **Gene Expression**: Spatial gene mapping
    - **Multiple Plots**: Side-by-side comparison
    - **Spatial Statistics**: Comprehensive analysis
    """)

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


# Squidpy Spatial Scatter Plot
if SQUIDPY_AVAILABLE:
    st.subheader("🔬 Squidpy Spatial Scatter Plot")
    
    # Plot settings above the figure
    st.markdown("### Plot Settings")
    
    col1, col2 = st.columns(2)
    
    with col1:
        # Color options for squidpy
        squidpy_color_options = ['celltype'] + [col for col in adata.obs.columns if col not in ['spatial_x', 'spatial_y', 'celltype']]
        squidpy_color = st.selectbox("Color by:", squidpy_color_options, key="squidpy_color")
    
    with col2:
        # Point size
        point_size = st.slider("Point size:", min_value=1, max_value=50, value=20, key="squidpy_size")
    
    # Generate plot automatically when data is loaded
    with st.spinner("Generating Squidpy spatial scatter plot..."):
        fig = spatial_viz_manager.plot_spatial_scatter(
            adata,
            color=squidpy_color,
            size=point_size,
            title=f"Squidpy Spatial Scatter - {selected_sample}"
        )
        if fig:
            st.pyplot(fig)
            plt.close(fig)
        else:
            st.warning("Could not generate spatial scatter plot. Check if spatial coordinates are available.")
else:
    st.info("💡 Install squidpy to enable advanced spatial analysis features: `pip install squidpy`")

# Neighborhood Enrichment Analysis
st.subheader("🔬 Neighborhood Enrichment Analysis")

# Squidpy Neighborhood Enrichment Analysis
if SQUIDPY_AVAILABLE:
    
    # Analysis settings above the figure
    st.markdown("### Analysis Settings")
    
    # Cluster key options
    cluster_key_options = ['celltype'] + [col for col in adata.obs.columns if col != 'celltype']
    selected_cluster_key = st.selectbox("Cluster key:", cluster_key_options, key="squidpy_cluster_key")
    
    # Generate analysis automatically when data is loaded
    with st.spinner("Computing neighborhood enrichment..."):
        fig = spatial_viz_manager.plot_neighborhood_enrichment(
            adata,
            cluster_key=selected_cluster_key
        )
        if fig:
            st.pyplot(fig)
            plt.close(fig)
            st.success(f"✅ Neighborhood enrichment computed")
        else:
            st.warning("Could not compute neighborhood enrichment. This may require specific data structure.")
else:
    st.info("💡 Install squidpy to enable neighborhood enrichment analysis: `pip install squidpy`")

# Tissuumap Webpage Embedding
st.subheader("🗺️ Tissuumap Interactive Visualization")

# Check if tissuumap files exist
tissuumap_path = Path(__file__).parent.parent / "tissuumaps"
tissuumap_index = tissuumap_path / "index.html"

if tissuumap_index.exists():
    st.success("✅ Tissuumap files found")
    
    # Display tissuumap information
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.markdown("""
        ### Interactive Spatial Visualization
        
        The embedded Tissuumap below provides an interactive visualization of your spatial data.
        You can:
        - **Zoom and Pan**: Navigate through the spatial tissue
        - **Select Cells**: Click on individual cells for details
        - **Toggle Layers**: Show/hide different data layers
        - **Export Views**: Save current visualization state
        """)
    
    with col2:
        st.markdown("### Tissuumap Features")
        st.markdown("""
        - 🎯 **Interactive Navigation**
        - 🔍 **Multi-scale Zooming**
        - 📊 **Data Overlays**
        - 💾 **Export Capabilities**
        - 🎨 **Custom Styling**
        """)
    
    # Embed the tissuumap webpage
    st.markdown("---")
    st.markdown("### 🗺️ Interactive Tissuumap")
    
    # Embed the tissuumap HTML file
    with open(tissuumap_index, 'r', encoding='utf-8') as f:
        html_content = f.read()
    
    # Display the tissuumap
    st.components.v1.html(html_content, height=600, scrolling=True)

# Gene Expression Spatial Plot
if SQUIDPY_AVAILABLE:
    st.subheader("🧬 Gene Expression Spatial Plot")
    
    # Gene selection settings above the figure
    st.markdown("### Gene Selection")
    
    col1, col2 = st.columns(2)
    
    with col1:
        # Get available genes
        available_genes = list(adata.var_names[:20])  # First 20 genes for performance
        selected_genes = st.multiselect(
            "Select genes:", 
            available_genes, 
            default=available_genes[:2] if len(available_genes) >= 2 else available_genes,
            key="spatial_genes"
        )
    
    with col2:
        # Point size
        gene_point_size = st.slider("Point size:", min_value=1, max_value=50, value=20, key="gene_point_size")
    
    # Generate plot automatically when genes are selected
    if selected_genes:
        with st.spinner("Generating gene expression spatial plot..."):
            fig = spatial_viz_manager.plot_gene_expression_spatial(
                adata,
                genes=selected_genes,
                size=gene_point_size,
                title=f"Gene Expression in Space - {selected_sample}"
            )
            if fig:
                st.pyplot(fig)
                plt.close(fig)
            else:
                st.warning("Could not generate gene expression plot. Check if spatial coordinates are available.")
    else:
        st.info("Please select at least one gene to generate the plot")
else:
    st.info("💡 Install squidpy to enable gene expression spatial analysis: `pip install squidpy`")

# Spatial Statistics Dashboard
if SQUIDPY_AVAILABLE:
    st.subheader("📊 Spatial Statistics Dashboard")
    
    # Statistics settings above the figure
    st.markdown("### Statistics Settings")
    
    # Cluster key for statistics
    stats_cluster_key = st.selectbox("Cluster key for statistics:", cluster_key_options, key="stats_cluster_key")
    
    # Generate statistics automatically when data is loaded
    with st.spinner("Computing spatial statistics..."):
        # Compute statistics
        stats = spatial_viz_manager.compute_spatial_statistics(adata, cluster_key=stats_cluster_key)
        if stats:
            st.json(stats)
        
        # Create statistics plot
        fig = spatial_viz_manager.plot_spatial_statistics(
            adata,
            cluster_key=stats_cluster_key,
            title=f"Spatial Statistics - {selected_sample}"
        )
        if fig:
            st.pyplot(fig)
            plt.close(fig)
        else:
            st.warning("Could not generate spatial statistics. This may require specific data structure.")

# Multiple Spatial Plots Comparison
if SQUIDPY_AVAILABLE:
    st.subheader("🔄 Multiple Spatial Plots Comparison")
    
    # Comparison settings above the figure
    st.markdown("### Comparison Settings")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        # Define color options for comparison (same as squidpy options)
        color_options = ['celltype'] + [col for col in adata.obs.columns if col not in ['spatial_x', 'spatial_y', 'celltype']]
        
        # Select multiple color options for comparison
        multiple_color_options = st.multiselect(
            "Select annotations to compare:", 
            color_options, 
            default=color_options[:3] if len(color_options) >= 3 else color_options,
            key="multiple_colors"
        )
    
    with col2:
        # Point size for multiple plots
        multi_point_size = st.slider("Point size:", min_value=1, max_value=50, value=15, key="multi_point_size")
    
    with col3:
        # Number of columns
        n_cols = st.selectbox("Number of columns:", [1, 2, 3], index=1, key="n_cols")
    
    # Generate comparison plots automatically when annotations are selected
    if multiple_color_options:
        with st.spinner("Generating multiple spatial plots..."):
            fig = spatial_viz_manager.plot_multiple_spatial_scatters(
                adata,
                color_options=multiple_color_options,
                size=multi_point_size,
                n_cols=n_cols,
                title=f"Multiple Spatial Visualizations - {selected_sample}"
            )
            if fig:
                st.pyplot(fig)
                plt.close(fig)
            else:
                st.warning("Could not generate comparison plots. Check if spatial coordinates are available.")
    else:
        st.info("Please select at least one annotation for comparison")

# Cell Type Proportion Plot
st.subheader("🥧 Cell Type Proportion")

# Sample selection above the figure
sample_options = ['E12', 'E14', 'E17']
selected_sample_day = st.selectbox("Select sample day:", sample_options, index=sample_options.index(selected_sample), key="spatial_proportion")

proportion_fig = viz_manager.plot_celltype_proportion(adata, selected_sample_day, title="Cell Type Proportion")
if proportion_fig:
    st.plotly_chart(proportion_fig, use_container_width=True)

# UMAP Plot
st.subheader("📊 UMAP Plot")

# Grouping selection above the figure
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
    
    # Show data source
    if selected_sample == "E14":
        st.write("**Data Source**: E14.5_2_Tangram.h5ad")
    
    # Spatial coordinates info
    if 'spatial_x' in adata.obs.columns and 'spatial_y' in adata.obs.columns:
        st.write("**Spatial Coordinates**: ✅ Available")
        st.write(f"**X Range**: {adata.obs['spatial_x'].min():.1f} - {adata.obs['spatial_x'].max():.1f}")
        st.write(f"**Y Range**: {adata.obs['spatial_y'].min():.1f} - {adata.obs['spatial_y'].max():.1f}")
    else:
        st.write("**Spatial Coordinates**: ❌ Not found")
    
    # Show available obsm keys
    if hasattr(adata, 'obsm') and adata.obsm:
        st.subheader("🗂️ Available Spatial Data")
        st.write("**obsm keys**:")
        for key in adata.obsm.keys():
            st.write(f"  - {key}: {adata.obsm[key].shape}")
        
        # Check for spatial coordinates specifically
        if 'spatial' in adata.obsm:
            st.success("✅ Spatial coordinates found in adata.obsm['spatial']")
        else:
            st.error("❌ No spatial coordinates in adata.obsm['spatial']")
            st.info("Squidpy requires spatial coordinates in adata.obsm['spatial']")
    
    # Show available obs columns
    st.subheader("📊 Available Annotations")
    st.write("**obs columns**:")
    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'object':
            unique_vals = adata.obs[col].nunique()
            st.write(f"  - {col}: {unique_vals} unique values")
        else:
            st.write(f"  - {col}: numeric")
    
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
    if 'total_counts' in adata.var.columns:
        top_10_genes = adata.var['total_counts'].nlargest(10)
        for i, (gene, counts) in enumerate(top_10_genes.items(), 1):
            st.write(f"{i}. {gene}: {counts:,.0f}")
    else:
        st.write("Count data not available")

# Additional spatial analysis options
st.subheader("🔬 Additional Spatial Analysis")

col1, col2 = st.columns(2)

with col1:
    st.markdown("""
    ### Available Analyses
    
    - **Spatial Scatter Plots**: Visualize cell distributions
    - **Neighborhood Enrichment**: Analyze spatial relationships
    - **Gene Expression Mapping**: Spatial gene visualization
    - **Multiple Plot Comparison**: Side-by-side analysis
    - **Spatial Statistics**: Comprehensive spatial metrics
    
    ### Export Options
    
    - **Spatial Coordinates**: Export cell positions
    - **Expression Matrix**: Export gene expression data
    - **Cell Metadata**: Export cell annotations
    - **Plots**: Save visualizations as images
    """)

with col2:
    st.markdown("""
    ### Analysis Parameters
    
    - **Point Size**: Adjust cell point sizes
    - **Color Schemes**: Customize plot colors
    - **Gene Selection**: Select genes for spatial analysis
    - **Cluster Keys**: Choose annotation columns
    
    ### Visualization Settings
    
    - **Plot Layout**: Configure subplot arrangements
    - **Title Customization**: Set plot titles
    - **Error Handling**: Comprehensive error checking
    - **Performance**: Optimized spatial computations
    """)
