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
import tifffile as tiff
from PIL import Image

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
st.sidebar.info(f"Data directory: {data_manager.spatial_data_dir}")

spatial_catalog = data_manager.get_spatial_catalog()
available_samples = data_manager.get_available_samples()

if not available_samples:
    st.sidebar.error("No spatial datasets detected in the catalog")
    st.info("Please add spatial .h5ad files to the data directory and reload the app.")
    st.stop()

sample_labels = [spatial_catalog.get(key, {}).get("display_name", key) for key in available_samples]
st.sidebar.caption("Available samples: " + ", ".join(sample_labels))

selected_sample = st.sidebar.selectbox(
    "Select Sample:",
    available_samples,
    format_func=lambda key: spatial_catalog.get(key, {}).get("display_name", key),
    help="Choose the spatial dataset to analyze"
)

selected_entry = spatial_catalog.get(selected_sample)
if not selected_entry:
    st.sidebar.error("Selected sample is missing from the catalog. Please reload the app.")
    st.stop()

sample_display = selected_entry.get("display_name", selected_sample)

# Data status and loading
st.sidebar.header("📊 Data Status")

# Show current data status
loaded_spatial = data_manager.get_spatial_dataset(selected_sample)
if loaded_spatial is not None:
    st.sidebar.success(f"✅ Spatial\n{loaded_spatial.n_obs} cells")
else:
    st.sidebar.info("⏳ Spatial\nNot loaded")

# Auto-load data if not loaded
if loaded_spatial is None:
    with st.spinner(f"Loading spatial data for {sample_display}..."):
        try:
            loaded_spatial = data_manager.load_spatial_data(selected_sample)
        except Exception as e:
            st.error(f"Error loading data: {e}")
            st.info("Please restart the app to use the updated DataManager")
            st.stop()

if loaded_spatial is None:
    st.error(f"Failed to load spatial data for {sample_display}")
    st.stop()

# Manual reload button
if st.sidebar.button("🔄 Reload Data"):
    with st.spinner(f"Reloading spatial data for {sample_display}..."):
        reloaded = data_manager.load_spatial_data(selected_sample)
        if reloaded is not None:
            st.sidebar.success(f"Reloaded {reloaded.n_obs} cells")
        else:
            st.sidebar.error("Reload failed. Check the dataset and try again.")

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

# Main content
adata = data_manager.get_spatial_dataset(selected_sample)

if adata is None:
    st.error(f"Failed to locate spatial data for {sample_display}")
    st.stop()

st.title(f"🗺️ Spatial Data - {sample_display}")

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

# DAPI image and Spatial Scatter side-by-side
st.subheader("DAPI image and Spatial Scatter")

if SQUIDPY_AVAILABLE:
    st.markdown("##### Plot Settings")

    set_col1, set_col2 = st.columns(2)

    with set_col1:
        # Color options for squidpy
        squidpy_color_options = ['celltype'] + [col for col in adata.obs.columns if col not in ['spatial_x', 'spatial_y', 'celltype']]
        squidpy_color = st.selectbox("Color by:", squidpy_color_options, key="squidpy_color")

    with set_col2:
        # Point size
        point_size = st.slider("Point size:", min_value=1, max_value=50, value=20, key="squidpy_size")

img_col, plot_col = st.columns(2, gap="medium")

with img_col:
    st.markdown("#### Masked DAPI image")

    # Load the matching DAPI image for the selected sample
    image_path = data_manager.get_dapi_image_path(selected_sample)

    if image_path and image_path.exists():
        try:
            arr = tiff.imread(str(image_path))

            # If 16-bit grayscale, normalize to 0–255 for display
            if arr.dtype == np.uint16 and arr.max() > 0:
                arr = (arr.astype(np.float32) / arr.max() * 255).astype(np.uint8)

            st.image(arr, caption=image_path.name, use_container_width=True, clamp=True)
        except Exception as e:
            st.write("tifffile could not read this TIFF:", e)
    else:
        st.info("No masked DAPI image found for this sample.")

with plot_col:
    if SQUIDPY_AVAILABLE:
        st.markdown("#### 🔬 Squidpy Spatial Scatter Plot")

        # Generate plot automatically when data is loaded
        with st.spinner("Generating Squidpy spatial scatter plot..."):
            fig = spatial_viz_manager.plot_spatial_scatter(
                adata,
                color=squidpy_color,
                size=point_size,
                title=f"Squidpy Spatial Scatter - {sample_display}"
            )
            if fig:
                fig.set_size_inches(6, 6)
                st.pyplot(fig, use_container_width=True)
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
    st.markdown("Cluster key fixed to **celltype** for neighborhood enrichment.")
    selected_cluster_key = 'celltype'
    
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
        """)
    
    # Embed the tissuumap webpage
    st.markdown("---")
    st.markdown("### 🗺️ Interactive Tissuumap")
    
    # Embed the tissuumap HTML file
    #with open(tissuumap_index, 'r', encoding='utf-8') as f:
    #    html_content = f.read()
    
    # Display the tissuumap
    #st.components.v1.html(html_content, height=600, scrolling=True)

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
                title=f"Gene Expression in Space - {sample_display}"
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

# Cell Type Proportion Plot
st.subheader("🥧 Cell Type Proportion")

proportion_fig = viz_manager.plot_celltype_proportion(adata, sample_display, title="Cell Type Proportion")
if proportion_fig:
    st.plotly_chart(proportion_fig, use_container_width=True)

# UMAP Plot
# st.subheader("📊 UMAP Plot")

# Grouping selection above the figure
# groupby_options = ['active.ident'] + [col for col in adata.obs.columns if col != 'active.ident']
#selected_groupby = st.selectbox("Group by:", groupby_options, key="spatial_umap_groupby")

#umap_fig = viz_manager.plot_umap(adata, groupby=selected_groupby, title=f"UMAP - {selected_sample}")
#if umap_fig:
#    st.plotly_chart(umap_fig, use_container_width=True)


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
