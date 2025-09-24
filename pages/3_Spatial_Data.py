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
from io import BytesIO

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
    st.sidebar.success(f"Squidpy version: {squidpy.__version__}")
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

# Sample selection
st.sidebar.header("📁 Sample Selection")

registered_catalog = data_manager.get_registered_catalog()
available_samples = data_manager.get_available_registered_samples()

def _clean_sample_label(label: str) -> str:
    """Return a concise version of the sample label without metadata suffixes."""
    if not label:
        return ""
    return label.split("(", 1)[0].strip()


def _get_sample_label(sample_key: str) -> str:
    entry = registered_catalog.get(sample_key, {})
    candidates = [
        entry.get("primary_sample"),
        entry.get("key"),
        entry.get("file_stem"),
        sample_key,
    ]
    for candidate in candidates:
        cleaned = _clean_sample_label(str(candidate)) if candidate else ""
        if cleaned:
            return cleaned
    return sample_key

if not available_samples:
    st.sidebar.error("No spatial datasets detected in the catalog")
    st.info("Please add spatial .h5ad files to the data directory and reload the app.")
    st.stop()

sample_labels = [_get_sample_label(key) for key in available_samples]

selected_sample = st.sidebar.selectbox(
    "Select Sample:",
    available_samples,
    format_func=_get_sample_label,
    help="Choose the spatial dataset to analyze"
)

selected_entry = registered_catalog.get(selected_sample)
if not selected_entry:
    st.sidebar.error("Selected sample is missing from the catalog. Please reload the app.")
    st.stop()

sample_display = _get_sample_label(selected_sample)

# Data status and loading
st.sidebar.header("📊 Data Status")

# Show current data status
loaded_spatial = data_manager.get_registered_dataset(selected_sample, finalized=True)
#if loaded_spatial is not None:
#    st.sidebar.success(f"✅ Spatial\n{loaded_spatial.n_obs} cells")
#else:
#    st.sidebar.info(" Spatial\nNot loaded")

# Auto-load data if not loaded
if loaded_spatial is None:
    with st.spinner(f"Loading spatial data for {sample_display}..."):
        try:
           loaded_spatial = data_manager.load_registered_data(selected_sample, finalize=True)
        except Exception as e:
            st.error(f"Error loading data: {e}")
            st.info("Please restart the app to use the updated DataManager")
            st.stop()

if loaded_spatial is None:
    st.error(f"Failed to load spatial data for {sample_display}")
    st.stop()

# Manual reload button
#if st.sidebar.button("🔄 Reload Data"):
#    with st.spinner(f"Reloading spatial data for {sample_display}..."):
#        reloaded = data_manager.load_registered_data(selected_sample, finalize=True)
#        if reloaded is not None:
#            st.sidebar.success(f"Reloaded {reloaded.n_obs} cells")
#        else:
#            st.sidebar.error("Reload failed. Check the dataset and try again.")

# Available Spatial Visualizations
if SQUIDPY_AVAILABLE:
    st.sidebar.header("🔬 Available Spatial Visualizations")
    st.sidebar.markdown("""
    - **Spatial Scatter Plots**: Color by annotations
    - **Neighborhood Enrichment**: Spatial relationships
    - **Gene Expression**: Spatial gene mapping
    """)

# Main content
adata = data_manager.get_registered_dataset(selected_sample, finalized=True)

if adata is None:
    st.error(f"Failed to locate spatial data for {sample_display}")
    st.stop()

st.title(f"🗺️ Spatial Data - {sample_display}")

# Summary statistics
#col1, col2, col3, col4 = st.columns(4)

#with col1:
#    st.metric("Cells", f"{adata.n_obs:,}")
#with col2:
#    st.metric("Genes", f"{adata.n_vars:,}")
#with col3:
#    st.metric("Avg Counts/Cell", f"{np.mean(adata.obs['total_counts']):.1f}")
#with col4:
#    if 'spatial_x' in adata.obs.columns:
#st.metric("Spatial Coords", "Available")
#    else:
#        st.metric("Spatial Coords", "Not found")

# DAPI image and Spatial Scatter side-by-side
st.subheader("DAPI image and Spatial Scatter")

if SQUIDPY_AVAILABLE:
    image_path = data_manager.get_dapi_image_path(selected_sample)
    dapi_image = None
    dapi_read_error = None
    if image_path and image_path.exists():
        try:
            dapi_image = tiff.imread(str(image_path))

            # If 16-bit grayscale, normalize to 0–255 for display
            if dapi_image.dtype == np.uint16 and dapi_image.max() > 0:
                dapi_image = (dapi_image.astype(np.float32) / dapi_image.max() * 255).astype(np.uint8)
        except Exception as exc:
            dapi_read_error = str(exc)

    dapi_width_key = f"dapi_width_{selected_sample}"
    dapi_default_width = 350
    current_width = int(max(100, min(st.session_state.get(dapi_width_key, dapi_default_width), 2000)))

    img_col, plot_col, setting_col = st.columns(3, gap="medium")
    
    with setting_col:
        st.markdown("#### Plot Settings")

        if dapi_image is not None:
            current_width = st.slider(
                "DAPI image width (px)",
                min_value=100,
                max_value=600,
                value=current_width,
                step=10,
                key=dapi_width_key
            )
        elif dapi_read_error:
            st.info("DAPI image detected but could not be read.")
        else:
            st.info("No DAPI image available for sizing.")

        # Color options for squidpy
        squidpy_color_options = ['celltype'] + [col for col in adata.obs.columns if col not in ['spatial_x', 'spatial_y', 'celltype']]
        squidpy_color = st.selectbox("Color by:", squidpy_color_options, key="squidpy_color")

        # Point size
        point_size = st.slider("Point size:", min_value=1, max_value=50, value=20, key="squidpy_size")

        # Figure dimensions in inches
        col_w, col_h = st.columns(2)
        with col_w:
            scatter_width = st.slider(
                "Plot width (inches)",
                min_value=3.0,
                max_value=12.0,
                value=6.0,
                step=0.5,
                key="squidpy_width"
            )
        with col_h:
            scatter_height = st.slider(
                "Plot height (inches)",
                min_value=3.0,
                max_value=12.0,
                value=8.0,
                step=0.5,
                key="squidpy_height"
            )
            
    with img_col:
        st.markdown("#### Masked DAPI image")

        if dapi_image is not None:
            display_width = st.session_state.get(dapi_width_key, current_width)
            display_width = int(max(100, min(display_width, 2000)))
            st.image(dapi_image, caption=image_path.name, width=display_width, clamp=True)
        elif image_path and image_path.exists():
            if dapi_read_error:
                st.info(f"Unable to display DAPI image: {dapi_read_error}")
            else:
                st.info("tifffile could not read this TIFF.")
        else:
            st.info("No masked DAPI image found for this sample.")
    
    
    with plot_col:
        st.markdown("#### Squidpy Spatial Scatter Plot")

        # Generate plot automatically when data is loaded
        with st.spinner("Generating Squidpy spatial scatter plot..."):
            fig = spatial_viz_manager.plot_spatial_scatter(
                adata,
                color=squidpy_color,
                size=point_size,
                title=f"Squidpy Spatial Scatter - {sample_display}",
                figsize=(scatter_width, scatter_height)
            )
            if fig:
                buffer = BytesIO()
                fig.savefig(buffer, format="png", dpi=fig.dpi, bbox_inches="tight")
                buffer.seek(0)
                st.image(buffer, caption="Spatial scatter", clamp=True)
                plt.close(fig)
            else:
                st.warning("Could not generate spatial scatter plot. Check if spatial coordinates are available.")
else:
    st.info("💡 Install squidpy to enable advanced spatial analysis features: `pip install squidpy`")



# Neighborhood Enrichment Analysis
st.subheader("🔬 Neighborhood Enrichment Analysis")

# Squidpy Neighborhood Enrichment Analysis
if SQUIDPY_AVAILABLE:

    # Set the cluster key
    selected_cluster_key = 'celltype'

    #col1, col2 = st.columns(2)
    #with col1:
    #    ne_width = st.slider(
    #        "Plot width (inches)",
    #        min_value=1.0,
    #        max_value=12.0,
    #        value=8.0,
    #        step=0.5,
    #        key="ne_plot_width"
    #    )
    #    ne_height = st.slider(
    #        "Plot height (inches)",
    #        min_value=1.0,
    #        max_value=12.0,
    #        value=8.0,
    #        step=0.5,
    #        key="ne_plot_height"
    #    )

    # Generate analysis automatically when data is loaded
    with st.spinner("Computing neighborhood enrichment..."):
        fig = spatial_viz_manager.plot_neighborhood_enrichment(
            adata,
            cluster_key=selected_cluster_key,
            figsize=(8,8)
        )
        if fig:
            buffer = BytesIO()
            fig.savefig(buffer, format="png", dpi=fig.dpi, bbox_inches="tight")
            buffer.seek(0)
            st.image(buffer, clamp=True)
            plt.close(fig)
            st.success("✅ Neighborhood enrichment computed")
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

    gene_col_w, gene_col_h = st.columns(2)
    with gene_col_w:
        gene_plot_width = st.slider(
            "Plot width (inches)",
            min_value=3.0,
            max_value=12.0,
            value=6.0,
            step=0.5,
            key="gene_plot_width"
        )
    with gene_col_h:
        gene_plot_height = st.slider(
            "Plot height (inches)",
            min_value=3.0,
            max_value=12.0,
            value=6.0,
            step=0.5,
            key="gene_plot_height"
        )
    
    # Generate plot automatically when genes are selected
    if selected_genes:
        with st.spinner("Generating gene expression spatial plot..."):
            fig = spatial_viz_manager.plot_gene_expression_spatial(
                adata,
                genes=selected_genes,
                size=gene_point_size,
                title=f"Gene Expression in Space - {sample_display}",
                figsize=(gene_plot_width, gene_plot_height)
            )
            if fig:
                buffer = BytesIO()
                fig.savefig(buffer, format="png", dpi=fig.dpi, bbox_inches="tight")
                buffer.seek(0)
                st.image(buffer, caption="Gene expression spatial plot", clamp=True)
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
