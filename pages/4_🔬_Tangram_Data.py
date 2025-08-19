#!/usr/bin/env python3
"""
Tangram Data Analysis Page with Tissuumaps Integration
"""

import streamlit as st
import pandas as pd
import numpy as np
import subprocess
import time
import webbrowser
from pathlib import Path
import sys

# Add the project root to the Python path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

# Import utilities
from utils.data_manager import DataManager
from utils.tissuumaps_manager import TissuumapsManager

st.set_page_config(
    page_title="Tangram Data - Spatial Transcriptomics Dashboard",
    page_icon="🔬",
    layout="wide"
)

# Initialize managers
@st.cache_resource
def get_managers():
    data_manager = DataManager()
    tissuumaps_manager = TissuumapsManager(data_manager.data_dir)
    return data_manager, tissuumaps_manager

data_manager, tissuumaps_manager = get_managers()

# Sidebar
st.sidebar.title("🔬 Tangram Data Analysis")

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
if selected_sample in data_manager.tangram_data:
    st.sidebar.success(f"✅ Tangram\n{data_manager.tangram_data[selected_sample].n_obs} cells")
else:
    st.sidebar.info("⏳ Tangram\nNot loaded")

# Auto-load data if not loaded
if selected_sample not in data_manager.tangram_data:
    with st.spinner(f"Loading Tangram data for {selected_sample}..."):
        data_manager.load_tangram_data(selected_sample)

# Manual reload button
if st.sidebar.button("🔄 Reload Data"):
    with st.spinner(f"Reloading Tangram data for {selected_sample}..."):
        data_manager.load_tangram_data(selected_sample)
        if selected_sample in data_manager.tangram_data:
            st.sidebar.success(f"Reloaded {data_manager.tangram_data[selected_sample].n_obs} cells")

# Check if data is available
if selected_sample not in data_manager.tangram_data:
    st.error(f"Failed to load Tangram data for {selected_sample}")
    st.info("Please check that Tangram data files exist for this sample")
    st.stop()

# Main content
adata = data_manager.tangram_data[selected_sample]

st.title(f"🔬 Tangram Data - {selected_sample}")

# Summary statistics
col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric("Cells", f"{adata.n_obs:,}")
with col2:
    st.metric("Genes", f"{adata.n_vars:,}")
with col3:
    st.metric("Avg Counts/Cell", f"{np.mean(adata.obs['total_counts']):.1f}")
with col4:
    st.metric("Total Counts", f"{np.sum(adata.X):,.0f}")

# Tissuumaps Integration
st.subheader("🌐 Tangram Object with DAPI Image")

# Server controls
col1, col2, col3 = st.columns(3)

with col1:
    if st.button("🚀 Start Tissuumaps Server", key="start_server"):
        with st.spinner("Starting Tissuumaps server..."):
            if tissuumaps_manager.start_tissuumaps_server():
                st.success("✅ Tissuumaps server started successfully!")
            else:
                st.error("❌ Failed to start server")

with col2:
    if st.button("🛑 Stop Server", key="stop_server"):
        tissuumaps_manager.stop_server()
        st.success("✅ Server stopped")

with col3:
    tissuumaps_url = tissuumaps_manager.get_tissuumaps_url()
    if st.button("🌐 Open in New Tab", key="open_browser"):
        try:
            webbrowser.open(tissuumaps_url)
            st.success("✅ Opened in browser")
        except Exception as e:
            st.error(f"❌ Failed to open browser: {e}")

# Display Tissuumaps URL
st.info(f"**Tissuumaps URL**: [{tissuumaps_url}]({tissuumaps_url})")

# Embed Tissuumaps in iframe
st.subheader("📊 Interactive Tissuumaps Viewer")

try:
    # Create iframe with Tissuumaps
    iframe_html = f"""
    <iframe 
        src="{tissuumaps_url}" 
        width="100%" 
        height="800" 
        frameborder="0"
        style="border: 2px solid #ddd; border-radius: 10px;"
    ></iframe>
    """
    st.components.v1.html(iframe_html, height=800)
except Exception as e:
    st.error(f"Failed to embed Tissuumaps: {e}")
    st.info("Please start the Tissuumaps server first using the button above.")

# Instructions
st.subheader("📖 Tissuumaps Instructions")

col1, col2 = st.columns(2)

with col1:
    st.markdown(f"""
    ### 🎯 Getting Started
    
    1. **Start the server** using the "Start Tissuumaps Server" button
    2. **Wait for server** to initialize (may take a few seconds)
    3. **Load your data** in Tissuumaps:
       - Click "Load Data" in the Tissuumaps interface
       - Select your DAPI image file (`.tif`) for {selected_sample}
       - Select the Tangram data file for {selected_sample}
    
    ### 🔧 Configuration
    
    - **X and Y selectors**: Set to spatial coordinates
    - **Color selector**: Choose cell type or gene expression
    - **Size selector**: Set to cell size or expression level
    - **Opacity**: Adjust for better visualization
    """)

with col2:
    st.markdown(f"""
    ### 🎨 Visualization Tips
    
    - **Zoom and Pan**: Use mouse to navigate the spatial view
    - **Layer Management**: Toggle different data layers
    - **Color Schemes**: Customize color palettes
    - **Export Options**: Save images and data
    
    ### 🔬 Analysis Features
    
    - **Spatial Clustering**: Identify spatial domains
    - **Gene Expression**: Map gene expression in space
    - **Cell Type Distribution**: Visualize cell type patterns
    - **Quantitative Analysis**: Measure spatial statistics
    """)

# Data preview
st.subheader("📋 Tangram Data Preview")

col1, col2 = st.columns(2)

with col1:
    st.write("**Sample of Tangram data:**")
    st.dataframe(adata.obs.head(10))
    
    # Show available columns
    st.write("**Available metadata columns:**")
    st.write(list(adata.obs.columns))

with col2:
    # Data statistics
    st.write("**Data Statistics:**")
    st.write(f"- **Shape**: {adata.shape}")
    st.write(f"- **Sparsity**: {(1 - np.count_nonzero(adata.X) / adata.X.size) * 100:.1f}%")
    st.write(f"- **Memory usage**: {adata.nbytes / 1024 / 1024:.1f} MB")
    
    # Quality metrics
    st.write("**Quality Metrics:**")
    st.write(f"- **Mean counts/cell**: {np.mean(adata.obs['total_counts']):.1f}")
    st.write(f"- **Median counts/cell**: {np.median(adata.obs['total_counts']):.1f}")
    st.write(f"- **Mean genes/cell**: {np.mean(adata.obs['n_genes_by_counts']):.1f}")
    st.write(f"- **Median genes/cell**: {np.median(adata.obs['n_genes_by_counts']):.1f}")

# File information
st.subheader("📁 File Information")

# Check for DAPI images
data_dir = Path(data_manager.data_dir)
dapi_files = list(data_dir.glob(f"*{selected_sample}*DAPI*.tif"))

if dapi_files:
    st.success(f"✅ Found {len(dapi_files)} DAPI image(s) for {selected_sample}:")
    for dapi_file in dapi_files:
        st.write(f"  - {dapi_file.name}")
else:
    st.warning(f"⚠️ No DAPI images found for {selected_sample}")
    st.info("Please ensure DAPI images follow the naming convention: `{sample}_DAPI.tif`")

# Check for Tangram data files
tangram_files = list(data_dir.glob(f"*{selected_sample}*tangram*.h5ad"))
if not tangram_files:
    tangram_files = list(data_dir.glob(f"*{selected_sample}*Tangram*.h5ad"))

if tangram_files:
    st.success(f"✅ Found {len(tangram_files)} Tangram data file(s) for {selected_sample}:")
    for tangram_file in tangram_files:
        st.write(f"  - {tangram_file.name}")
else:
    st.warning(f"⚠️ No Tangram data files found for {selected_sample}")

# Troubleshooting
st.subheader("🔧 Troubleshooting")

col1, col2 = st.columns(2)

with col1:
    st.markdown("""
    ### Common Issues
    
    **Server won't start:**
    - Check if port 5100 is available
    - Ensure Tissuumaps is installed
    - Check firewall settings
    
    **Data won't load:**
    - Verify file formats (.h5ad, .tif)
    - Check file permissions
    - Ensure data is in correct directory
    """)

with col2:
    st.markdown("""
    ### Performance Tips
    
    **For large datasets:**
    - Use data subsetting
    - Reduce image resolution
    - Close other applications
    
    **For better visualization:**
    - Adjust point sizes
    - Use appropriate color schemes
    - Enable hardware acceleration
    """)

# Footer
st.markdown("---")
st.markdown("*Tissuumaps integration powered by the Tissuumaps web viewer*")
