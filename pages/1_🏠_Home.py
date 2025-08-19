#!/usr/bin/env python3
"""
Home page for the Spatial Transcriptomics Dashboard
"""

import streamlit as st

st.set_page_config(
    page_title="Home - Spatial Transcriptomics Dashboard",
    page_icon="🏠",
    layout="wide"
)

st.title("🧬 Spatial Transcriptomics Dashboard")
st.markdown("---")

# Introduction
st.header("📖 Welcome to the Spatial Transcriptomics Dashboard")

st.markdown("""
This comprehensive dashboard provides an integrated platform for analyzing and visualizing 
spatial transcriptomics data across multiple developmental stages (E12, E14, E17). 
The dashboard combines single-cell RNA sequencing analysis, spatial data visualization, 
and Tangram mapping results in an interactive web interface.
""")

# Features overview
st.header("🚀 Key Features")

col1, col2, col3 = st.columns(3)

with col1:
    st.markdown("""
    ### 🧬 ScRNA-seq Analysis
    - **UMAP Visualization**: Interactive dimensionality reduction plots
    - **Dot Plots**: Gene expression across cell types
    - **Feature Plots**: Gene expression on UMAP coordinates
    - **Cell Type Proportions**: Composition analysis across samples
    """)

with col2:
    st.markdown("""
    ### 🗺️ Spatial Data Analysis
    - **Spatial Scatter Plots**: Cell distribution in tissue space
    - **Neighborhood Enrichment**: Spatial interaction analysis
    - **UMAP Integration**: Combined spatial and transcriptomic views
    - **Cell Type Mapping**: Spatial distribution of cell types
    """)

with col3:
    st.markdown("""
    ### 🔬 Tangram Results
    - **Tissuumaps Integration**: Interactive spatial visualization
    - **DAPI Image Overlay**: Tangram mapping with tissue images
    - **Sample-specific Analysis**: E12, E14, E17 developmental stages
    - **Real-time Exploration**: Dynamic data exploration
    """)

# Sample information
st.header("📊 Available Samples")

st.markdown("""
The dashboard supports analysis of three developmental stages:

| Stage | Description | Key Features |
|-------|-------------|--------------|
| **E12** | Early embryonic development | Initial cell type specification |
| **E14** | Mid-embryonic development | Tissue patterning and growth |
| **E17** | Late embryonic development | Mature tissue organization |
""")

# Data workflow
st.header("🔄 Analysis Workflow")

st.markdown("""
1. **Sample Selection**: Choose your developmental stage (E12, E14, E17)
2. **Data Loading**: Load scRNA-seq, spatial, or Tangram data
3. **Visualization**: Explore interactive plots and analyses
4. **Integration**: Compare results across data types and samples
""")

# Technical information
st.header("⚙️ Technical Details")

col1, col2 = st.columns(2)

with col1:
    st.markdown("""
    ### Data Formats
    - **H5AD files**: AnnData objects for single-cell data
    - **TIF images**: DAPI staining for spatial reference
    - **CSV files**: Metadata and expression matrices
    
    ### Analysis Tools
    - **Scanpy**: Single-cell analysis pipeline
    - **Tissuumaps**: Spatial data visualization
    - **Tangram**: Spatial mapping algorithms
    """)

with col2:
    st.markdown("""
    ### Visualization Libraries
    - **Plotly**: Interactive plots and dashboards
    - **Matplotlib**: Static scientific plots
    - **Seaborn**: Statistical data visualization
    
    ### Platform
    - **Streamlit**: Web application framework
    - **Python**: Scientific computing environment
    """)

# Getting started
st.header("🎯 Getting Started")

st.markdown("""
### Quick Start Guide

1. **Navigate to Analysis Pages**: Use the sidebar to select your analysis type
2. **Select Sample**: Choose E12, E14, or E17 from the sample dropdown
3. **Load Data**: Click the appropriate "Load" button for your data type
4. **Explore Visualizations**: Use interactive controls to customize your analysis
5. **Compare Results**: Switch between samples and data types for comprehensive analysis

### Tips for Best Experience

- **Start with scRNA-seq**: Load single-cell data first for baseline analysis
- **Use Tissuumaps**: Explore spatial relationships with the embedded viewer
- **Compare Stages**: Analyze differences across developmental timepoints
- **Save Results**: Export plots and data for further analysis
""")

# Contact and support
st.header("📞 Support and Documentation")

st.markdown("""
For technical support, documentation, or feature requests, please refer to:

- **Documentation**: Check the README files in the project directory
- **Issues**: Report bugs or request features through the project repository
- **Data Format**: Ensure your data follows the expected H5AD format
- **Dependencies**: Verify all required packages are installed

### Data Requirements

Ensure your data files follow the naming convention:
- `E12_scrna_data.h5ad` for scRNA-seq data
- `E12_spatial_data.h5ad` for spatial data  
- `E12_tangram_data.h5ad` for Tangram results
- `E12_DAPI.tif` for DAPI images
""")

# Footer
st.markdown("---")
st.markdown("*Spatial Transcriptomics Dashboard - Powered by Streamlit, Scanpy, and Tissuumaps*")
