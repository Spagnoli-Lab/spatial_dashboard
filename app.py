#!/usr/bin/env python3
"""
Main Application for Spatial Transcriptomics Dashboard
Uses Streamlit's native multi-page navigation
"""

import streamlit as st

# Page configuration
st.set_page_config(
    page_title="Spatial Transcriptomics Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Hide default navigation elements
st.markdown("""
<style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    header {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# Main content - this will be the home page
st.title("🧬 Spatial Transcriptomics Dashboard")
st.markdown("---")

# Welcome message
st.markdown("""
## 🎯 Welcome to the Spatial Transcriptomics Dashboard

This comprehensive dashboard provides an integrated platform for analyzing and visualizing 
spatial transcriptomics data across multiple developmental stages (E12, E14, E17).

### 🚀 Quick Start

1. **Select a page** from the sidebar navigation above
2. **Choose your sample** (E12, E14, or E17) 
3. **Explore the data** with interactive visualizations

### 📊 Available Pages

- **🏠 Home**: Introduction and overview
- **🧬 ScRNA-seq Data**: Single-cell RNA sequencing analysis
- **🗺️ Spatial Data**: Spatial transcriptomics analysis  
- **🔬 Tangram Data**: Tangram mapping results with Tissuumaps

### 📁 Data Requirements

Ensure your data files follow the naming convention:
- `E12_scrna_data.h5ad` for scRNA-seq data
- `E12_spatial_data.h5ad` for spatial data  
- `E12_tangram_data.h5ad` for Tangram results
- `E12_DAPI.tif` for DAPI images

### 🔧 Features

- **Automatic Data Loading**: Data loads automatically when you select a page
- **Interactive Visualizations**: UMAP, dot plots, feature plots, and more
- **Tissuumaps Integration**: Embedded spatial visualization
- **Sample Comparison**: Analyze E12, E14, and E17 developmental stages

---

*Use the sidebar navigation to explore different analysis pages*
""")

# Footer
st.markdown("---")
st.markdown("*Spatial Transcriptomics Dashboard - Powered by Streamlit, Scanpy, and Tissuumaps*")
