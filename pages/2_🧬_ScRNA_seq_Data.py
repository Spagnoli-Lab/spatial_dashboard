#!/usr/bin/env python3
"""
ScRNA-seq Data Analysis Page
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
    page_title="ScRNA-seq Data - Spatial Transcriptomics Dashboard",
    page_icon="🧬",
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
st.sidebar.title("🧬 ScRNA-seq Analysis")

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
if selected_sample in data_manager.scrna_data:
    st.sidebar.success(f"✅ scRNA-seq\n{data_manager.scrna_data[selected_sample].n_obs} cells")
else:
    st.sidebar.info("⏳ scRNA-seq\nNot loaded")

# Auto-load data if not loaded
if selected_sample not in data_manager.scrna_data:
    with st.spinner(f"Loading scRNA-seq data for {selected_sample}..."):
        data_manager.load_scrna_data(selected_sample)

# Manual reload button
if st.sidebar.button("🔄 Reload Data"):
    with st.spinner(f"Reloading scRNA-seq data for {selected_sample}..."):
        data_manager.load_scrna_data(selected_sample)
        if selected_sample in data_manager.scrna_data:
            st.sidebar.success(f"Reloaded {data_manager.scrna_data[selected_sample].n_obs} cells")

# Check if data is available
if selected_sample not in data_manager.scrna_data:
    st.error(f"Failed to load scRNA-seq data for {selected_sample}")
    st.info("Please check that scRNA-seq data files exist for this sample")
    st.stop()

# Main content
adata = data_manager.scrna_data[selected_sample]

st.title(f"🧬 ScRNA-seq Data - {selected_sample}")

# Summary statistics
col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric("Cells", f"{adata.n_obs:,}")
with col2:
    st.metric("Genes", f"{adata.n_vars:,}")
with col3:
    st.metric("Avg Counts/Cell", f"{np.mean(adata.obs['total_counts']):.1f}")
with col4:
    st.metric("Avg Genes/Cell", f"{np.mean(adata.obs['n_genes_by_counts']):.1f}")

# UMAP Plot
st.subheader("📊 UMAP Plot")

col1, col2 = st.columns([3, 1])

with col1:
    # Get available grouping variables
    groupby_options = ['active.ident'] + [col for col in adata.obs.columns if col != 'active.ident']
    selected_groupby = st.selectbox("Group by:", groupby_options, key="umap_groupby")
    
    umap_fig = viz_manager.plot_umap(adata, groupby=selected_groupby, title=f"UMAP - {selected_sample}")
    if umap_fig:
        st.plotly_chart(umap_fig, use_container_width=True)

# Dot Plot
st.subheader("🔴 Dot Plot")

col1, col2 = st.columns([3, 1])

with col1:
    # Get available grouping variables
    groupby_options = ['active.ident'] + [col for col in adata.obs.columns if col != 'active.ident']
    selected_groupby = st.selectbox("Group by:", groupby_options, key="dot_groupby")
    
    # Gene selection
    top_genes = adata.var['total_counts'].nlargest(20).index.tolist()
    selected_genes = st.multiselect("Select genes:", top_genes, default=top_genes[:5])
    
    dot_fig = viz_manager.plot_dotplot(adata, groupby=selected_groupby, genes=selected_genes, title=f"Dot Plot - {selected_sample}")
    if dot_fig:
        st.plotly_chart(dot_fig, use_container_width=True)

# Feature Plot
st.subheader("🎨 Feature Plot")

col1, col2 = st.columns([3, 1])

with col1:
    # Gene selection
    top_genes = adata.var['total_counts'].nlargest(20).index.tolist()
    selected_genes = st.multiselect("Select genes:", top_genes, default=top_genes[:3], key="feature_genes")
    
    feature_fig = viz_manager.plot_feature_plot(adata, genes=selected_genes, title=f"Feature Plot - {selected_sample}")
    if feature_fig:
        st.plotly_chart(feature_fig, use_container_width=True)

# Cell Type Proportion Plot
st.subheader("🥧 Cell Type Proportion")

col1, col2 = st.columns([3, 1])

with col1:
    sample_options = ['E12', 'E14', 'E17']
    selected_sample_day = st.selectbox("Select sample day:", sample_options, index=sample_options.index(selected_sample))
    
    proportion_fig = viz_manager.plot_celltype_proportion(adata, selected_sample_day, title="Cell Type Proportion")
    if proportion_fig:
        st.plotly_chart(proportion_fig, use_container_width=True)

# Data information sidebar
with col2:
    st.subheader("📋 Data Info")
    st.write(f"**Sample**: {selected_sample}")
    st.write(f"**Cells**: {adata.n_obs:,}")
    st.write(f"**Genes**: {adata.n_vars:,}")
    
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
    
    # Top genes
    st.subheader("🔝 Top Expressed Genes")
    top_10_genes = adata.var['total_counts'].nlargest(10)
    for i, (gene, counts) in enumerate(top_10_genes.items(), 1):
        st.write(f"{i}. {gene}: {counts:,.0f}")
