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
import matplotlib.pyplot as plt

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
available_samples = data_manager.get_scrna_sample_options()

if not available_samples:
    st.sidebar.error("No sample data found in data directory")
    st.info("Please ensure your data files contain sample identifiers")
    st.stop()

# Sample selection with checkboxes
st.sidebar.subheader("Selected Samples:")

# Quick selection buttons
if st.button("✅ Select All"):
    st.session_state.selected_samples = available_samples.copy()
    st.rerun()

if st.button("❌ Deselect All"):
    st.session_state.selected_samples = []
    st.rerun()

# Initialize selected_samples in session state if not exists
if 'selected_samples' not in st.session_state:
    st.session_state.selected_samples = []

# Individual sample checkboxes
#st.sidebar.write("**Individual Samples:**")
for sample in available_samples:
    if st.checkbox(
        sample, 
        value=sample in st.session_state.selected_samples,
        key=f"sample_{sample}"
    ):
        if sample not in st.session_state.selected_samples:
            st.session_state.selected_samples.append(sample)
    else:
        if sample in st.session_state.selected_samples:
            st.session_state.selected_samples.remove(sample)

# Get the selected samples
selected_samples = st.session_state.selected_samples

# Show selected samples
if selected_samples:
    st.sidebar.success(f"✅ {', '.join(selected_samples)}")
else:
    st.sidebar.warning("⚠️ No samples selected")
    st.info("Please select at least one sample to continue")
    st.stop()

# Data status and loading
st.sidebar.header("📊 Data Status")

# Load data for all selected samples
loaded_samples = []
total_cells = 0

for sample in selected_samples:
    if sample not in data_manager.scrna_data:
        with st.spinner(f"Loading scRNA-seq data for {sample}..."):
            data_manager.load_scrna_data(sample)
    
    if sample in data_manager.scrna_data:
        loaded_samples.append(sample)
        total_cells += data_manager.scrna_data[sample].n_obs

# Show data status
if loaded_samples:
    st.sidebar.success(f"✅ scRNA-seq\n{len(loaded_samples)} samples, {total_cells:,} total cells")
else:
    st.sidebar.error("❌ No data loaded")

# Show information about the data source
st.sidebar.info(f"📁 Data source: Cartana.h5ad\n🔍 Filtered by: orig.ident")

# Manual reload button
if st.sidebar.button("🔄 Reload All Data"):
    with st.spinner("Reloading all selected samples..."):
        for sample in selected_samples:
            data_manager.load_scrna_data(sample)
        st.rerun()

# Check if data is available
if not loaded_samples:
    st.error("Failed to load scRNA-seq data for selected samples")
    st.info("Please check that scRNA-seq data files exist for the selected samples")
    st.stop()

# Combine data from all selected samples
if len(loaded_samples) == 1:
    # Single sample - use directly
    adata = data_manager.scrna_data[loaded_samples[0]]
    sample_display_name = loaded_samples[0]
else:
    # Multiple samples - combine them
    import anndata
    adata_list = [data_manager.scrna_data[sample] for sample in loaded_samples]
    adata = anndata.concat(adata_list, join='outer', index_unique=None)
    sample_display_name = f"{len(loaded_samples)} samples ({', '.join(loaded_samples)})"

st.title(f"🧬 ScRNA-seq Data - {sample_display_name}")

# Summary statistics
col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric("Cells", f"{adata.n_obs:,}")
with col2:
    st.metric("Genes", f"{adata.n_vars:,}")
with col3:
    if 'total_counts' in adata.obs.columns:
        st.metric("Avg Counts/Cell", f"{np.mean(adata.obs['total_counts']):.1f}")
    else:
        st.metric("Avg Counts/Cell", "N/A")
with col4:
    if 'n_genes_by_counts' in adata.obs.columns:
        st.metric("Avg Genes/Cell", f"{np.mean(adata.obs['n_genes_by_counts']):.1f}")
    else:
        st.metric("Avg Genes/Cell", "N/A")

# Debug information (can be removed later)
with st.expander("🔍 Debug Info"):
    st.write("**Available obs columns:**", list(adata.obs.columns))
    st.write("**Available var columns:**", list(adata.var.columns))
    st.write("**Data shape:**", adata.shape)
    st.write("**Data type:**", type(adata.X))

# UMAP Plot
st.subheader("📊 UMAP Plot")

# Get available grouping variables
groupby_options = ['active.ident'] + [col for col in adata.obs.columns if col != 'active.ident']
selected_groupby = st.selectbox("Group by:", groupby_options, key="umap_groupby")

umap_fig = viz_manager.plot_umap(adata, groupby=selected_groupby, title=f"UMAP - {sample_display_name}")
if umap_fig:
    st.plotly_chart(umap_fig, use_container_width=True)

# Dot Plot
st.subheader("🔴 Dot Plot")

# Get available grouping variables
groupby_options = ['active.ident'] + [col for col in adata.obs.columns if col != 'active.ident']
selected_groupby = st.selectbox("Group by:", groupby_options, key="dot_groupby")

# Get available genes (all genes in the dataset)
all_genes = adata.var_names.tolist()

# Show top genes by expression if available, otherwise first 100 genes
if 'total_counts' in adata.var.columns:
    display_genes = adata.var['total_counts'].nlargest(100).index.tolist()
else:
    display_genes = all_genes[:100]

# Use multiselect for gene selection
gene_list = st.multiselect(
    "Select genes to include in dot plot:",
    options=display_genes,
    default=display_genes[:5] if len(display_genes) >= 5 else display_genes,
    help="Select the genes you want to include in the dot plot"
)

# Show selected genes count
if gene_list:
    st.success(f"✅ Selected {len(gene_list)} genes for dot plot")
else:
    st.warning("⚠️ No genes selected for dot plot")

# Create dot plot using scanpy
if gene_list and selected_groupby in adata.obs.columns:
    
    # Import scanpy for plotting
    import scanpy as sc
    
    # Create the dot plot
    fig, ax = plt.subplots(figsize=(12, 8))
    sc.pl.dotplot(adata, gene_list, groupby=selected_groupby, ax=ax, show=False)
    
    # Display the plot
    st.pyplot(fig)
    plt.close(fig)
    
    # Show plot info
    st.info(f"Dot plot showing {len(gene_list)} genes grouped by '{selected_groupby}'")
elif not gene_list:
    st.info("Please select genes to create the dot plot")
else:
    st.error(f"Grouping variable '{selected_groupby}' not found in data")

# Feature Plot
st.subheader("🎨 Feature Plot")

# Gene selection - check if total_counts exists, otherwise use gene names
if 'total_counts' in adata.var.columns:
    top_genes = adata.var['total_counts'].nlargest(20).index.tolist()
else:
    # Fallback: use first 20 genes by name
    top_genes = adata.var_names[:20].tolist()

selected_genes = st.multiselect("Select genes:", top_genes, default=top_genes[:3], key="feature_genes")

feature_fig = viz_manager.plot_feature_plot(adata, genes=selected_genes, title=f"Feature Plot - {sample_display_name}")
if feature_fig:
    st.plotly_chart(feature_fig, use_container_width=True)

# Cell Type Proportion Plot
st.subheader("🧫 Cell Type Proportion")

# Use fixed grouping variables
selected_cell_type = 'active.ident'  # Fixed cell type grouping
selected_sample_group = 'orig.ident'  # This represents the samples users selected

# Create the stacked bar chart
if selected_cell_type in adata.obs.columns and selected_sample_group in adata.obs.columns:
    
    # Create the crosstab and stacked bar plot
    fig, ax = plt.subplots(figsize=(10, 6))
    tmp = pd.crosstab(adata.obs[selected_cell_type], adata.obs[selected_sample_group], normalize='columns').T.plot(kind='bar', stacked=True, ax=ax)
    
    # Adjust the legend position
    tmp.legend(title= "Cell Type", bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Set labels and title
    #plt.xlabel(selected_sample_group)
    plt.ylabel('Proportion')
    plt.title(f'Cell Type Proportion - {sample_display_name}')
    
    # Adjust layout to prevent legend overlap
    plt.tight_layout()
    
    # Display the plot
    st.pyplot(fig)
    plt.close(fig)
    
    # Show plot info
    st.info(f"📊 Cell type proportion showing '{selected_cell_type}' grouped by '{selected_sample_group}'")
else:
    st.error(f"Selected grouping variables not found in data")

# Data Information
st.subheader("📋 Data Info")

col1, col2, col3 = st.columns(3)

with col1:
    st.subheader("📊 Summary")
    st.write(f"**Samples**: {sample_display_name}")
    st.write(f"**Cells**: {adata.n_obs:,}")
    st.write(f"**Genes**: {adata.n_vars:,}")
    
    if 'active.ident' in adata.obs.columns:
        cell_types = adata.obs['active.ident'].value_counts()
        st.write("**Cell Types**:")
        for ct, count in cell_types.head(5).items():
            st.write(f"  - {ct}: {count}")
        if len(cell_types) > 5:
            st.write(f"  - ... and {len(cell_types) - 5} more")

with col2:
    # Quality metrics
    st.subheader("📊 Quality Metrics")
    if 'total_counts' in adata.obs.columns:
        st.write(f"**Mean Counts/Cell**: {np.mean(adata.obs['total_counts']):.1f}")
        st.write(f"**Median Counts/Cell**: {np.median(adata.obs['total_counts']):.1f}")
    else:
        st.write("**Mean Counts/Cell**: N/A")
        st.write("**Median Counts/Cell**: N/A")
    
    if 'n_genes_by_counts' in adata.obs.columns:
        st.write(f"**Mean Genes/Cell**: {np.mean(adata.obs['n_genes_by_counts']):.1f}")
        st.write(f"**Median Genes/Cell**: {np.median(adata.obs['n_genes_by_counts']):.1f}")
    else:
        st.write("**Mean Genes/Cell**: N/A")
        st.write("**Median Genes/Cell**: N/A")

with col3:
    # Top genes
    st.subheader("🔝 Top Expressed Genes")
    if 'total_counts' in adata.var.columns:
        top_10_genes = adata.var['total_counts'].nlargest(10)
        for i, (gene, counts) in enumerate(top_10_genes.items(), 1):
            st.write(f"{i}. {gene}: {counts:,.0f}")
    else:
        # Show first 10 gene names if total_counts not available
        for i, gene in enumerate(adata.var_names[:10], 1):
            st.write(f"{i}. {gene}")
