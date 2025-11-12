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

# Friendly labels for dropdown menus
FRIENDLY_NAMES = {
    "active.ident": "Cell Type",
    "panel2_active.ident": "Cell type with mesenchyme subclusters",
    "orig.ident": "Original sample ID",
    "old.ident": "Previous clustering ID",
    "nCount_RNA": "Total RNA counts per cell",
    "nFeature_RNA": "Number of detected genes (RNA)",
    "percent.mt": "% mitochondrial reads",
    "S.Score": "Cell cycle S phase score",
    "G2M.Score": "Cell cycle G2/M phase score",
    "Phase": "Cell cycle phase",
    "CC.Difference": "Cell cycle score difference",
    "nCount_SCT": "Total counts (SCT normalized)",
    "nFeature_SCT": "Number of detected genes (SCT normalized)",
    "total_counts": "Total counts per cell",
    "n_genes_by_counts": "Number of detected genes"
}


def friendly_label(name):
    """Return a user-friendly label for known metadata fields."""
    return FRIENDLY_NAMES.get(name, str(name))


PRIMARY_CELLTYPE_COLUMN = "panel2_active.ident"
FALLBACK_CELLTYPE_COLUMN = "active.ident"


def resolve_celltype_column(columns):
    """Choose the preferred cell type column with a sensible fallback."""
    columns_list = list(columns)
    if PRIMARY_CELLTYPE_COLUMN in columns_list:
        return PRIMARY_CELLTYPE_COLUMN
    if FALLBACK_CELLTYPE_COLUMN in columns_list:
        return FALLBACK_CELLTYPE_COLUMN
    return columns_list[0] if columns_list else None


def order_groupby_columns(columns, primary):
    """Ensure dropdown options prioritise the preferred cell type columns."""
    columns_list = list(columns)
    ordered = []
    if primary and primary in columns_list:
        ordered.append(primary)
    if FALLBACK_CELLTYPE_COLUMN in columns_list and FALLBACK_CELLTYPE_COLUMN not in ordered:
        ordered.append(FALLBACK_CELLTYPE_COLUMN)
    ordered.extend([col for col in columns_list if col not in ordered])
    return ordered


# Initialize managers
@st.cache_resource
def get_managers():
    data_manager = DataManager()
    viz_manager = VisualizationManager()
    return data_manager, viz_manager

data_manager, viz_manager = get_managers()

# Sidebar
st.sidebar.title("🧬 ScRNA-seq Analysis")

st.sidebar.header("🔬 Available Visualizations")
st.sidebar.markdown("""
- **UMAP Plot**: Explore clusters
- **Violin Plot**: Compare feature distributions
- **Dot Plot**: Summarise gene expression
- **Feature Plot**: Visualize genes on embeddings
- **Cell Type Proportion**: View group abundances
""")

# Sample selection (removed UI; use all available samples)
available_samples = data_manager.get_scrna_sample_options()

if not available_samples:
    st.sidebar.error("No sample data found in data directory")
    st.info("Please ensure your data files contain sample identifiers")
    st.stop()

# Data loading
adata = data_manager.load_all_scrna_data()
if adata is None:
    st.error(
        "Unable to load the scRNA-seq reference dataset. "
        "Set `DASHBOARD_DATA_DIR` (or copy `Cartana_simplified_fixname.h5ad` "
        "into `test_data/scRNA-seq`) before launching the dashboard."
    )
    st.stop()

sample_count = len(available_samples)
total_cells = adata.n_obs
obs_columns = adata.obs.columns
primary_celltype = resolve_celltype_column(obs_columns)

# Manual reload button
#if st.sidebar.button("🔄 Reload Data"):
#    data_manager.refresh_scrna_data()
#    with st.spinner("Reloading scRNA-seq data..."):
#        data_manager.load_all_scrna_data()
#    st.rerun()

st.title(f"🧬 ScRNA-seq Data")

# Background
with st.expander("Data Origin and scRNA-seq Processing", expanded=False):
    st.markdown("""
    **Data Origin**  
    - **Newly generated Smart-seq2 data:** GFP⁺ cells were FACS-isolated from E12.5 Nkx2.5-Cre;R26mTmG, Nkx3.2-Cre;R26mTmG and Tg(Prox1-GFP) mouse embryonic pancreata.  
    - **Public embryonic datasets:** E12.5, E14.5 and E17.5 (GSE101099) from the Gene Expression Omnibus (GEO).  
    - **Adult pancreatic mesenchyme datasets:** GEO (GSE125588, GSE176063) and ArrayExpress (E-MTAB-8483).  
    - **Human data:** OMIX database (OMIX001616).  

    **Processing of Smart-seq2 Data**  
    - Sequencing libraries prepared with Smart-seq2 and sequenced on Illumina NextSeq 500.  
    - Reads aligned to mouse genome (GRCm38) using **Kallisto pseudo-alignment** (v0.46.1).  
    - Transcript counts converted to gene counts with **Tximport** (v1.18.0).  
    - **Quality control filters:**  
      - Remove genes expressed in <3 cells.  
      - Exclude cells with <1,000 genes or <100,000 counts.  
      - Exclude cells with >30% mitochondrial counts.  
    - **Normalization and regression:** Log normalization and regression of cell-cycle differences (ScaleData).  
    - **Dimensionality reduction:** PCA, JackStraw to select significant PCs, then **UMAP** visualization.  
    - **Clustering:** Shared nearest neighbor (SNN) graph; FindNeighbors/FindClusters at resolution 1.9.  
    - **Annotation:** Clusters manually curated using marker genes.  

    **Integration of Public Datasets**  
    - Each batch processed independently with **SCTransform** (3,000 features; regress cell-cycle).  
    - Integrated using **Seurat integration workflow**.  
    - Embryonic datasets clustered with first 20 PCs (resolution 1); subclusters at resolution 0.5 for M-VSM analysis.  
    - Adult datasets clustered with 30 PCs (resolution 0.7); integrated with embryonic mesenchyme using **Harmony** (50 PCs, resolution 0.5).  
    - Human-mouse integration performed with **Mousipy** orthologue mapping followed by **MultiMAP** (strength 0.7 and 0.3).  

    **Final Output**  
    - A harmonized, high-coverage **single-cell atlas of embryonic, adult, and human pancreatic cells**.  
    - Provides the reference dataset for spatial mapping and downstream analyses (Tangram, Squidpy, CellChat, Matricom).
    """)

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

# UMAP Plot
st.subheader("📊 UMAP Plot")

# Get available grouping variables
groupby_options = order_groupby_columns(obs_columns, primary_celltype)
selected_groupby = st.selectbox(
    "Group by:",
    groupby_options,
    key="umap_groupby",
    format_func=friendly_label
)

umap_fig = viz_manager.plot_umap(adata, groupby=selected_groupby)
if umap_fig:
    st.plotly_chart(umap_fig, use_container_width=True)

# Violin Plot
st.subheader("🎻 Violin Plot")

violin_controls_col, violin_plot_col = st.columns([1, 2])

# Restrict options to numeric features and categorical groupings
obs_df = adata.obs
numeric_feature_options = obs_df.select_dtypes(include=[np.number]).columns.tolist()
categorical_groupby_options = obs_df.select_dtypes(include=["object", "category", "bool"]).columns.tolist()
categorical_groupby_options = order_groupby_columns(categorical_groupby_options, primary_celltype)

if 'nCount_SCT' in obs_df.columns and 'nCount_SCT' not in numeric_feature_options:
    numeric_feature_options.insert(0, 'nCount_SCT')

selected_feature = None
selected_groupby = None

with violin_controls_col:
    if not numeric_feature_options:
        st.warning("No numeric features available for the violin plot.")
    else:
        selected_feature = st.selectbox(
            "Feature:",
            numeric_feature_options,
            key="violin_feature",
            format_func=friendly_label
        )

    if not categorical_groupby_options:
        st.warning("No categorical annotations available for violin grouping.")
    else:
        selected_groupby = st.selectbox(
            "Group by:",
            categorical_groupby_options,
            key="violin_groupby",
            format_func=friendly_label
        )

with violin_plot_col:
    if selected_feature and selected_groupby:
        viz_manager.plot_violinplot(
            adata,
            groupby=selected_groupby,
            feature=selected_feature,
        )

# Dot Plot
st.subheader("🔴 Dot Plot")

# Get available grouping variables
groupby_options = order_groupby_columns(obs_columns, primary_celltype)
selected_groupby = st.selectbox(
    "Group by:",
    groupby_options,
    key="dot_groupby",
    format_func=friendly_label
)

# Get available genes (all genes in the dataset)
all_genes = adata.var_names.tolist()

# Show top genes by expression for defaults; still expose full list for selection
if 'total_counts' in adata.var.columns:
    display_genes = adata.var['total_counts'].nlargest(100).index.tolist()
else:
    display_genes = all_genes[:100]

default_genes = [gene for gene in display_genes[:5] if gene in all_genes]
if not default_genes:
    default_genes = all_genes[:5]

# Use multiselect for gene selection
gene_list = st.multiselect(
    "Select genes to include in dot plot:",
    options=all_genes,
    default=default_genes,
    help="Select the genes you want to include in the dot plot"
)

# Show selected genes count
if gene_list:
    st.success(f"✅ Selected {len(gene_list)} genes for dot plot")
else:
    st.warning("⚠️ No genes selected for dot plot")

# Create Plotly dot plot (VisualizationManager)
if gene_list and selected_groupby in adata.obs.columns:
    plotly_fig = viz_manager.plot_dotplot(
        adata,
        groupby=selected_groupby,
        genes=gene_list
    )
    if plotly_fig:
        st.plotly_chart(plotly_fig, use_container_width=True)
    #st.info(f"Dot plot shows {len(gene_list)} genes grouped by '{selected_groupby}'")
elif not gene_list:
    st.info("Please select genes to create the dot plot")
else:
    st.error(f"Grouping variable '{selected_groupby}' not found in data")

# Feature Plot
# Gene selection - use top genes for defaults while exposing full list
st.subheader("🎨 Feature Plot")

if 'total_counts' in adata.var.columns:
    top_genes = adata.var['total_counts'].nlargest(20).index.tolist()
else:
    # Fallback: use first 20 genes by name
    top_genes = adata.var_names[:20].tolist()

default_feature_genes = [gene for gene in top_genes[:3] if gene in all_genes]
if not default_feature_genes:
    default_feature_genes = all_genes[:3]

selected_genes = st.multiselect(
    "Select genes:",
    options=all_genes,
    default=default_feature_genes,
    key="feature_genes"
)



if selected_genes:
    scanpy_fig = viz_manager.plot_feature_plot_scanpy(adata, genes=selected_genes)
    if scanpy_fig:
        st.pyplot(scanpy_fig)
        plt.close(scanpy_fig)
else:
    st.info("Select at least one gene to render scanpy feature plots.")

# Cell Type Proportion Plot
st.subheader("🧫 Cell Type Proportion")

# Use fixed grouping variables
selected_cell_type = primary_celltype  # Fixed cell type grouping
selected_sample_group = 'orig.ident'  # This represents the samples users selected

# Create the stacked bar chart
if selected_cell_type in adata.obs.columns and selected_sample_group in adata.obs.columns:
    
    # Create the crosstab and stacked bar plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Compute proportions per sample and adjust axis labels for display only
    plot_df = pd.crosstab(
        adata.obs[selected_cell_type],
        adata.obs[selected_sample_group],
        normalize='columns'
    ).T

    tmp = plot_df.plot(kind='bar', stacked=True, ax=ax)
    
    # Adjust the legend position
    tmp.legend(title= "Cell Type", bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Set labels and title
    #plt.xlabel(selected_sample_group)
    plt.ylabel('Proportion')
    #plt.title(f'Cell Type Proportion')
    
    # Adjust layout to prevent legend overlap
    plt.tight_layout()
    
    # Display the plot
    st.pyplot(fig)
    plt.close(fig)
    
    # Show plot info
    #st.info(f"📊 Cell type proportion showing '{selected_cell_type}' grouped by '{selected_sample_group}'")
else:
    st.error(f"Selected grouping variables not found in data")
