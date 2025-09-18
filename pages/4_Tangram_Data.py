#!/usr/bin/env python3
"""
Tangram Data Analysis Page grouped by developmental stage.
"""

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import tangram as tg
from pathlib import Path
import sys
import re

# Add the project root to the Python path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

# Import utilities
from utils.data_manager import DataManager

st.set_page_config(
    page_title="Tangram Data - Spatial Transcriptomics Dashboard",
    page_icon="🔬",
    layout="wide"
)

# Initialize data manager
@st.cache_resource
def get_data_manager():
    return DataManager()


data_manager = get_data_manager()

# Sidebar
st.sidebar.title("🔬 Tangram Data Analysis")
st.sidebar.header("Tangram Dataset Selection")

tangram_catalog = data_manager.get_tangram_catalog()

stage_groups = {}
for entry in tangram_catalog.values():
    dataset_key = entry.get("key")
    if not dataset_key:
        continue
    primary_value = (entry.get("primary_sample") or dataset_key).upper()
    stage_match = re.match(r"(E\d+)", primary_value)
    stage_label = stage_match.group(1) if stage_match else primary_value
    stage_groups.setdefault(stage_label, []).append(entry)

if not stage_groups:
    st.sidebar.error("No Tangram datasets found in the data directory")
    st.info("Add .h5ad files whose names end with Tangram.h5ad to the data directory.")
    st.stop()

stage_options = sorted(stage_groups.keys())
selected_stage = st.sidebar.selectbox(
    "Select developmental day:",
    stage_options,
    help="Tangram replicates load automatically for the selected developmental stage."
)

selected_entries = sorted(
    stage_groups[selected_stage],
    key=lambda entry: (entry.get("display_name") or entry.get("file_name") or entry["key"]).lower()
)

reload_stage = st.sidebar.button("🔄 Reload Stage Data")

stage_datasets = []
failed_entries = []

for entry in selected_entries:
    dataset_key = entry["key"]
    display_name = entry.get("display_name") or dataset_key
    needs_reload = reload_stage or dataset_key not in data_manager.tangram_data
    if needs_reload:
        with st.spinner(f"Loading Tangram data for {display_name}..."):
            data_manager.load_tangram_data(dataset_key)
    adata_obj = data_manager.tangram_data.get(dataset_key)
    if adata_obj is None:
        failed_entries.append(display_name)
        continue
    stage_datasets.append((entry, adata_obj))

if not stage_datasets:
    st.error(f"Failed to load Tangram data for {selected_stage}")
    st.info("Please check that Tangram data files exist for this stage.")
    st.stop()


def _format_entry_label(entry):
    display_name = entry.get("display_name") or entry["key"]
    file_name = entry.get("file_name")
    if file_name and file_name not in display_name:
        return f"{display_name} • {file_name}"
    return display_name


labeled_stage_datasets = []
used_labels = set()
for entry, adata_obj in stage_datasets:
    label = _format_entry_label(entry)
    if label in used_labels:
        label = f"{label} [{entry['key']}]"
    used_labels.add(label)
    labeled_stage_datasets.append((label, entry, adata_obj))

st.sidebar.header("📊 Data Status")
for label, entry, adata_obj in labeled_stage_datasets:
    st.sidebar.success(f"{label}\n{adata_obj.n_obs} cells")
for display_name in failed_entries:
    st.sidebar.error(f"{display_name}\nFailed to load")

# Aggregate metrics across replicates
stage_total_cells = 0
stage_total_counts = 0.0
gene_sets = []

for _, adata_obj in stage_datasets:
    gene_sets.append({str(gene) for gene in adata_obj.var_names})
    if "total_counts" in adata_obj.obs:
        counts = np.asarray(adata_obj.obs["total_counts"], dtype=float)
    else:
        counts = np.asarray(np.sum(adata_obj.X, axis=1)).ravel()
    stage_total_cells += int(adata_obj.n_obs)
    stage_total_counts += float(np.sum(counts))

if gene_sets:
    if len(gene_sets) == 1:
        common_genes = sorted(next(iter(gene_sets)), key=str.lower)
    else:
        common_genes = sorted(set.intersection(*gene_sets), key=str.lower)
else:
    common_genes = []

if common_genes:
    available_genes = common_genes
else:
    union_genes = set()
    for gene_set in gene_sets:
        union_genes.update(gene_set)
    available_genes = sorted(union_genes, key=str.lower)

avg_counts = stage_total_counts / stage_total_cells if stage_total_cells else 0.0
shared_gene_count = len(common_genes) if common_genes else len(available_genes)

st.title(f"🔬 Tangram Data - {selected_stage}")

col1, col2, col3, col4 = st.columns(4)
with col1:
    st.metric("Replicates", f"{len(stage_datasets)}")
with col2:
    st.metric("Cells", f"{stage_total_cells:,}")
with col3:
    st.metric("Shared Genes", f"{shared_gene_count:,}")
with col4:
    st.metric("Avg Counts/Cell", f"{avg_counts:.1f}")

st.caption(f"Total counts across replicates: {stage_total_counts:,.0f}")


def _sanitize_gene_name(name: str) -> str:
    return "".join(ch for ch in name.lower() if ch.isalnum())


def _match_default_genes(defaults, candidates):
    lookup = {_sanitize_gene_name(gene): gene for gene in candidates}
    matched = []
    for gene in defaults:
        key = _sanitize_gene_name(gene)
        value = lookup.get(key)
        if value and value not in matched:
            matched.append(value)
    return matched


default_gene_candidates = ["Ptf1a", "Gcg", "Col6a1", "Hoxb6", "Mki67"]
default_genes = _match_default_genes(default_gene_candidates, available_genes)
if not default_genes and available_genes:
    default_genes = available_genes[: min(5, len(available_genes))]

st.subheader("Gene Selection")
selected_genes = st.multiselect(
    "Genes to visualise",
    options=available_genes,
    default=default_genes,
    help="Select one or more genes shared across the loaded Tangram replicates."
)

plot_kwargs = dict(spot_size=50, scale_factor=0.1, perc=0.01, cmap="inferno")


def _resolve_measured_dataset(entry):
    candidate_ids = [
        entry.get("primary_sample"),
        entry.get("key"),
        entry.get("file_stem"),
        entry.get("file_name"),
    ]
    for candidate in candidate_ids:
        candidate_key = (candidate or "").strip()
        if not candidate_key:
            continue
        if candidate_key in data_manager.spatial_data:
            return data_manager.spatial_data[candidate_key]
        loaded_spatial = data_manager.load_spatial_data(candidate_key)
        if loaded_spatial is not None:
            return loaded_spatial
    return None


if selected_genes:
    st.subheader("Gene Expression Across Replicates")
    tab_labels = [label for label, _, _ in labeled_stage_datasets]
    if len(tab_labels) > 1:
        tab_contexts = list(zip(st.tabs(tab_labels), labeled_stage_datasets))
    else:
        tab_contexts = [(st.container(), labeled_stage_datasets[0])]
    measured_cache = {}
    for tab, (label, entry, adata_obj) in tab_contexts:
        with tab:
            adata_measured = measured_cache.get(label)
            if adata_measured is None:
                adata_measured = _resolve_measured_dataset(entry) or adata_obj
                measured_cache[label] = adata_measured
            try:
                fig = tg.plot_genes_sc(
                    selected_genes,
                    adata_measured=adata_measured,
                    adata_predicted=adata_obj,
                    return_figure=True,
                    **plot_kwargs,
                )
                st.pyplot(fig)
                plt.close(fig)
            except Exception as exc:
                st.warning(f"Unable to render gene panel for {label}: {exc}")
else:
    st.info("Select at least one gene to display expression plots.")

st.subheader("Tangram Data Preview")
preview_labels = [label for label, _, _ in labeled_stage_datasets]
preview_map = {label: (entry, adata) for label, entry, adata in labeled_stage_datasets}
selected_preview = st.selectbox(
    "Dataset to preview:",
    preview_labels,
    help="Inspect metadata and file information for a specific Tangram replicate.",
)

preview_entry, preview_adata = preview_map[selected_preview]

col1, col2 = st.columns(2)

with col1:
    st.write("**Sample of Tangram data:**")
    st.dataframe(preview_adata.obs.head(10))
    st.write("**Available metadata columns:**")
    st.write(list(preview_adata.obs.columns))

with col2:
    st.write("**Data Statistics:**")
    st.write(f"- **Shape**: {preview_adata.shape}")
    if preview_adata.X.size:
        sparsity = (1 - np.count_nonzero(preview_adata.X) / preview_adata.X.size) * 100
    else:
        sparsity = 0.0
    st.write(f"- **Sparsity**: {sparsity:.1f}%")
    st.write(f"- **Memory usage**: {preview_adata.nbytes / 1024 / 1024:.1f} MB")
    st.write("**Quality Metrics:**")
    if {"total_counts", "n_genes_by_counts"}.issubset(preview_adata.obs.columns):
        st.write(f"- **Mean counts/cell**: {np.mean(preview_adata.obs['total_counts']):.1f}")
        st.write(f"- **Median counts/cell**: {np.median(preview_adata.obs['total_counts']):.1f}")
        st.write(f"- **Mean genes/cell**: {np.mean(preview_adata.obs['n_genes_by_counts']):.1f}")
        st.write(f"- **Median genes/cell**: {np.median(preview_adata.obs['n_genes_by_counts']):.1f}")
    else:
        st.write("- Total count metrics are not available for this dataset.")

st.subheader("File Information")
data_dir = Path(data_manager.data_dir)
preview_key = preview_entry.get("primary_sample") or preview_entry.get("key")
file_stem = preview_entry.get("file_stem") or preview_entry.get("key")
target_label = preview_key or file_stem or "selected dataset"

if preview_key:
    dapi_files = list(data_dir.glob(f"*{preview_key}*DAPI*.tif"))
else:
    dapi_files = []

if dapi_files:
    st.success(f"✅ Found {len(dapi_files)} DAPI image(s) for {target_label}:")
    for dapi_file in dapi_files:
        st.write(f"  - {dapi_file.name}")
else:
    st.warning(f"⚠️ No DAPI images found for {target_label}")
    st.info("Ensure DAPI images follow the naming convention: `{sample}_DAPI.tif`")

tangram_patterns = []
if file_stem:
    tangram_patterns.append(f"*{file_stem}*tangram*.h5ad")
    tangram_patterns.append(f"*{file_stem}*Tangram*.h5ad")
elif preview_key:
    tangram_patterns.append(f"*{preview_key}*tangram*.h5ad")
    tangram_patterns.append(f"*{preview_key}*Tangram*.h5ad")

tangram_files = []
for pattern in tangram_patterns:
    tangram_files.extend(data_dir.glob(pattern))

if tangram_files:
    st.success(f"✅ Found {len(tangram_files)} Tangram data file(s) for {target_label}:")
    for tangram_file in tangram_files:
        st.write(f"  - {tangram_file.name}")
else:
    st.warning(f"⚠️ No Tangram data files found for {target_label}")
