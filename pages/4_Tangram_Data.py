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

registered_catalog = data_manager.get_registered_catalog()

stage_groups = {}
for entry in registered_catalog.values():
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

def _entry_sample_label(entry):
    """Return a concise sample label for UI display."""
    return (
        entry.get("primary_sample")
        or entry.get("key")
        or entry.get("file_stem")
        or entry.get("file_name")
    )

selected_entries = sorted(
    stage_groups[selected_stage],
    key=lambda entry: _entry_sample_label(entry).lower()
)

reload_stage = st.sidebar.button("🔄 Reload Stage Data")

stage_datasets = []
failed_entries = []


for entry in selected_entries:
    dataset_key = entry["key"]
    display_label = _entry_sample_label(entry)
    needs_reload = reload_stage or dataset_key not in data_manager.registered_data
    if needs_reload:
        with st.spinner(f"Loading Tangram data for {display_label}..."):
            data_manager.load_registered_data(dataset_key)
    adata_obj = data_manager.registered_data.get(dataset_key)
    if adata_obj is None:
        failed_entries.append(display_label)
        continue
    stage_datasets.append((entry, adata_obj))

if not stage_datasets:
    st.error(f"Failed to load Tangram data for {selected_stage}")
    st.info("Please check that Tangram data files exist for this stage.")
    st.stop()


def _format_entry_label(entry):
    return _entry_sample_label(entry)


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
        if candidate_key in data_manager.finalized_data:
            return data_manager.finalized_data[candidate_key]
        loaded_spatial = data_manager.load_registered_data(candidate_key, finalize=True)
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
