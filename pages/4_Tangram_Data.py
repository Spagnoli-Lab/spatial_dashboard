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

# ---- Helper utilities -----------------------------------------------------

def _group_catalog_by_stage(catalog):
    """Return registered entries grouped by developmental stage label."""
    stage_groups = {}
    for entry in catalog.values():
        dataset_key = entry.get("key")
        if not dataset_key:
            continue
        primary_value = (entry.get("primary_sample") or dataset_key).upper()
        stage_match = re.match(r"(E\d+)", primary_value)
        stage_label = stage_match.group(1) if stage_match else primary_value
        stage_groups.setdefault(stage_label, []).append(entry)
    return stage_groups


def _clean_sample_label(label: str) -> str:
    if not label:
        return ""
    return label.split("(", 1)[0].strip()


def _entry_sample_label(entry):
    """Return a concise sample label for UI display."""
    candidates = [
        entry.get("primary_sample"),
        entry.get("key"),
        entry.get("file_stem"),
        entry.get("file_name"),
    ]
    for candidate in candidates:
        cleaned = _clean_sample_label(str(candidate)) if candidate else ""
        if cleaned:
            return cleaned
    fallback = _clean_sample_label(str(entry.get("key", "")))
    return fallback or str(entry.get("key", ""))


def _load_stage_datasets(selected_entries, should_reload):
    """Load Tangram data for the selected replicates."""
    stage_datasets = []
    failed_entries = []

    for entry in selected_entries:
        dataset_key = entry["key"]
        display_label = _entry_sample_label(entry)
        needs_reload = should_reload or dataset_key not in data_manager.registered_data
        if needs_reload:
            with st.spinner(f"Loading Tangram data for {display_label}..."):
                data_manager.load_registered_data(dataset_key)
        adata_obj = data_manager.registered_data.get(dataset_key)
        if adata_obj is None:
            failed_entries.append(display_label)
            continue
        stage_datasets.append((entry, adata_obj))

    return stage_datasets, failed_entries


def _label_stage_datasets(stage_datasets):
    """Attach unique UI labels to each replicate entry."""
    labeled_stage_datasets = []
    used_labels = set()
    for entry, adata_obj in stage_datasets:
        label = _entry_sample_label(entry)
        if label in used_labels:
            label = f"{label} [{entry['key']}]"
        used_labels.add(label)
        labeled_stage_datasets.append((label, entry, adata_obj))
    return labeled_stage_datasets


def _sanitize_gene_name(name: str) -> str:
    """Return a normalised gene identifier for lookups."""
    return "".join(ch for ch in name.lower() if ch.isalnum())


def _match_default_genes(defaults, candidates):
    """Cross-reference preferred defaults against available genes."""
    lookup = {_sanitize_gene_name(gene): gene for gene in candidates}
    matched = []
    for gene in defaults:
        key = _sanitize_gene_name(gene)
        value = lookup.get(key)
        if value and value not in matched:
            matched.append(value)
    return matched


def _resolve_measured_dataset(entry):
    """Best-effort lookup of the measured dataset associated with a replicate."""
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


# ---- Sidebar configuration -------------------------------------------------

st.sidebar.title("🔬 Tangram Data Analysis")
st.sidebar.header("Tangram Dataset Selection")

registered_catalog = data_manager.get_registered_catalog()
stage_groups = _group_catalog_by_stage(registered_catalog)

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
    key=lambda entry: _entry_sample_label(entry).lower()
)

reload_stage = st.sidebar.button("🔄 Reload Stage Data")

stage_datasets, failed_entries = _load_stage_datasets(selected_entries, reload_stage)

if not stage_datasets:
    st.error(f"Failed to load Tangram data for {selected_stage}")
    st.info("Please check that Tangram data files exist for this stage.")
    st.stop()

labeled_stage_datasets = _label_stage_datasets(stage_datasets)

# ---- Stage-level summaries -------------------------------------------------

stage_total_cells = 0
stage_total_counts = 0.0
gene_sets = []

# Track shared gene sets and sequencing depth across the loaded replicates.
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

# Derived metrics retained for potential sidebar or summary displays.
avg_counts = stage_total_counts / stage_total_cells if stage_total_cells else 0.0
shared_gene_count = len(common_genes) if common_genes else len(available_genes)

st.title(f"🔬 Tangram Data - {selected_stage}")

# ---- Gene selection UI ----------------------------------------------------

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


# ---- Gene expression plots -------------------------------------------------

if selected_genes:
    st.subheader("Gene Expression Across Replicates")
    tab_labels = [label for label, _, _ in labeled_stage_datasets]
    if len(tab_labels) > 1:
        tab_contexts = list(zip(st.tabs(tab_labels), labeled_stage_datasets))
    else:
        tab_contexts = [(st.container(), labeled_stage_datasets[0])]
    # Cache measured datasets per tab to keep the UI responsive on reruns.
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
