#!/usr/bin/env python3
"""Tangram gene expression explorer for early developmental samples."""

from __future__ import annotations

import re
from pathlib import Path
import sys
from dataclasses import dataclass
from typing import List, Optional

import anndata
import matplotlib.pyplot as plt
import streamlit as st
import tangram as tg

# Ensure project utilities are importable when page is executed directly.
PROJECT_ROOT = Path(__file__).parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

from utils.data_manager import DataManager, format_sample_label  # noqa: E402

st.set_page_config(
    page_title="Tangram Data - Spatial Transcriptomics Dashboard",
    page_icon="🧬",
    layout="wide",
)

TARGET_STAGE_PREFIXES = ("E14.5", "E12.5")


@st.cache_resource
def get_data_manager() -> DataManager:
    return DataManager()


data_manager = get_data_manager()
tangram_root = data_manager.spatial_data_dir


@st.cache_resource
def _load_predicted_dataset(path_str: str) -> anndata.AnnData:
    """Read and memoise a Tangram predicted dataset."""
    return anndata.read_h5ad(path_str)


@st.cache_resource
def _load_measured_dataset(path_str: str) -> anndata.AnnData:
    """Read and memoise the measured Tangram dataset."""
    return anndata.read_h5ad(path_str)


@dataclass(frozen=True)
class SampleOption:
    """UI representation of a Tangram dataset and its prediction pair."""

    key: str
    label: str
    measured_path: Path
    trained_path: Path


def _matches_stage_prefix(value: str) -> bool:
    """Return True when a value begins with one of the target stage prefixes."""
    if not value:
        return False
    upper_value = value.upper()
    return any(upper_value.startswith(prefix.upper()) for prefix in TARGET_STAGE_PREFIXES)


def _infer_label_from_name(value: str) -> str:
    formatted = format_sample_label(value)
    return formatted or value


def _find_measured_partner(trained_path: Path) -> Optional[Path]:
    """Locate the measured Tangram dataset corresponding to a trained file."""
    if not trained_path.exists():
        return None

    folder = trained_path.parent
    stem = trained_path.stem
    base_name = re.sub(r"(?i)_trained$", "", stem)

    candidate_names = (
        f"{base_name}_Tangram.h5ad",
        f"{base_name}_tangram.h5ad",
        stem.replace("trained", "Tangram") + ".h5ad",
        stem.replace("trained", "tangram") + ".h5ad",
        base_name + "_Tangram.h5ad",
    )

    for name in candidate_names:
        candidate = folder / name
        if candidate.is_file() and "tangram" in candidate.name.lower():
            return candidate

    base_fragment = base_name.lower()
    for file_path in sorted(folder.glob("*Tangram*.h5ad")):
        if not file_path.is_file():
            continue
        if base_fragment and base_fragment not in file_path.stem.lower():
            continue
        return file_path

    return None


def _trained_files() -> List[Path]:
    """Return *_trained.h5ad files under the tangram root directory."""
    if not tangram_root.exists():
        return []

    trained_files = []
    for file_path in sorted(tangram_root.rglob("*.h5ad")):
        name_lower = file_path.name.lower()
        if not name_lower.endswith("trained.h5ad"):
            continue
        if any(part.lower() == "scrna-seq" for part in file_path.parts):
            continue
        trained_files.append(file_path)

    return trained_files


def _build_sample_options() -> List[SampleOption]:
    """Return sidebar options constrained to the target developmental stage."""
    options: List[SampleOption] = []
    for trained_path in _trained_files():
        base_name = re.sub(r"(?i)_trained$", "", trained_path.stem)
        if not _matches_stage_prefix(base_name):
            continue

        measured_path = _find_measured_partner(trained_path)
        if measured_path is None:
            continue

        option_key = base_name or trained_path.stem
        option_label = _infer_label_from_name(base_name or trained_path.stem)

        options.append(
            SampleOption(
                key=option_key,
                label=option_label,
                measured_path=measured_path,
                trained_path=trained_path,
            )
        )

    options.sort(key=lambda item: item.label.lower())
    return options


# ---- Sidebar --------------------------------------------------------------

st.sidebar.title("Tangram Data")
sample_options = _build_sample_options()

if not sample_options:
    stages_display = ", ".join(TARGET_STAGE_PREFIXES)
    st.sidebar.error(f"No Tangram datasets found for stages: {stages_display}.")
    st.stop()

st.sidebar.header("🔬 Available Visualizations")
st.sidebar.markdown("""
- **Spatial gene expression**: Compare predicted vs measured signal
""")

st.title("Tangram Data Explorer")

sample_labels = [item.label for item in sample_options]
selection_col, gene_col = st.columns([1, 1])

with selection_col:
    selected_label = st.selectbox("Select sample", sample_labels)

selected_sample = next(item for item in sample_options if item.label == selected_label)

# ---- Load Tangram data ----------------------------------------------------

with st.spinner(f"Loading Tangram measurements for {selected_label}..."):
    adata_measured = _load_measured_dataset(str(selected_sample.measured_path))

if adata_measured is None:
    st.error(f"Unable to load Tangram data for {selected_label}.")
    st.stop()

with st.spinner(f"Loading Tangram predictions for {selected_label}..."):
    adata_predicted = _load_predicted_dataset(str(selected_sample.trained_path))

# ---- Gene selection -------------------------------------------------------

measured_genes = {str(gene) for gene in adata_measured.var_names}
predicted_genes = [str(gene) for gene in adata_predicted.var_names]
available_genes = sorted(predicted_genes, key=str.lower)

if not available_genes:
    st.error("No genes were found in the trained Tangram dataset.")
    st.stop()

preferred_defaults = ["gcg", "ptf1a", "nkx6.1", "Rbpjl"]
available_lookup = {gene.lower(): gene for gene in available_genes}
default_genes = [
    available_lookup[name.lower()]
    for name in preferred_defaults
    if name.lower() in available_lookup
]

if not default_genes:
    default_genes = available_genes[:4]
with gene_col:
    selected_genes = st.multiselect(
        "Select genes",
        available_genes,
        default=default_genes,
        help="Pick genes to visualise across measured and predicted spaces.",
    )

if not selected_genes:
    st.info("Select at least one gene to generate a plot.")
    st.stop()

missing_genes = sorted({gene for gene in selected_genes if gene not in measured_genes}, key=str.lower)
if missing_genes:
    missing_display = ", ".join(missing_genes)
    st.warning(
        "Measured Tangram data do not contain: "
        f"{missing_display}. Predicted expression will still be shown where possible."
    )

# ---- Plotting -------------------------------------------------------------

st.subheader(f"Spatial gene expression for {selected_label}")

try:
    plt.close("all")
    tg.plot_genes_sc(
        selected_genes,
        adata_measured=adata_measured,
        adata_predicted=adata_predicted,
        spot_size=50,
        scale_factor=0.1,
        perc=0.01,
        return_figure=False,
        cmap="inferno",
    )
    fig = plt.gcf()
    plot_col, _ = st.columns([0.7, 0.3])
    with plot_col:
        st.pyplot(fig, use_container_width=True)
finally:
    plt.close("all")

# Provide context statistics for users exploring genes.
#st.markdown(
#    f"**Loaded cells:** {adata_measured.n_obs:,} · "
#    f"**Genes available:** {len(available_genes):,}"
#)
