#!/usr/bin/env python3
"""Tangram gene expression explorer for early developmental samples."""

from __future__ import annotations

import re
from pathlib import Path
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional

import anndata
import matplotlib.pyplot as plt
import streamlit as st
import tangram as tg

# Ensure project utilities are importable when page is executed directly.
PROJECT_ROOT = Path(__file__).parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

from utils.data_manager import DataManager  # noqa: E402

st.set_page_config(
    page_title="Tangram Gene Explorer",
    page_icon="🧬",
    layout="wide",
)

TARGET_STAGE_PREFIXES = ("E14.5", "E12.5")


@st.cache_resource
def get_data_manager() -> DataManager:
    return DataManager()


data_manager = get_data_manager()


@st.cache_resource
def _load_predicted_dataset(path_str: str) -> anndata.AnnData:
    """Read and memoise a Tangram predicted dataset."""
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


def _is_target_stage(entry: Dict[str, object]) -> bool:
    """Return True when the catalog entry belongs to one of the target stages."""
    candidates = [
        entry.get("primary_sample"),
        entry.get("file_stem"),
        entry.get("display_name"),
    ]
    return any(_matches_stage_prefix(str(candidate or "")) for candidate in candidates)


def _shorten_label(value: str) -> str:
    """Return a concise label without auxiliary annotations."""
    label = str(value or "").strip()
    if not label:
        return ""
    if " (" in label:
        label = label.split(" (", 1)[0].strip()
    return label


def _find_trained_path(tangram_path: Path) -> Optional[Path]:
    """Locate the trained predictions that correspond to a Tangram dataset."""
    if not tangram_path.exists():
        return None

    folder = tangram_path.parent
    stem = tangram_path.stem
    base_name = re.sub(r"(?i)_tangram$", "", stem)

    candidate_names = {
        f"{base_name}_trained.h5ad",
        f"{base_name}_Trained.h5ad",
        stem.replace("Tangram", "trained") + ".h5ad",
        stem.replace("tangram", "trained") + ".h5ad",
    }

    for name in candidate_names:
        candidate = folder / name
        if candidate.exists():
            return candidate

    base_fragment = base_name.lower()
    fallback_match: Optional[Path] = None
    for file_path in sorted(folder.glob("*.h5ad")):
        lowered = file_path.name.lower()
        if "trained" not in lowered:
            continue
        if base_fragment and base_fragment in lowered:
            return file_path
        if fallback_match is None:
            fallback_match = file_path

    if fallback_match is not None:
        return fallback_match

    return None


def _build_sample_options() -> List[SampleOption]:
    """Return sidebar options constrained to the target developmental stage."""
    catalog = data_manager.get_registered_catalog() or {}
    options: List[SampleOption] = []

    for entry in catalog.values():
        if not _is_target_stage(entry):
            continue

        entry_path = Path(entry.get("path", ""))
        if not entry_path.exists() or "tangram" not in entry_path.name.lower():
            continue

        trained_path = _find_trained_path(entry_path)
        if trained_path is None or "trained" not in trained_path.name.lower():
            continue

        option_key = str(entry.get("key") or entry_path.stem)
        option_label = _shorten_label(entry.get("display_name")) or _shorten_label(option_key)

        options.append(
            SampleOption(
                key=option_key,
                label=option_label,
                measured_path=entry_path,
                trained_path=trained_path,
            )
        )

    options.sort(key=lambda item: item.label.lower())
    return options


# ---- Sidebar --------------------------------------------------------------

st.sidebar.title("Tangram Gene Explorer")
sample_options = _build_sample_options()

if not sample_options:
    stages_display = ", ".join(TARGET_STAGE_PREFIXES)
    st.sidebar.error(f"No Tangram datasets found for stages: {stages_display}.")
    st.stop()

sample_labels = [item.label for item in sample_options]
selected_label = st.sidebar.selectbox("Select sample", sample_labels)
selected_sample = next(item for item in sample_options if item.label == selected_label)

# ---- Load Tangram data ----------------------------------------------------

with st.spinner(f"Loading Tangram measurements for {selected_label}..."):
    adata_measured = data_manager.load_registered_data(selected_sample.key)

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
selected_genes = st.sidebar.multiselect(
    "Select genes",
    available_genes,
    default=default_genes,
    help="Pick genes to visualise across measured and predicted spaces.",
)

if not selected_genes:
    st.info("Select at least one gene from the sidebar to generate a plot.")
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
    st.pyplot(fig, use_container_width=True)
finally:
    plt.close("all")

# Provide context statistics for users exploring genes.
st.markdown(
    f"**Loaded cells:** {adata_measured.n_obs:,} · "
    f"**Genes available:** {len(available_genes):,}"
)

#st.caption(
#    f"Measured: {selected_sample.measured_path.name} · "
#    f"Predicted: {selected_sample.trained_path.name}"
#)
