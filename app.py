#!/usr/bin/env python3
"""Main entry point for the Spatial Transcriptomics Dashboard."""

import streamlit as st
from streamlit import session_state as ss

from utils.home_page import render_home_page


def _initialize_sidebar_state() -> None:
    """Ensure the sidebar starts expanded for first-time visitors."""
    if "sidebar_state" not in ss:
        ss.sidebar_state = "expanded"


def main() -> None:
    """Configure Streamlit and render the combined home page."""
    st.set_page_config(
        page_title="Spatial Transcriptomics Dashboard",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    _initialize_sidebar_state()
    render_home_page()


if __name__ == "__main__":
    main()
