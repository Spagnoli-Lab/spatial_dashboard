#!/usr/bin/env python3
"""Utilities for rendering the shared home page content."""

from pathlib import Path

import streamlit as st

LOGO_PATH = Path(__file__).resolve().with_name("King's_College_London_logo.svg.png")


def render_home_page() -> None:
    
    # --- Dashboard Home Page ---
    st.title("🧬 Pancreatic Spatial & Single-Cell Transcriptomics Atlas")

    st.markdown("""
    ## **About this Study**
    This dashboard presents a high-resolution, single-cell and spatial transcriptomic map of the mouse developing pancreas.  
    By integrating scRNA-seq data with spatially resolved transcriptomics, we chart how diverse pancreatic cell types—especially poorly defined mesenchymal populations—organize across space and time during development.  
    This approach reveals the multicellular “niches” that shape endocrine and exocrine tissue architecture and provides a framework to guide in-vitro organogenesis and tissue engineering for pancreatic disease.

    ---

    ### **What We Achieve**
    Using **Tangram**, we project single-cell gene expression profiles onto spatial coordinates at single-cell resolution. This enables us to:

    - Reconstruct the spatial location of thousands of genes beyond the original targeted panel.  
    - Identify epithelial–mesenchymal units and their niche-specific ECM components.  
    - Highlight conserved pro-endocrine environments, such as the M-II mesenchyme enriched with Collagen VI, in both mouse and human pancreas.  

    Overall, we connect **cell identity** to **spatial position**, providing unprecedented insight into how developing tissues self-organize.

    ---

    ### **Bioinformatics and Data Analysis in Brief**
    - **Single-cell RNA-seq integration:** Smart-seq2 plus public embryonic and adult pancreas datasets processed in Seurat (SCTransform, Harmony, MultiMAP).  
    - **Spatial transcriptomics:** High-plex Hybridization-based in situ sequencing (HybISS) analysed by segmentation-free (SSAM) and segmentation-based pipelines.  
    - **Gene expression mapping:** Tangram used to impute full transcriptomes onto spatial maps; Squidpy used for neighbourhood enrichment and tissue-domain analysis.  
    
    Together these methods create a reproducible workflow that moves from raw sequencing reads to an interactive, spatially resolved atlas.
    
    ---
    """)

    # Available pages
    #st.markdown(
    """
    ### Available Pages
    - **ScRNA-seq Data**: Single-cell RNA sequencing analysis  
    - **Spatial Data**: Spatial transcriptomics exploration  
    - **Tangram Data**: Tangram mapping results
    ---
    """
    #)

    # Key Features
    st.markdown(
    """
    ## Key Features

    

    ### ScRNA-seq Analysis
    - **UMAP Visualisation**: Interactive dimensionality reduction plots  
    - **Dot Plots**: Gene expression across cell types  
    - **Feature Plots**: Gene expression on UMAP coordinates  
    - **Cell Type Proportions**: Composition analysis across samples
    """
    )

    st.markdown(
    """
    ### Spatial Data Analysis
    - **DAPI Images**
    - **Spatial Scatter Plots**: Cell distribution in tissue space  
    - **Neighbourhood Enrichment**: Spatial interaction analysis 
    - **Cell Type Mapping**: Spatial distribution of cell types
    """
    )

    st.markdown(
    """
    ### Tangram Mapping Results
    - **Spatial Scatter Plots**: Real spatial gene expression and Tangram-predicted gene expression

    ---
    """
    )

    """
    ## Technical Stack
    ### Analysis Tools
    - **Scanpy**: Single-cell analysis pipeline  
    - **Squidpy**: Spatial data visualisation  
    - **Tangram**: Spatial mapping algorithms
    
    ### Visualisation Libraries
    - **Plotly**  
    - **Matplotlib**  
    - **Seaborn**  

    ### Platform
    - **Streamlit**: Web application framework  
    - **Python**: Scientific computing environment
    
    ---
    """
    #)


    st.markdown("---")

    # --- GitHub Link Section ---
    st.header("GitHub Repository")
    st.markdown(
    """
    Explore the full source code and documentation of this project:  
    [https://github.com/Spagnoli-Lab/spatial-pancreas-communities](https://github.com/Spagnoli-Lab/spatial-pancreas-communities)
    The repository of this dashboard:
    [https://github.com/Spagnoli-Lab/spatial_dashboard](https://github.com/Spagnoli-Lab/spatial_dashboard)
    """
    )

    st.markdown("---")

    # --- Citation, Author, Maker at the END ---
    st.markdown(
    """
    ## Citation
    Alejo Torres-Cano et al. ,Spatially organized cellular communities shape functional tissue architecture in the pancreas.Sci. Adv.11,eadx5791(2025).DOI:10.1126/sciadv.adx5791
    """
    )

    st.subheader("BibTeX Citation")
    st.code("""
    @article{torres2025spatially,
    title={Spatially organized cellular communities shape functional tissue architecture in the pancreas},
    author={Torres-Cano, Alejo and Darrigrand, Jean-Francois and Herrera-Oropeza, Gabriel and Goss, Georgina and Willnow, David and Salowka, Anna and Ma, Siwanart and Chitnis, Debashish and Rouault, Morgane and Vigilante, Alessandra and others},
    journal={Science Advances},
    volume={11},
    number={46},
    pages={eadx5791},
    year={2025},
    publisher={American Association for the Advancement of Science}
    }
    """, language="bibtex")

    # Dashboard Created By
    st.markdown(
    """
    ## Dashboard Created By
    **Siwanart Ma**
    
    ---
    """)

    # --- Add a logo/image at the end on the right ---
    left_col, right_col = st.columns([6, 1])  # left column bigger, right column smaller

    with left_col:
        st.markdown("*Powered by Streamlit, Scanpy, Squidpy, Tangram, Plotly, Matplotlib, and Seaborn*")

    with right_col:
        if LOGO_PATH.exists():
            st.image(str(LOGO_PATH), use_container_width=True)
