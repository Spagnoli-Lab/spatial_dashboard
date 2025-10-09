#!/usr/bin/env python3
"""Utilities for rendering the shared home page content."""

import streamlit as st

def render_home_page() -> None:
    
    # --- Dashboard Home Page ---
    st.title("🧬 Pancreatic Spatial & Single-Cell Transcriptomics Atlas")

    st.markdown("""
    ## **About this Study**
    This dashboard presents the first high-resolution, single-cell and spatial transcriptomic map of the developing pancreas.  
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
    - **Cell–cell and ECM signalling:** CellChat and Matricom identify ligand–receptor and ECM interactions between mesenchymal niches and pancreatic progenitors.  

    Together these methods create a reproducible workflow that moves from raw sequencing reads to an interactive, spatially resolved atlas.
    
    ---
    """)

    # Available pages
    #st.markdown(
    """
    ### Available Pages
    - **ScRNA-seq Data**: Single-cell RNA sequencing analysis  
    - **Spatial Data**: Spatial transcriptomics exploration  
    - **Tangram Data**: Tangram mapping with Tissuumaps integration
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
    - **Spatial Scatter Plots**: Cell distribution in tissue space  
    - **Neighbourhood Enrichment**: Spatial interaction analysis  
    - **UMAP Integration**: Combined spatial and transcriptomic views  
    - **Cell Type Mapping**: Spatial distribution of cell types
    """
    )

    st.markdown(
    """
    ### Tangram Results
    - **Tissuumaps Integration**: Interactive spatial visualisation  
    - **DAPI Image Overlay**: Tangram mapping with tissue images  
    - **Sample-specific Analysis**: E12, E14, E17 developmental stages  
    - **Real-time Exploration**: Dynamic data exploration

    ---
    """
    )

    st.subheader("Available Samples")
    st.markdown(
    """
    | Stage | Description | Key Features |
    |-------|-------------|--------------|
    | **E12.5** | Early embryonic development | Initial cell type specification |
    | **E14.5** | Mid-embryonic development | Tissue patterning and growth |
    ---
    """
    )

    # Technical Stack
    #st.markdown(
    
    
    ### Data Formats
    #- **H5AD**: AnnData objects for single-cell data  
    #- **TIF**: DAPI staining for spatial reference  
    #- **CSV**: Metadata and expression matrices  
    """
    ## Technical Stack
    ### Analysis Tools
    - **Scanpy**: Single-cell analysis pipeline  
    - **Tissuumaps**: Spatial data visualisation  
    - **Tangram**: Spatial mapping algorithms
    
    ### Visualisation Libraries
    - **Plotly**: Interactive plots and dashboards  
    - **Matplotlib**: Static scientific figures  
    - **Seaborn**: Statistical visualisation  

    ### Platform
    - **Streamlit**: Web application framework  
    - **Python**: Scientific computing environment
    
    ---
    """
    #)

    st.header("Getting Started")
    st.markdown(
    """
    ### Quick Start Guide
    1. Navigate to a page using the sidebar  
    2. Select your sample (E12.5, E14.5)  
    3. Load the relevant datasets  
    4. Explore the interactive visualisations  
    5. Compare stages by switching between samples  

    ### Tips for Best Experience
    - Start with scRNA-seq data for baseline insights  
    - Use Tissuumaps to investigate spatial context  
    - Compare developmental stages to spot trends  
    - Export figures and data for downstream analysis
    """
    )


    st.markdown("---")

    # --- GitHub Link Section ---
    st.header("GitHub Repository")
    st.markdown(
    """
    Explore the full source code and documentation on GitHub:  
    [https://github.com/Spagnoli-Lab/spatial-pancreas-communities](https://github.com/Spagnoli-Lab/spatial-pancreas-communities)
    """
    )

    st.markdown("---")

    # --- Citation, Author, Maker at the END ---
    st.markdown(
    """
    ## Citation
    ### **Spatially organized cellular communities shape functional tissue architecture in the pancreas**  
    Alejo Torres-Cano¹, Jean-Francois Darrigrand¹, Gabriel Herrera-Oropeza², Georgina Goss¹, David Willnow¹³, Anna Salowka¹, Siwanart Ma¹, Debashish Chitnis⁴, Morgane Rouault⁴, Alessandra Vigilante¹⁵, Francesca M. Spagnoli¹*  

    ¹Centre for Gene Therapy and Regenerative Medicine, King’s College London; Great Maze Pond, London SE1 9RT, UK.  
    ²Centre for Neurodevelopmental Biology, Institute of Psychiatry, Psychology & Neuroscience, King’s College London; London SE1 1UL, UK.  
    ³Present address: The Francis Crick Institute; London, NW1 1AT, UK.  
    ⁴10X Genomics, Pleasanton, CA, USA.  
    ⁵Hub for Applied Bioinformatics (HAB) King’s College London; Newcomen St, London SE1 1UL, UK.  

    *Corresponding author:* **Francesca M. Spagnoli** — [francesca.spagnoli@kcl.ac.uk](mailto:francesca.spagnoli@kcl.ac.uk)
    """
    )

    st.subheader("BibTeX Citation")
    st.code("""
    @article{TorresCano2025,
    title={Spatially organized cellular communities shape functional tissue architecture in the pancreas},
    author={Torres-Cano, Alejo and Darrigrand, Jean-Francois and Herrera-Oropeza, Gabriel and Goss, Georgina and Willnow, David and Salowka, Anna and Ma, Siwanart and Chitnis, Debashish and Rouault, Morgane and Vigilante, Alessandra and Spagnoli, Francesca M.},
    journal={},
    year={2025},
    note={*Corresponding author: Francesca M. Spagnoli (francesca.spagnoli@kcl.ac.uk)}
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
        st.markdown("*Powered by Streamlit, Scanpy, Squidpy, Tangram, Plotly, Matplotlib, Seaborn, and Tissuumaps*")

    with right_col:
        st.image(
            "/Users/mayongzhi/Desktop/FS_lab/dashboard/utils/King's_College_London_logo.svg.png",
            use_container_width=True
        )

