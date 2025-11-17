
# Spatial pancreatic atlas dashboard

In our paper [Spatially organized cellular communities shape functional tissue architecture in the pancreas](https://eur03.safelinks.protection.outlook.com/?url=https%3A%2F%2Fdoi.org%2F10.1126%2Fsciadv.adx5791&data=05%7C02%7Csiwanart.ma%40kcl.ac.uk%7C56758459d04241b2a25608de2364125e%7C8370cf1416f34c16b83c724071654356%7C0%7C0%7C638987109879200955%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C&sdata=jCy0uAt07%2FR%2BR5zHWO90O5GzUj9dgboSy0aCoZ3meI4%3D&reserved=0), we project single-cell-resolution gene expression profile to the spatial transcriptomics data. This dashboard provides interactive visualisation to explore our findings.  

Visit our dashboard here: [https://spatial-pancreatic-maps.sites.er.kcl.ac.uk/](https://spatial-pancreatic-maps.sites.er.kcl.ac.uk/)


## Features and Demo

### Single-cell transcriptomics
This scRNA-seq page provides flexible group-by and genes options. You can explore our dataset with UMAP plots, violin plots, dot plots, feature plots, and a cell composition plot.  

<p align="center">
  <img src="https://github.com/user-attachments/assets/1976c595-dec2-4bf7-9d99-45be58cbd9b0"
       alt="scRNA-seq UMAP"
       width="50.5%" />
  <img src="https://github.com/user-attachments/assets/c3034b30-21f6-4822-bde6-e1a43761cf76"
       alt="scRNA-seq Dot Plot"
       width="48%" />
</p>



### Spatial transcriptomics
This page has DAPI microscopic images, spatial scatter plots, neighborhood enrichment plot, and a cell type composition plot. 

<p align="center">
  <img src="https://github.com/user-attachments/assets/f2f2192f-83e8-4db8-a5fe-692ef9a53c86"
       alt="spatial scatter plot"
       width="70%" />
</p>



### Single-cell-resolution gene expression profile projection using Tangram
Cannot find the genes of interest in the Spatial transcriptomics page?  
This is because in our spatial transcriptomic assay, we customised the probes for 50 pancreatic development-related genes, which is far less than the number of genes that scRNA-seq can find.  
However, using Tangram, we project single-cell resolution gene expression profiles onto spatial transcriptomics data.  
Now you can see the predicted spatial location of genes that were not quantified in our spatial transcriptomics assay!  

<p align="center">
  <img src="https://github.com/user-attachments/assets/24a07d83-8355-4a59-9053-a3176b39b081"
       alt="tangram"
       width="80%" />
</p>


## Tech Stack

- Single-cell transcriptomics: `scanpy`
- Spatial transcriptomics: `squidpy`, `Tangram`
- App development: `Streamlit`

## Cite our paper
Torres-Cano, Alejo, et al. "Spatially organized cellular communities shape functional tissue architecture in the pancreas." Science Advances 11.46 (2025): eadx5791.

## Acknowledgements
King's College London. (2022). King's Computational Research, Engineering and Technology Environment (CREATE). Retrieved March 2, 2022, from https://doi.org/10.18742/rnvf-m076

## Feedback
If you have any feedback or questions, please reach out to us by francesca.spagnoli@kcl.ac.uk

## Future plan
We plan to extend the visualisation options with `tissuumaps`
