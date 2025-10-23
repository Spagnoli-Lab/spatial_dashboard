
# Spatial pancreatic atlas dashboard

In our paper, we project single-cell-resolution gene expression profile to the spatial transcriptomics data. This dashboard provides interactive visualisation to explore our findings. 


## Demo

- Select samples of interest
- Select annotation
- Adjust image size and point size
- Choose genes of interest through the drop-down menu or typing

## Features

- Single-cell transcriptomics: 
  - UMAP
  - Violin plots
  - Dot plot
  - Feature plot
  - Cell type proportion plot
- Spatial transcriptomics: 
  - DAPI image
  - Spatial scatter plot
  - Cell type proportion plot
- Single-cell-resolution gene expression profile projection:
  - Spatial scatter plot


## Tech Stack

- Single-cell transcriptomics: `scanpy`

- Spatial transcriptomics: `squidpy`, `Tangram`

- App development: `Streamlit`


## Acknowledgements

 - King's College London. (2022). King's Computational Research, Engineering and Technology Environment (CREATE). Retrieved March 2, 2022, from https://doi.org/10.18742/rnvf-m076

## Authors

- [@SiwanartMa](https://github.com/SiwanartMa)


## Feedback

If you have any feedback, please reach out to us by francesca.spagnoli@kcl.ac.uk


## Future plan
We plan to extend the visualisation options with `tissuumaps`