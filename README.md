
# Spatial pancreatic atlas dashboard

In our paper [Spatially organized cellular communities shape functional tissue architecture in the pancreas](https://eur03.safelinks.protection.outlook.com/?url=https%3A%2F%2Fdoi.org%2F10.1126%2Fsciadv.adx5791&data=05%7C02%7Csiwanart.ma%40kcl.ac.uk%7C56758459d04241b2a25608de2364125e%7C8370cf1416f34c16b83c724071654356%7C0%7C0%7C638987109879200955%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C&sdata=jCy0uAt07%2FR%2BR5zHWO90O5GzUj9dgboSy0aCoZ3meI4%3D&reserved=0), we project single-cell-resolution gene expression profile to the spatial transcriptomics data. This dashboard provides interactive visualisation to explore our findings. 


## Features and Demo

- Single-cell transcriptomics: 

![scRNA-seq UMAP](https://github.com/user-attachments/assets/1976c595-dec2-4bf7-9d99-45be58cbd9b0)

![scRNA-seq Violin Plot](https://github.com/user-attachments/assets/ef52cb8a-85d9-4907-8696-ab9d21e194f0)

![scRNA-seq Dot Plot](https://github.com/user-attachments/assets/c3034b30-21f6-4822-bde6-e1a43761cf76)

![scRNA-seq Feature Plot](https://github.com/user-attachments/assets/1cf3f688-1e01-46fb-9fc3-4dd802525d44)

- Spatial transcriptomics: 

![spatial DAPI and Spatial Scatter Plot](https://github.com/user-attachments/assets/f2f2192f-83e8-4db8-a5fe-692ef9a53c86)

- Single-cell-resolution gene expression profile projection using Tangram:

![tangram](https://github.com/user-attachments/assets/24a07d83-8355-4a59-9053-a3176b39b081)

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
