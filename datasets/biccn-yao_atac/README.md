# Yao data curation (BICCN)

## Project details

Last updated: Fri Jan 15th 2021

Author: Maximilian Lombardo

Project title: An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types

DOI: [https://doi.org/10.1101/2020.02.29.970558](https://doi.org/10.1101/2020.02.29.970558)

Reqs:

- scanpy
- pandas
- anndata
- cellxgene with schema functions

Note:
- Data from this study were taken from Li2020 (only a subset) and used for further integration with additional transcriptomic data - this section of the repo describes how to obtain the correct subset of sNuc profiles from the Li dataset and create objects that were reported in the Yao Integration paper. For more information check out the `./datasets/biccn-li` directory in this repo

## Curation steps

Clone this repo and download Li objects (GABA, GLUTA, NonN) from (Locaiton placeholder)

### Obtain Yao metadata

Navigate to the `./code/` directory, and run the download script to obtain the appropriate metadata (specifically the Yao embeddings for UMAP and tSNE). You can then move the downloaded file to the appropriate data folder (also a child of this main directory)

```bash
bash metadataDownload.sh
```



## Object Creation and annotation

You can then run the `YaoCuration.py` script which accomplishes a few things:

- reconstructs the barcodes for the Yao dataset (`sample` + `cellID`) so that they match the naming scheme of the Li objects
- check that there are merged cellIDs from the Yao metadata that match the IDs in the Li datasets
- Find the subsets from each of the Li objects that have entries in the Yao metadata and merge those into a single `.h5ad` object
- Make observation names unique
- Update the merged `.h5ad` object with embeddings from the Yao paper
- Do some basic plotting to make sure that the object reflects what was reported in the paper
- write out the object


## Apply cellxgene schema to newly created objects

Cellxgene schema can be applied by using the provided config file (found in `./datasets/biccn-yao_atac/objects/config.yaml`) and by running the following commands in the`object` directory

```bash
cellxgene schema apply --source-h5ad Yao2020ATAC_allCells_uniqueObsNames.h5ad --remix-config config.yaml --output-filename Yao2020ATAC_allCells_remixed.h5ad 
cellxgene schema validate Yao2020ATAC_allCells_remixed.h5ad
```


