# Zheng data curation (BICCN)

## Project details

Last updated: Tue Feb  2 11:12:11 2021

Author: Pablo Garcia-Nieto

Project title: Molecular, spatial and projection diversity of neurons in primary motor cortex revealed by in situ single-cell transcriptomics

DOI: [https://doi.org/10.1101/2020.06.04.105700](https://doi.org/10.1101/2020.06.04.105700)

Reqs:

- scanpy
- cellxgene with schema functions
- R and yaml package

Notes:


## Curation steps

Clone this repo and execute the steps below from the downladed path

- Download data to `./data/original/`:

```bash
bash scripts/all_download_data.sh
```


- Reformat to h5ad and save to `./data/transformed/1_reformatted`

```bash
mkdir -p ./data/transformed/1_reformatted/
python3 ./scripts/reformat_add_metadata.py ./data/original/counts.h5ad ./data/original/cell_metadata.csv ./data/original/umap_embedding.csv ./data/transformed/1_reformatted/zheng_biccn_merfish.h5ad
```

- Append BICCN ontologies and add rename columns

```bash

mkdir -p ./data/transformed/2_renamed/

python3 ./scripts/append_biccn_ontology.py ./data/transformed/1_reformatted/zheng_biccn_merfish.h5ad ./data/misc/ontology_biccn.tsv ./data/transformed/2_renamed/temp_biccn.h5ad
python3 ./scripts/rename_columns_anndata.py subclass:BICCN_subclass_label,class_label:BICCN_class_label,label:BICCN_cluster_label ./data/transformed/2_renamed/temp_biccn.h5ad ./data/transformed/2_renamed/zheng_biccn_merfish.h5ad

rm ./data/transformed/2_renamed/temp*
```

- Create schema yaml file and apply it to data

```bash
mkdir -p ./data/transformed/3_remixed/

Rscript scripts/create_yaml.R ./data/misc/ontology_lookup_cell_type.tsv ./schema.yml

cellxgene schema apply --source-h5ad ./data/transformed/2_renamed/zheng_biccn_merfish.h5ad --remix-config ./schema.yml --output-filename ./data/transformed/3_remixed/zheng_biccn_merfish.h5ad

cellxgene schema validate ./data/transformed/3_remixed/zheng_biccn_merfish.h5ad
```
