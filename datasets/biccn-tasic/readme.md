# Kozareva data curation (BICCN)

## Project details

Last updated: Fri Jan  8 20:40:24 2021

Author: Pablo Garcia-Nieto

Project title: Adult mouse cortical cell taxonomy revealed by single cell transcriptomics

DOI: [https://doi.org/10.1038/nn.4216](https://doi.org/10.1038/nn.4216)

Reqs:

- scanpy
- cellxgene with schema functions
- R and yaml package

Notes:

BICCN ontologies:

- This is a study that precedes BICCN, that combined with the limmited capabilities of the BICNN ontology led to a significant misalignment of terms:
    - All L5 glutamatergic neurons were mapped to ILX:0770160 (L5 IT neuron)
    - All L6 glutamatergic neurons were mapped to ILX:0770158 (L6 IT neuron)
    - All Sst neurons (including that and other markers) were mapped to ILX:0770152 (Sst neuron)

## Curation steps

Clone this repo and execute the steps below from the downladed path

- Download data to `./data/original/`:

```bash
bash scripts/all_download_data.sh
```


- Reformat to h5ad and save to `./data/transformed/1_reformatted`

```bash
mkdir -p ./data/transformed/1_reformatted/
python3 ./scripts/reformat_add_metadata.py ./data/original/genes_counts.csv ./data/original/genes_rpkm.csv ./data/original/cell_metadata.csv ./data/transformed/1_reformatted/tassic_Adult_mouse_cortical_cell_taxonomy.h5ad
```

- Create ontology tables, append biccn ontology terms, and create UMAP coordinates

```bash
mkdir -p ./data/transformed/2_normalized
mkdir -p ./data/misc/

python3 ./scripts/create_ontology_lookup_cell_type.py ./data/transformed/1_reformatted/tassic_Adult_mouse_cortical_cell_taxonomy.h5ad ./data/misc/ontology_lookup_cell_type.tsv
python3 ./scripts/create_biccn_ontology_table.py ./data/transformed/1_reformatted/tassic_Adult_mouse_cortical_cell_taxonomy.h5ad ./data/misc/ontology_biccn.txt
python3 ./scripts/append_biccn_ontology.py ./data/transformed/1_reformatted/tassic_Adult_mouse_cortical_cell_taxonomy.h5ad ./data/misc/ontology_biccn.txt ./data/transformed/2_normalized/temp_biccn.h5ad
python3 ./scripts/rename_columns_anndata.py sub_class:BICCN_cluster_label,major_class:BICCN_class_label ./data/transformed/2_normalized/temp_biccn.h5ad ./data/transformed/2_normalized/temp_biccn_renamed.h5ad
python3 ./scripts/umap_anndata.py ./data/transformed/2_normalized/temp_biccn_renamed.h5ad ./data/transformed/2_normalized/tassic_Adult_mouse_cortical_cell_taxonomy.h5ad

rm ./data/transformed/2_normalized/temp*

```

- Create schema yaml file and apply it to data

```bash
mkdir -p ./data/transformed/3_remixed/

Rscript scripts/create_yaml.R ./data/misc/ontology_lookup_cell_type.tsv ./schema.yml

cellxgene schema apply --source-h5ad ./data/transformed/2_normalized/tassic_Adult_mouse_cortical_cell_taxonomy.h5ad --remix-config ./schema.yml --output-filename ./data/transformed/3_remixed/tassic_Adult_mouse_cortical_cell_taxonomy.h5ad

cellxgene schema validate ./data/transformed/3_remixed/tassic_Adult_mouse_cortical_cell_taxonomy.h5ad 
```
