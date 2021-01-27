# Kozareva data curation (BICCN)

## Project details

Last updated: 

Author: Pablo Garcia-Nieto

Project title: 

DOI: [https://doi.org/10.1038/s41586-018-0654-5](https://doi.org/10.1038/s41586-018-0654-5)

Reqs:

- scanpy
- cellxgene with schema functions
- R and yaml package

Notes:
- Low quality cells were dropped (using the information in the column "class")

## Curation steps

Clone this repo and execute the steps below from the downladed path

- Download data to `./data/original/`:

```bash
bash scripts/all_download_data.sh
```


- Reformat to h5ad and save to `./data/transformed/1_reformatted`

```bash
mkdir -p ./data/transformed/1_reformatted/
python3 ./scripts/reformat_add_metadata.py ./data/original/GSE115746_cells_exon_counts.csv.gz ./data/original/Supplementary_Table_10_Full_Metadata2.txt ./data/original/tsne_tassic_2018_no_NA.csv ./data/transformed/1_reformatted/tassic_2018.h5ad
```

- Create ontology tables, append biccn ontology terms, and create normalize dat

```bash
mkdir -p ./data/transformed/2_normalized
mkdir -p ./data/misc/

python3 ./scripts/create_ontology_lookup_cell_type.py ./data/transformed/1_reformatted/tassic_2018.h5ad ./data/misc/ontology_lookup_cell_type.tsv
python3 ./scripts/create_ontology_lookup_tissue.py ./data/transformed/1_reformatted/tassic_2018.h5ad ./data/misc/ontology_lookup_tissue.tsv
python3 ./scripts/create_biccn_ontology_table.py ./data/transformed/1_reformatted/tassic_2018.h5ad ./data/misc/ontology_biccn.txt

python3 ./scripts/append_biccn_ontology.py ./data/transformed/1_reformatted/tassic_2018.h5ad ./data/misc/ontology_biccn.txt ./data/transformed/2_normalized/temp_biccn.h5ad
python3 ./scripts/rename_columns_anndata.py subclass:BICCN_subclass_label,class:BICCN_class_label ./data/transformed/2_normalized/temp_biccn.h5ad ./data/transformed/2_normalized/temp_biccn_renamed.h5ad
python3 ./scripts/normalize.py ./data/transformed/2_normalized/temp_biccn_renamed.h5ad ./data/transformed/2_normalized/tassic_2018.h5ad

rm ./data/transformed/2_normalized/temp*

```

- Create schema yaml file and apply it to data

```bash
mkdir -p ./data/transformed/3_remixed/

Rscript scripts/create_yaml.R ./data/misc/ontology_lookup_cell_type.tsv ./data/misc/ontology_lookup_tissue.tsv ./schema.yml

cellxgene schema apply --source-h5ad ./data/transformed/2_normalized/tassic_2018.h5ad --remix-config ./schema.yml --output-filename ./data/transformed/3_remixed/tassic_2018.h5ad

cellxgene schema validate ./data/transformed/3_remixed/tassic_2018.h5ad
```
