# Kozareva data curation (BICCN)

## Project details

Last updated: Wed Dec  9 21:57:30 2020

Author: Pablo Garcia-Nieto

Project title: An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types

DOI: [https://doi.org/10.1101/2020.02.29.970558](https://doi.org/10.1101/2020.02.29.970558)

Reqs:
- scanpy
- cellxgene with schema functions


## Curation steps

Clone this repo and execute the step below from the repo path

- Download data to `./data/original/` from [https://singlecell.broadinstitute.org/single_cell/study/SCP795/a-transcriptomic-atlas-of-the-mouse-cerebellum#study-summary](https://singlecell.broadinstitute.org/single_cell/study/SCP795/a-transcriptomic-atlas-of-the-mouse-cerebellum#study-summary)


- Reformat to h5ad and save to `./data/transformed/1_reformatted`

```bash
mkdir -p ./data/transformed/1_reformatted
python3 ./scripts/reformat_add_metadata.py
```

- Normalize data, create ontology tables and append biccn ontology terms

```bash
mkdir -p ./data/transformed/2_normalized
mkdir -p ./data/misc/

python3 ./scripts/normalize.py ./data/transformed/2_normalized/kozareva_transcriptomic_atlas_of_the_mouse_cerebellum.h5ad ./data/transformed/2_normalized/temp_norm.h5ad

python3 ./scripts/create_ontology_lookup_cell_type.py ./data/transformed/2_normalized/temp_norm.h5ad ./data/misc/ontology_lookup_cell_type.tsv
python3 ./scripts/create_biccn_ontology_table.py ./data/transformed/2_normalized/temp_norm.h5ad ./data/misc/ontology_biccn.txt
python3 ./scripts/append_biccn_ontology.py ./data/transformed/2_normalized/temp_norm.h5ad ./data/misc/ontology_biccn.txt ./data/transformed/2_normalized/temp_biccn.h5ad
python3 ./scripts/rename_columns_anndata.py cell_type__custom:BICCN_subclass_label ./data/transformed/2_normalized/temp_biccn.h5ad ./data/transformed/2_normalized/kozareva_transcriptomic_atlas_of_the_mouse_cerebellum.h5ad

rm ./data/transformed/2_normalized/temp*

```

- Create schema yaml file and apply it to data

```bash
mkdir -p ./data/transformed/3_remixed/

Rscript scripts/create_yaml.R ./data/misc/ontology_lookup_cell_type.tsv ./schema.yml

cellxgene schema apply --source-h5ad ./data/transformed/2_normalized/kozareva_transcriptomic_atlas_of_the_mouse_cerebellum.h5ad --remix-config ./schema.yml --output-filename ./data/transformed/3_remixed/kozareva_transcriptomic_atlas_of_the_mouse_cerebellum.h5ad

cellxgene schema validate ./data/transformed/3_remixed/kozareva_transcriptomic_atlas_of_the_mouse_cerebellum.h5ad 
```
