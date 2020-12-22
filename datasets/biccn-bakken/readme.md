# Bakken data curation (BICCN)

## Project details

Last updated: Mon Dec 21 22:26:10 2020

Author: Pablo Garcia-Nieto

Project title: Evolution of cellular diversity in primary motor cortex of human, marmoset monkey, and mouse.

DOI: [https://doi.org/10.1101/2020.03.31.016972](https://doi.org/10.1101/2020.03.31.016972)

Reqs:

- sceasy
- scanpy
- cellxgene with schema functions
- R and yaml package

Note:
- These are 4 datasets, each contains a different major cell type and has integrated counts from human, mouse, marmoset, and macaque.
- The integration of multiple datasets has raised concerns about our current cellxgene schema -- detailed [here](https://github.com/chanzuckerberg/single-cell-curation/issues/7)

## Curation steps

Clone this repo and execute the steps below from the downladed path

- Download data. It will be saved to `./data/original/` 

```bash
bash scripts/all_download_data.sh
```

- Reformat from Seurat to h5ad and save to `./data/transformed/1_reformatted`

```bash
bash scripts/all_reformat_data.sh
```

- Append information, add ontolgies, and rename columns. This will:

    1. Create ontology mapping tables later used to create yaml file for cellxgene
    2. Appends biccn ontolgies
    3. Corrects sex labels
    4. Appends organism labels and ontologies
    5. Renames certain columns to adhere to cellxgene schema
    
Resulting datasets will be on `./data/transformed/2_clean`

```bash
scripts/all_ontology.bash
```

- Create schema yaml file and apply it to data. Remixed files will be on `./data/transformed/3_remixed/`

```bash
bash scripts/all_apply_schema.sh
```
