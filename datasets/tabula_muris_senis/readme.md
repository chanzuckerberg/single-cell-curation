# Tabula muris senis  data curation

## Project details

Last updated: Mon Feb  8 15:48:37 2021

Author: Pablo Garcia-Nieto

Project title: A single-cell transcriptomic atlas characterizes ageing tissues in the mouse

DOI: [https://doi.org/10.1038/s41586-020-2496-1](https://doi.org/10.1038/s41586-020-2496-1)

Reqs:

- scanpy
- aws cli
- cellxgene with schema functions
- R and yaml package

Notes:
- Cell ontology ids from the original h5ads seem to be incorrect. I programtically looked up the ids using `cell_ontology_class` from `anndata.obs`

## Curation steps

Clone this repo and execute the steps below from the downladed path

- Download data to `./data/original/`:

```bash
mkdir -p ./data/original

aws s3 cp s3://czb-tabula-muris-senis/Data-objects/ ./data/original --recursive --exclude "*" --include "*h5ad"

```


- Look up ontology ids for cell types and tissues. The file `./data/misc/ontology_lookup_tissue_DO_MANUAL_CURATION.txt` created in this step has to be manually curated, i.e. double-check that the values from the first column correspond to the ontology id, if not correct accordingly using EBI's ontology look up service. The final file should be `./data/misc/ontology_lookup_tissue.txt`

```bash
mkdir -p ./data/misc/

# Tissue
python3 ./scripts/create_ontology_lookup_table.py tissue tissue ./data/original/tabula-muris-senis-facs-processed-official-annotations.h5ad ./data/misc/ontology_lookup_tissue_DO_MANUAL_CURATION.tsv

# Cell types
python3 ./scripts/create_ontology_lookup_table.py cell_ontology_class cell_type ./data/original/tabula-muris-senis-droplet-processed-official-annotations.h5ad ./data/misc/droplet_ontology_lookup_cell_type.tsv
python3 ./scripts/create_ontology_lookup_table.py cell_ontology_class cell_type ./data/original/tabula-muris-senis-facs-processed-official-annotations.h5ad ./data/misc/facs_ontology_lookup_cell_type.tsv

cat <(tail -n +2 ./data/misc/droplet_ontology_lookup_cell_type.tsv) <(tail -n +2 ./data/misc/facs_ontology_lookup_cell_type.tsv) | sort | uniq > ./data/misc/temp
cat <(head -n 1 ./data/misc/droplet_ontology_lookup_cell_type.tsv) ./data/misc/temp > ./data/misc/ontology_lookup_cell_type.tsv

rm ./data/misc/temp ./data/misc/droplet_ontology_lookup_cell_type.tsv ./data/misc/facs_ontology_lookup_cell_type.tsv

```

- Create yaml with schema mapping and apply it to files

```bash
mkdir -p ./data/transformed/1_remixed/

Rscript scripts/create_yaml.R ./data/misc/ontology_lookup_cell_type.tsv ./data/misc/ontology_lookup_tissue.tsv ./data/misc/ontology_lookup_developmental_stage.tsv ./schema_droplet.yml
Rscript scripts/create_yaml_facs.R ./data/misc/ontology_lookup_cell_type.tsv ./data/misc/ontology_lookup_tissue.tsv ./data/misc/ontology_lookup_developmental_stage.tsv ./schema_facs.yml

for i in ./data/original/*droplet*h5ad
do
    echo $i 
    cellxgene schema apply --source-h5ad $i --remix-config ./schema_droplet.yml --output-filename ./data/transformed/1_remixed/$(basename $i)
done

for i in ./data/original/*facs*h5ad
do
    echo $i 
    cellxgene schema apply --source-h5ad $i --remix-config ./schema_facs.yml --output-filename ./data/transformed/1_remixed/$(basename $i)
done

# Remove one column repeated ontology column and validate
for i in ./data/transformed/1_remixed/*h5ad
do
    echo $i 
    python3 ./scripts/remove_obs_columns.py cell_ontology_class,cell_ontology_id $i 
    cellxgene schema validate $i
done
```
