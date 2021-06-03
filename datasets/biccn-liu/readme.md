# BICCN Liu et al. epigenomic data 

## Project details

Last updated: Thu Jun  3 16:33:12 2021

Author: Pablo Garcia-Nieto

Project title: DNA Methylation Atlas of the Mouse Brain at Single-Cell Resolution

DOI: [https://doi.org/10.1101/2020.04.30.069377](https://doi.org/10.1101/2020.04.30.069377)

Reqs:

- scanpy
- cellxgene with schema functions
- R and yaml package
- gdown


## Curation steps

Clone this repo and execute the steps below from the downladed path

- Download data to `./data/original/`:

```bash
bash scripts/download.sh
```



- Append BICCN ontology terms

```bash
mkdir -p ./data/transformed/2_reformatted

for i in ./data/original/*h5ad
do
    python3 ./scripts/append_biccn_ontology.py $i ./data/misc/ontology_biccn.txt ./data/transformed/2_reformatted/temp_biccn.h5ad
    python3 ./scripts/remove_nan_fix_embeddings.py ./data/transformed/2_reformatted/temp_biccn.h5ad ./data/transformed/2_reformatted/temp_biccn_2.h5ad
    python3 ./scripts/rename_columns_anndata.py MajorType:BICCN_subclass_label,CellClass:BICCN_class_label,SubType:BICCN_cluster_label ./data/transformed/2_reformatted/temp_biccn_2.h5ad ./data/transformed/2_reformatted/$(basename $i)
    rm ./data/transformed/2_reformatted/temp*
done

```

- Create schema yaml file and apply it to data

```bash
mkdir -p ./data/transformed/3_remixed/

Rscript scripts/create_yaml.R "DNA Methylation (CGN) Atlas of the Mouse Brain at Single-Cell Resolution" ./data/misc/ontology_lookup_* ./schema_CGN.yml
Rscript scripts/create_yaml.R "DNA Methylation (CHN) Atlas of the Mouse Brain at Single-Cell Resolution" ./data/misc/ontology_lookup_* ./schema_CHN.yml


cellxgene schema apply --source-h5ad ./data/transformed/2_reformatted/MouseBrainMethylome.mCG.h5ad --remix-config ./schema_CGN.yml --output-filename ./data/transformed/3_remixed/MouseBrainMethylome.mCG.h5ad

cellxgene schema apply --source-h5ad ./data/transformed/2_reformatted/MouseBrainMethylome.mCH.h5ad --remix-config ./schema_CHN.yml --output-filename ./data/transformed/3_remixed/MouseBrainMethylome.mCH.h5ad

for i in ./data/transformed/3_remixed/*h5ad
do
    cellxgene schema validate $i
done 

```
