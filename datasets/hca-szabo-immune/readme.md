# Tasic 2018 data curation (BICCN)

## Project details

Last updated: Fri Feb 26 16:22:28 2021

Author: Pablo Garcia-Nieto

Project title: Single-cell transcriptomics of human T cells reveals tissue and activation signatures in health and disease

DOI: [https://doi.org/10.1038/s41467-019-12464-3](https://doi.org/10.1038/s41467-019-12464-3)

Reqs:

- scanpy
- cellxgene-schema CLI
- R and yaml package

Notes:
- No normalization on counts performed on original dataset for 2D-projection

## Curation steps

Clone this repo and execute the steps below from the downladed path

- Download data to `./data/original/`:

    - Count matrices are available from GEO [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126030](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126030)
    - Metadata is in sheet `Fig 6` from the paper source data. A processed file is provided with this repos, so no need to redownload. [https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-12464-3/MediaObjects/41467_2019_12464_MOESM9_ESM.xlsx](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-12464-3/MediaObjects/41467_2019_12464_MOESM9_ESM.xlsx)
    
    
- Uncompress matrices and join them into a single one:

```bash
tar -zvf GSE126030_RAW.tar
python3 ./scripts/join_tables.py
```


- Reformat to h5ad, append metadata and save to `./data/transformed/1_reformatted`

```bash
mkdir -p ./data/transformed/1_reformatted/
python3 scripts/reformat_add_metadata.py
```

- Normalize to log(CPM+1)

```bash
mkdir -p ./data/transformed/2_normalized
python3 scripts/normalize.py ./data/transformed/1_reformatted/hca_immune_szabo.h5ad ./data/transformed/2_normalized/hca_immune_szabo.h5ad
```

- Create schema yaml file and apply it to data

```bash
mkdir -p ./data/transformed/3_remixed/

Rscript scripts/create_yaml.R ./data/misc/ontology_lookup_* ./schema.yml

cellxgene schema apply --source-h5ad ./data/transformed/2_normalized/hca_immune_szabo.h5ad --remix-config ./schema.yml --output-filename ./data/transformed/3_remixed/hca_immune_szabo.h5ad

cellxgene schema validate ./data/transformed/3_remixed/hca_immune_szabo.h5ad
```
