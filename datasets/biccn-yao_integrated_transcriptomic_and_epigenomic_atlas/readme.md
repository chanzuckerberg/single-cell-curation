# Yao et al. data curation (BICCN)

## Project details

Last updated: Tue Dec  8 00:11:37 2020

Author: Pablo Garcia-Nieto

Project title: An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types

DOI: [https://doi.org/10.1101/2020.02.29.970558](https://doi.org/10.1101/2020.02.29.970558)

Reqs:
- scanpy
- cellxgene with schema functions


Notes:
- There's RNA-seq, ATAC-seq and methylation sequencing data
- Getting all the infomation needed to curate the RNA-seq data required manual extraction of infomation from paper, getting in touch with authors, and getting in touch with collaborators at the Allen Institute to clarify brain ontology questions.
- This study was a used a reference to start a [brain-centric ontology](https://github.com/SciCrunch/NIF-Ontology/blob/master/docs/Neurons.md)

## Curation steps

Clone this repo and execute the steps below from the downloaded path

- Download data to `./data/original/`

```bash
bash scripts/all_download_data.sh
```

- Reformat to h5ad and save to `./data/transformed/1_reformatted`

```bash
bash scripts/all_reformat_data.sh
```

- Clean data (eliminate low-quality cells and rename BICCN columns). Results are saved to `./data/transformed/2_cleaned`

```bash
bash scripts/all_clean.sh
```

- Collect cell type information for later use

```bash
bash scripts/all_collect_cell_types.sh
```

- Add brain centric ontology to `./data/transformed/2_1_cleaned_ontolgy`

```bash
bash scripts/all_append_biccn_ontology.sh
```

- Normalize data and save to `./data/transformed/3_prepared`

```bash
scripts/all_prepare.sh
```

- Create ontology mapping, make schema yaml files and apply it to files, save to `./data/transformed/4_remixed`

```
bash scripts/all_apply_schema.sh
```

- Concatenate into a single file `./data/transformed/4_remixed/10X_cells_v3_AIBS_filtered.h5ad

```
mkdir -p ./data/transformed/5_merged/
python3 scripts/all_concatenate.py ./data/transformed/4_remixed/ ./data/transformed/5_merged/yao_biccn_mop_merged.h5ad
```
