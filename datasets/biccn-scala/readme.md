# BICCN Scala et al. Patch-seq

## Project details

Last updated: Wed Jun  2 15:01:22 2021

Author: Pablo Garcia-Nieto

Project title: Phenotypic variation of transcriptomic cell types in mouse motor cortex

DOI: [https://doi.org/10.1038/s41586-020-2907-3](https://doi.org/10.1038/s41586-020-2907-3)

Reqs:

- scanpy
- cellxgene with schema functions
- R and yaml package


Notes:
- This is Patch-seq dataset.
- To annotate the cells, the data was integrated with 10X data from the flagship BICCN data. The integrated is processed here as well.


## Curation steps

Clone this repo and execute the steps below from the downladed path

- Download data to `./data/original/`:

```bash
bash scripts/download.sh
```

- Arrange patch-seq data into h5ads, this will produce 6 files, 3 of them are the cells separeted based on patch metadata (morphological, electophysiological, and both), the other 3 are the integrated subsets with 10X data (glutamatergic, and two GABAergic subsets).

Files will be saved at `./data/transformed/1_refformated/`

```
mkdir -p ./data/transformed/1_refformated/
python3 ./scripts/make_h5ad_patchseq.py
```

- Arrange 10X data from flagship paper, this will produce 3 files corresponding to the patch-seq files (glutamatergic, and two GABAergic subsets)

```
mkdir -p ./data/transformed/1_refformated/
mkdir -p ./data/transformed/2_fixed_genes/
python3 ./scripts/make_h5ad_allen.py

for i in ./data/transformed/1_refformated/10X_cells*
do
    cellxgene-schema fixup-genes --source-h5ad $i --remix-config ./schema_files/10x.yml --output-filename ./data/transformed/2_fixed_genes/$(basename $i)
done
```


- Append BICCN ontology terms to the 6 patch-seq datasets, files will be saved to `./data/transformed/2_biccn_ontologies`

```bash
mkdir -p ./data/transformed/2_biccn_ontologies

for i in ./data/transformed/1_refformated/patchseq*
do
    python3 ./scripts/append_biccn_ontology.py $i ./data/misc/ontology_biccn.txt ./data/transformed/2_biccn_ontologies/temp_biccn.h5ad
    python3 ./scripts/add_biccn_subclass_convert_type.py ./data/transformed/2_biccn_ontologies/temp_biccn.h5ad ./data/misc/ontology_biccn.txt ./data/transformed/2_biccn_ontologies/temp_biccn2.h5ad
    python3 ./scripts/rename_columns_anndata.py "RNA type:BICCN_cluster_label" ./data/transformed/2_biccn_ontologies/temp_biccn2.h5ad ./data/transformed/2_biccn_ontologies/$(basename $i)
    rm ./data/transformed/2_biccn_ontologies/temp*
done
```

- Create schema yaml files and apply it to data

```bash
mkdir -p ./schema_files
mkdir -p ./data/transformed/3_remixed/

Rscript scripts/create_yaml.R "CGE-derived interneurons integrated with 10X sequencing MOp data" ./data/misc/ontology_lookup_* ./schema_files/patchseq_nonMorpho_Lamp5Vip.h5ad.yml
Rscript scripts/create_yaml.R "MGE-derived interneurons integrated with 10X sequencing MOp data" ./data/misc/ontology_lookup_* ./schema_files/patchseq_nonMorpho_PvalbSst.h5ad.yml
Rscript scripts/create_yaml.R "Excitatory neurons integrated with 10X sequencing MOp data" ./data/misc/ontology_lookup_* ./schema_files/patchseq_nonMorpho_Excitatory.h5ad.yml

Rscript scripts/create_yaml.R "All cells with electrophysiological recordings" ./data/misc/ontology_lookup_* ./schema_files/patchseq_electro_phys.h5ad.yml
Rscript scripts/create_yaml.R "All cells with morphological recordings" ./data/misc/ontology_lookup_* ./schema_files/patchseq_morphological.h5ad.yml
Rscript scripts/create_yaml.R "All cells with electrophysiological and morphological recordings" ./data/misc/ontology_lookup_* ./schema_files/patchseq_morpholphys.h5ad.yml

for i in ./data/transformed/2_biccn_ontologies/*
do
    cellxgene-schema apply --source-h5ad $i --remix-config ./schema_files/$(basename $i).yml --output-filename ./data/transformed/3_remixed/$(basename $i)
done
```

- Merge 10x data to corresponding patch seq
```bash
python3 ./scripts/merge_patchseq_10x.py ./data/transformed/3_remixed/patchseq_nonMorpho_Excitatory.h5ad ./data/transformed/2_fixed_genes/10X_cells_v2_AIBS_exc.h5ad ./data/misc/ontology_biccn_yao.txt ./data/transformed/3_remixed/patchseq_nonMorpho_Excitatory_with10x.h5ad

python3 ./scripts/merge_patchseq_10x.py ./data/transformed/3_remixed/patchseq_nonMorpho_Lamp5Vip.h5ad ./data/transformed/2_fixed_genes/10X_cells_v2_AIBS_viplamp.h5ad ./data/misc/ontology_biccn_yao.txt ./data/transformed/3_remixed/patchseq_nonMorpho_Lamp5Vip_with10x.h5ad

python3 ./scripts/merge_patchseq_10x.py ./data/transformed/3_remixed/patchseq_nonMorpho_PvalbSst.h5ad ./data/transformed/2_fixed_genes/10X_cells_v2_AIBS_pvsst.h5ad ./data/misc/ontology_biccn_yao.txt ./data/transformed/3_remixed/patchseq_nonMorpho_PvalbSst_with10x.h5ad
```
