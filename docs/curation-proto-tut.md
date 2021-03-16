# Cellxgene Curation Tutorial

![image](https://user-images.githubusercontent.com/25663501/111377133-b07a4c80-8676-11eb-8eb8-07ca4d7a77e9.png)


## Introduction
purpose: Motivate submission to the data portal give context

Cellxgene's publishing platform and interactive single cell data explorer is a system which is optimized for easy access, interrogation and sharing/reuse of single cell data. In order to meet these standards and to facilitate easy navigation/aggregation among the growing number of single cell datasets, the cellxgene publishing platform currently accepts curated anndata objects adhering to a [comprehensive data schema (@TODO add link)](https://github.com/chanzuckerberg/single-cell-curation/main/docs/corpora_schema.md). 

In this tutorial, you will gain all the essential information for doing the following: 
- curating one of your own single cell datasets in a streamlined manner
- using CZI curation tools to perform most of the curation work
- uploading data in the cellxgene data portal
- navigating the portal and ensuring public or private access to your datasets
- basic data exploration in the cellxgene explorer view

Beyond this we hope that you see the value of using this portal as a publishing platform and also welcome any comments that pertain improving the portal and curation experience. You can contact us via [support email]()

## General data requirements
- count matrix (required (expcept certain assayss))
- normalized expression matrix used for visulization (required)
- cell level metadata (required)
- variable/feature level metadata (not required)
- embedding (at least one required, UMAP, tSNE, spatial)


## Anndata object structure

<!--- ![image](https://user-images.githubusercontent.com/25663501/111377611-3c8c7400-8677-11eb-8176-cf9de8f64c70.png) --->

<!--- resized anndata image below --->
<img src="https://user-images.githubusercontent.com/25663501/111377611-3c8c7400-8677-11eb-8176-cf9de8f64c70.png" width="500">

- Talk about object structure
- limitations 
- what are the features/aspects of the anndata object that cellxgene data portal expects (essential slots/entries in the anndata object) - talk about where schema information is stored (schema info stored in obs and uns)
- obs
- uns
- obsm ("X_example")
- X
- raw.X or raw layer
- layers

### Conversion to Anndata from other object types (seurat, loom, sce)

Talk about conversion - sceasy, other tools....



## Accepted data types

Talk about the data types that we are able to support (it is mostly ambiguous, but anndata does have some limitations that can be mentioned about here)

Guidelines for curating the following assays. In each section, also include a link to curation of an object of the specified assay type:
- ATAC/mC
- RNAseq
- CITEseq
- Spatial/Visium

## Schema information

As discussed before, schema information included in the uns slot gives us dataset level metadata (ex: doi, paper title). These metadata are stored in a dictionary of key value pairs. The required keys and constraints on their respective values are detailed in the table below:

### Required/suggested dataset level schema fields (i.e. what goes in the `uns` slot?)

- 'layer_descriptions'
- 'organism'
- 'organism_ontology_term_id',
- 'publication_doi'
- 'title'
- 'version'

| Key                         | Value         | Value description  |
| --------------------------- |:-------------:| ------------------:|
| 'version'                   |               |                    |
| 'title'                     |               |                    |
| 'publication_doi'           |               |                    |
| 'organism'                  |               |                    |
| 'organism_ontology_term_id' |               |                    |
| 'layer_descriptions'        |               |                    |


In the obs slot, we store schema information at the individual cell level (important to note here that while original metdata fields here are not altered or replaced, only augmeneted by additional schema specific columns). Each schema field that we require gets two columns in the obs dataframe  (one for the ontology term id, and one for the term itself)

### Required cell level schema fields (i.e. what goes in `obs`?)

| `obs` column name                    | Expected Value| Description        |
| ------------------------------------ |:-------------:| ------------------:|
| 'tissue'                             |               |                    |
| 'tissue_ontology_term_id'            |               |                    |
| 'assay'                              |               |                    |
| 'assay_ontology_term_id'             |               |                    |
| 'disease'                            |               |                    |
| 'disease_ontology_term_id'           |               |                    |
| 'cell_type'                          |               |                    |
| 'cell_type_ontology_term_id'         |               |                    |
| 'sex'                                |               |                    |
| 'ethnicity'                          |               |                    |
| 'ethnicity_ontology_term_id	'        |               |                    |
| 'development_stage'                  |               |                    |
| 'development_stage_ontology_term_id	'|               |                    |
