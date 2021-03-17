# Cellxgene Curation Tutorial

<!--- ![image](https://user-images.githubusercontent.com/25663501/111377133-b07a4c80-8676-11eb-8eb8-07ca4d7a77e9.png) --->

<img align="right" src="https://user-images.githubusercontent.com/25663501/111377133-b07a4c80-8676-11eb-8eb8-07ca4d7a77e9.png">




<!--- ## Overview --->
<!--- purpose: Motivate submission to the data portal give context --->

Cellxgene's publishing platform and interactive single cell data explorer is a system which is optimized for easy access, interrogation and sharing/reuse of single cell data. In order to meet these standards and to facilitate easy navigation/aggregation among the growing number of single cell datasets, the cellxgene publishing platform currently accepts curated anndata objects adhering to a [comprehensive data schema (@TODO add link)](https://github.com/chanzuckerberg/single-cell-curation/main/docs/corpora_schema.md). In this tutorial, you will gain all the essential information for doing the following: 
- curating one of your own single cell datasets in a streamlined manner
- using CZI curation tools to perform most of the curation work
- uploading data in the cellxgene data portal
- navigating the portal and ensuring public or private access to your datasets
- basic data exploration in the cellxgene explorer view

Beyond this we hope that you see the value of using this portal as a publishing platform and also welcome any comments that pertain improving the portal and curation experience. You can contact us via [support email]()

<details>
<summary> Table of contents </summary>

- [Why contribute?](#why-contribute)
- [Basic data requirements](#basic-data-requirements)
- [Anndata object requirements](#anndata-object-requirements)
- [Schema definition](#schema-definition)
- [Cellxgene curation tools](#cellxgene-curation-tools)

</details>

<br/>

---

## Why contribute?

<br/>

This is some example text about why to contribute...
<!--- ![image](https://user-images.githubusercontent.com/25663501/111403110-d1589700-86a2-11eb-99dd-bf52e348d43a.png) --->
<br/>

<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/111403110-d1589700-86a2-11eb-99dd-bf52e348d43a.png" width="1000">
</p>

<p align="center">
  <b> Figure: <b> Collection View of the Cellxgene Data Portal
</p>
<br/>

Reason to contribute:
 - link permanence, allows you to reference data portal links in publication with guarantee that the link will always be there
 - sharing of private datasets with collaborators
 - easy for readers to explore your dataset (with no extra work on your side - goodbye shiny!)
 - accesibility of your dataset through many of the popular single cell toolchains (link to an instance of the cellxgene explorer app)

<br/>
 
## Basic data requirements

<br/>

Some basic requirements...
- count matrix (required (expcept certain assays))
- normalized expression matrix used for visulization (required)
- cell level metadata (required)
- variable/feature level metadata (not required)
- embedding (at least one required, UMAP, tSNE, spatial)

<br/>


## Anndata object requirements

<br/>

Introductory text about anndata objects - link to anndata/scanpy

<br/>

<details>
<summary>Anndata description</summary>

<p align="center">
  <img width="500" src="https://user-images.githubusercontent.com/25663501/111377611-3c8c7400-8677-11eb-8176-cf9de8f64c70.png">
</p>

<p align="center">
  <b> Figure: Anndata components <b>
</p>

This is some text introducting anndata


- Talk about object structure
- limitations 
- what are the features/aspects of the anndata object that cellxgene data portal expects (essential slots/entries in the anndata object) - talk about where schema information is stored (schema info stored in obs and uns)
- obs
- uns
- obsm ("X_example")
- X
- raw.X or raw layer
- layers

</details>

<br/>
 
<details>
<summary>Format conversion (i.e. Seurat => Anndata) up</summary>
 
Talk about conversion - sceasy, other tools....
 
</details>

<br/>
 
<details>
<summary>Alternative assays</summary>

Talk about the data types that we are able to support (the tool is mostly ambiguous, but anndata does have some limitations that can be mentioned about here)

Guidelines for curating the following assays. In each section, also include a link to curation of an object of the specified assay type:
- ATAC/mC
- RNAseq
- CITEseq
- Spatial/Visium

</details>

<br/>

## Schema definition

<br/>

The purpose of the cellxgene schema is to support the construction of a data corpus that facilitates data integration across multiple tissues and experiments. This goal requires that we collect a standardized set of metadata about single cell datasets that are uploaded to the cellxgene data portal. To make this process easy to adhere to, we only require a few fields (detailed below) to support easy search and integration across datasets. To access a more comprehensive decription about our schema requirents, you can refer to the [following resource](corpora_schema.md).



<br/>

#### `uns`

`adata.uns` is a dictionary of key values pairs that we use to store dataset level metadata. This is also where you store information on the different data representations that you are sharing (i.e you are sharing a raw counts matrix, normalized expression matrix, scaled and centered expression matrix)

<details>
<summary>click to see required dataset level metadata</summary>

- 'layer_descriptions'
- 'organism'
- 'organism_ontology_term_id',
- 'publication_doi'
- 'title'
- 'version'

| Key                         | Value      | Description        | Example |
| --------------------------- |:----------:| ------------------:| ---|
| 'version'                   | string     |                    | `` |
| 'title'                     | string     |                    | `` |
| 'publication_doi'           | string     |                    | `` |
| 'organism'                  | string     |                    | `` |
| 'organism_ontology_term_id' | string     |                    | `` |
| 'layer_descriptions'        | dictionary | See table below    | `` |


In the obs slot, we store schema information at the individual cell level (important to note here that while original metdata fields here are not altered or replaced, only augmeneted by additional schema specific columns). Each schema field that we require gets two columns in the obs dataframe  (one for the ontology term id, and one for the term itself)

</details>

<br/>

#### `obs`

`adata.obs` is a dataframe that is used to store cell level metadata. In addition to experiemental metadata, the following additional field are required:

<details>
<summary>click to see required cell level metadata</summary>

Table will also include scenrios where ontologies can be relaxed

| `obs` column name                    | Type   | Descritpion                                                   | Example         |
| :----------------------------------- |:------:| -------------------------------------------------------------:|:------------------------------------------|
| 'tissue'                             | string | UBERON term                                                   | `area postrema`                           |
| 'tissue_ontology_term_id'            | string | UBERON term id                                                | `UBERON:0002162`                          |
| 'assay'                              | string | EFO term                                                      | `scRNA-seq`                               |   
| 'assay_ontology_term_id'             | string | EFO term id                                                   | `EFO:0008913`                             |
| 'disease'                            | string | MONDO term or `normal`                                        | `kuru`                                    |
| 'disease_ontology_term_id'           | string | MONDO term id or `PATO:0000461`                               | `MONDO:0006825`                           |
| 'cell_type'                          | string | CL term                                                       | `excitatory neuron`                       |
| 'cell_type_ontology_term_id'         | string | CL term id                                                    | `CL:0008030`                              |
| 'sex'                                | string | `male`, `female`, `mixed`, `unknown`, or `other`              | `mixed`                                   |
| 'ethnicity'                          | string | HANCESTRO term, `na` if non-human, `unknown` if not available | `genetically isolated population`         |
| 'ethnicity_ontology_term_id	'        | string | HANCESTRO term id, `na` if non-human                          | `HANCESTRO:0290`                          |
| 'development_stage'                  | string | HsapDv term, `unknown` if not available                       | `9th week post-fertilization human stage` |
| 'development_stage_ontology_term_id	'| string | HsapDv term id if human, child of `EFO:0000399` otherwise     | `HsapDv:0000046`                          |

</details>

<br/>

## Cellxgene curation tools

<br/>

#### `config.yaml`

<details>
 
<summary> structuring your config file </summary>

What is the config YAML file and how does it modify your anndata object?

</details>

<br/>


#### `cellxgene schema apply`
 
<details>
 
<summary> applying the schema </summary>

Necessary arguments
- input h5ad
- config file
- output h5ad
- warnings and errors

</details>

<br/>

#### `cellxgene schema validate`

<details>
 
<summary> validating the schema up</summary>

Neccesary arguments:
- h5ad to validate

What does a succesful output look like?

</details>

<br/>

### Test `h5ad` in local version of cellxgene

## Uploading data to the cellxgene data portal

Introduce collection view, dataset view, and explorer view

### Portal Sign in

sign in options for the portal (gmail, github, etc)

### Create a collection

`screenshot of creating a collection`

### Add a dataset to a collection

`screenshot adding a dataset to collection`

### Remove dataset from a collection

`screenshot removing dataset`

### Share uploaded datasets with private links
Explain that people can share private datasets 


### Publish Collection to the portal

`Screenshot of Explorer mode`
