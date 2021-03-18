# Cellxgene Curation Tutorial

<!--- ![image](https://user-images.githubusercontent.com/25663501/111377133-b07a4c80-8676-11eb-8eb8-07ca4d7a77e9.png) --->

<img align="right" src="https://user-images.githubusercontent.com/25663501/111377133-b07a4c80-8676-11eb-8eb8-07ca4d7a77e9.png">




<!--- ## Overview --->
<!--- purpose: Motivate submission to the data portal give context --->

Cellxgene's publishing platform and interactive single cell data explorer is a system which is optimized for easy access, interrogation and sharing/reuse of single cell data. In order to meet these standards and to facilitate easy navigation/aggregation among the growing number of single cell datasets, the cellxgene publishing platform currently accepts curated anndata objects adhering to a [comprehensive data schema](corpora_schema.md). In this tutorial, you will gain all the essential information for doing the following: 
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

<!--- ![image](https://user-images.githubusercontent.com/25663501/111403110-d1589700-86a2-11eb-99dd-bf52e348d43a.png) --->
<br/>

<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/111403110-d1589700-86a2-11eb-99dd-bf52e348d43a.png" width="1000">
</p>

<p align="center">
  <b> Figure: </b> Collection View of the Cellxgene Data Portal
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

The following elements of your single cell dataset are required. One of the prerequisites for contributing data to the portal is that you have these components represented as an `anndata` object (more on that in the next section). Note that, in addition to the following, other representations of the expression matrix (alternative normalization, SCTransform, corrected counts from SCTransform or background correction) can all be stored as `layers` your anndata object.

- count matrix (except for certain assays)
- normalized expression matrix used for visulization
- cell level metadata
- variable/feature level metadata (not required)
- embedding (at least one required, UMAP, tSNE, spatial, PCA)



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
  <b> Figure:</b> Anndata components
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

The purpose of the cellxgene schema is to support the construction of a data corpus that facilitates data integration across multiple tissues and experiments. This goal requires that we collect a standardized set of metadata about single cell datasets that are uploaded to the cellxgene data portal. To make this process easy to adhere to, we only require a few fields (detailed below) to support easy search and integration across datasets. To access a more comprehensive decription about our schema requirements, you can refer to the [official schema definition](corpora_schema.md).

#### Basic requirements ([expanded version](corpora_schema.md#basic-requirements))
- **Unique observation identifiers:** Each observation (usually a cell) must have an id that is unique within the dataset.
- **Unique feature identifiers:** Every feature (usually a gene or transcript) also needs a unique identifier.
- **No PII:** No metadata can be personally identifiable: no names, dates of birth, specific locations, etc. There's a [list](https://docs.google.com/document/d/1nlUuRiQ8Awo_QsiVFi8OWdhGONjeaXE32z5GblVCfcI/edit?usp=sharing).


<br/>

#### `uns`

`adata.uns` is a dictionary of key values pairs that we use to store dataset level metadata. This is also where you store information on the different data representations that you are sharing (i.e you are sharing a raw counts matrix, normalized expression matrix, scaled and centered expression matrix)

<details>
<summary>required and suggested dataset level metadata</summary>

<br/>

<div align="center">

| Key                         | Value      | Description        | Example |
| --------------------------- |:----------:| ------------------:| ---|
| 'version'                   | string     |                    | `` |
| 'title'                     | string     |                    | `` |
| 'publication_doi'           | string     |                    | `` |
| 'organism'                  | string     |                    | `` |
| 'organism_ontology_term_id' | string     |                    | `` |
| 'layer_descriptions'        | dictionary | See table below    | `` |

</div>

<div align="center">
  <b>Table: </b> Dataset level metadata
</div>

<br/>

In the above table we see what type on information is necessary at the dataset level. Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum

<br/>

<div align="center">

| Key                         | Value      | Description        | Required                               |
| --------------------------- |:----------:| ------------------:| -------------------------------------- |
| 'raw'                       | string     |                    | required (unless `raw.X` is specified) |
| 'raw.X'                     | string     |                    | required (unless `raw`) is specified   |
| 'X'                         | string     |                    | required (used for visualization)      |
| 'corrected.counts'          | string     |                    | optional                               |
| 'scale.data'                | string     |                    | optional                               |
| 'sct'                       | string     |                    | optional                               |

</div>

<div align="center">
  <b>Table: </b> Required and suggested data layers
</div>

<br/>

In our second table we, see how our different included representations can be described according to the schema standard. In this way, cellxgene explorer and data portal know how to access different parts of the object for visualization and access.


</details>

<br/>

#### `obs`

`adata.obs` is a dataframe that is used to store cell level metadata. In addition to experiemental metadata, the following additional field are required:

<details>
<summary>click to see required cell level metadata</summary>

Table will also include scenrios where ontologies can be relaxed

| `obs` column name                    | Type   | Descritpion                                                       | Example                                   |
| :----------------------------------- |:------:| -----------------------------------------------------------------:|:------------------------------------------|
| 'tissue'                             | string | UBERON term                                                       | `area postrema`                           |
| 'tissue_ontology_term_id'            | string | UBERON term id                                                    | `UBERON:0002162`                          |
| 'assay'                              | string | EFO term                                                          | `scRNA-seq`                               |   
| 'assay_ontology_term_id'             | string | EFO term id                                                       | `EFO:0008913`                             |
| 'disease'                            | string | MONDO term or `normal`                                            | `kuru`                                    |
| 'disease_ontology_term_id'           | string | MONDO term id or `PATO:0000461`                                   | `MONDO:0006825`                           |
| 'cell_type'                          | string | CL term (or UBERON term when no CL term is available)             | `excitatory neuron`                       |
| 'cell_type_ontology_term_id'         | string | CL term id (or UBERON term id when no CL term id is is available) | `CL:0008030`                              |
| 'sex'                                | string | `male`, `female`, `mixed`, `unknown`, or `other`                  | `mixed`                                   |
| 'ethnicity'                          | string | HANCESTRO term, `na` if non-human, `unknown` if not available     | `genetically isolated population`         |
| 'ethnicity_ontology_term_id	'        | string | HANCESTRO term id, `na` if non-human                              | `HANCESTRO:0290`                          |
| 'development_stage'                  | string | HsapDv term, `unknown` if not available                           | `9th week post-fertilization human stage` |
| 'development_stage_ontology_term_id	'| string | HsapDv term id if human, child of `EFO:0000399` otherwise         | `HsapDv:0000046`                          |

<div align="center">
  <b>Table: </b> Required cell level metadata
</div>

<br/>

In the above table, the types of cell level metadata that is collected for the schema. In particular, it is imporant to exapnd on how to annotate `cell_type` field. In general, we try to find most specific ontology term that accurately represents a given cell type in your dataset. In the event that your cell type is not described accurately by any term in the `CL` ontology, then it is sufficient to report the most accurate `UBERON` term and id instead.

</details>

<br/>

## Cellxgene curation tools

<br/>

The cellxgene curation tools include functions that can make the curation process easier for you. Essentially the workflow looks like this:
1. create `config.yaml` that will specify how the schema is to be applied to your anndata object (see example [here](example_config.yaml))
2. run `cellxgene schema apply`, which takes your config and source anndata file to produce a curated anndata object.
3. use `cellxgene schema validate` to ensure that your curated object meets the schema requirements

You can check out the full documenatation for this tooling [here](#schema_guide.md), but a quick introduction to the tools is below:

<br/>

#### Installation

`pip install cellxgene-schema`

<br/>

#### `config.yaml`

<details>
 
<summary> structuring your config file </summary>

<br/>

Your [`config.yaml`](example_config.yaml) file is used to update values and columns in `adata.uns` and `adata.obs` (repectively) with required schema information. Here is an example/mock of a config file:

##### `uns` section

```
uns:
    version:
        corpora_schema_version: 1.1.0                           #(ex: schema_version_number)
        corpora_encoding_version: 1.1.0                         #(ex: encoding version)
    title: publication_title                                    #(free text field)
    layer_descriptions:
        raw.X: raw                                              #(it is essential for at one the layer_descriptions to be set to 'raw')
        X: log1p                                                #(free text description of normalization method )
    organism: Homo sapiens                                      #(Specify organism name)
    organism_ontology_term_id: NCBITaxon:9606                   #(#Specfiy organism NCBI Taxon ID)

```


<br/>

In the config file snippet above, our 0 level indendentation specifies that we are modifying the `uns` slot. The next level of indentation specifies a key name to add to `uns`. The key's corresponding value will the be string following the colon, or if the key is follow by more lines that are further indented, then the corresponding value will be a dictonary containing the key-value pairs specified in the subsequent lines. More concretely, `adata.uns['layer_descriptions']`, will return a dictionary with the key value pairs `{'raw.X': 'raw', 'X': 'log1p'}` after the schema config has been applied.


**Note:** at least one of the keys in `layer_descriptions` must return the value 'raw'

<br/>

##### `obs` section
```
obs:
    tissue_ontology_term_id: UBERON:0000006                                 #UBERON tissue term
    assay_ontology_term_id: EFO:0010961                                     #EFO assay term
    disease_ontology_term_id: PATO:0000461                                  #MONDO disease term or PATO:0000461
    sex: male                                                               #male, female, mixed, unknown, or other
    ethnicity_ontology_term_id: na                                          #HANCESTRO term
    development_stage_ontology_term_id: HsapDv:0000160                      #HsapDv term
    cell_type_ontology_term_id:
      cell_label:                                                           #this is column name in your obs dataframe (the field that specifies cell type)
        acinar: CL:0002064                                                  #mapping between cell types (specified in anndata object) and cl ontology id
        alpha: CL:0000171
        delta: CL:0000173
        endothelial: CL:0000115
        epsilon: CL:0005019

```


<br/>

In the config file snippet above, the 0 level indentation specifies that we are modifying `obs`. At the first level of indentation, we start adding schema fields to `obs`. If the first level indentation fields ar follow by an ontology value (i.e. `tissue_ontology_term_id: UBERON:0000006`), then a column called `tissue_ontology_term_id` will created in the obs dataframe, and all of its values will be `UBERON:0000006`. Additionally, another column named `tissue` will be created with the corresponding UBERON term (`islet of Langerhans`).

On other hand, if the schema field is followed by lines with have further indentation, then the next line in the config file specify a column whose values should be mapped to ontology terms (with these mappings specified by the subsequent lines). Two new columns will be added to the `obs` data frame `cell_type_ontology_term_id` and `cell_type`.

<br/>

@TODO - the fixup gene symbols section

</details>

<br/>


#### `cellxgene schema apply`
 
<details>
 
<summary> applying the schema </summary>

<br/>

In order to use the `cellxgene schema apply` command, you will need to pass the following required arguments:
- `--source-h5ad` your original anndata file
- `--remix-config` the `config.yaml` file that we specified above 
- `--output-filename` the name of the resulting anndata that is consistent with the cellxgene schema

<br/>

The next step will be to validate the resulting object.

</details>

<br/>

#### `cellxgene schema validate`

<details>
 
<summary> validating the schema</summary>

<br/>

In order to validate the remixed object, needs to simply run `cellxgene schema validate remixed_anndata.h5ad`. If there has been no terminal output from the function, then your object has been validated successfuly and is ready for upload!

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
