# Cellxgene Curation Tutorial

<!--- ![image](https://user-images.githubusercontent.com/25663501/111377133-b07a4c80-8676-11eb-8eb8-07ca4d7a77e9.png) --->

<img align="right" src="https://user-images.githubusercontent.com/25663501/111377133-b07a4c80-8676-11eb-8eb8-07ca4d7a77e9.png">


<!--- ## Overview --->
<!--- purpose: Motivate submission to the data portal give context --->

[Cellxgene's publishing platform](https://cellxgene.cziscience.com/) and [interactive single cell data explorer](https://github.com/chanzuckerberg/cellxgene) is a system which is optimized for access, exploration and reuse of single cell data. In order to achieve these goals, the cellxgene platform currently accepts curated [anndata](https://anndata.readthedocs.io/en/latest/#) objects adhering to a [succinct data schema](corpora_schema.md). Adherence to a standardized data schema allows for efficient navigation and integration of the growing number of single cell datasets that are becoming available. In this tutorial, you will learn the essential information for curating a single cell dataset using CZI's curation tools and uploading to the data portal. Hosting your data on the cellxgene portal will offer the following benefits:
 - link permanence (you can reference in your publication without ever worrying about dead links)
 - sharing of private datasets with collaborators (keep the data private unitl it is ready for publication)
 - no barrier for readers to explore your dataset (and no need for you to build your own single cell data explorer)
 - accesibility of your dataset in the major single cell data formats (including `anndata`, `seurat`, and `loom`)

This tutorial will consist of an explanation of 1) how to create and structure an `anndata` with your single cell data 2) how to augment this object with information that is specific to the cellxgene schema and 3) how to upload this object to the cellxgene data portal. If you run into any issues during this tutorial, or have any suggestions on how to improve the portal and curation experience, you can contact us via [cellxgene@chanzuckerberg.com](cellxgene@chanzuckerberg.com).

<!---
You can also check out the following links for more information on the cellxgene ecosytem:

- [cellxgene data portal](https://cellxgene.cziscience.com/)
- [cellxgene explorer github](https://github.com/chanzuckerberg/cellxgene)
- [cellxgene data portal github](https://github.com/chanzuckerberg/corpora-data-portal)
- [corpora schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/docs/corpora_schema.md)
- [schema guide](https://github.com/chanzuckerberg/single-cell-curation/blob/main/docs/schema_guide.md)

--->

### Table of Contents

- [Quick start](#quick-start)
- [Data requirements and anndata structure](#required-data-and-anndata-structure)
  - [Format conversion](#format-conversion)
  - [Alternative assays](#alternative-assays) 
- [Schema definition](#schema-definition)
  - [Basic requirements](#basic-requirements-expanded-version)
  - [`uns`](#uns)
  - [`obs`](#obs) 
- [Cellxgene curation tools](#cellxgene-curation-tools)
  - [Installation](#installation)
  - [`config.yaml`](#configyaml)
    - [`uns` section](#uns-section)
    - [`obs` section](#obs-section)
    - [`fixup_gene_symbols` section](#fixup_gene_symbols-section-gene-symbol-harmonization)
  - [`cellxgene schema apply`](#cellxgene-schema-apply)
  - [`cellxgene schema validate`](#cellxgene-schema-validate)
  - [Local testing and dropbox upload](#testing-locally-optional-and-upload-to-dropbox)

<br/>

---

## Quick start

<br/>

If you are already familiar with cellxgene, anndata, and the cellxgene data schema, then give curation a shot with this [quick start guide](https://github.com/chanzuckerberg/single-cell-curation#quick-start). Otherwise, continue reading!

---

<!--- ![image](https://user-images.githubusercontent.com/25663501/111403110-d1589700-86a2-11eb-99dd-bf52e348d43a.png) --->

<!---
<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/111403110-d1589700-86a2-11eb-99dd-bf52e348d43a.png" width="1000">
</p>

<p align="center">
  <b> Figure: </b> Collection View of the Cellxgene Data Portal
</p>

<br/>


Reason to contribute:
 - link permanence
 - sharing of private datasets with collaborators
 - no barrier for readers to explore your dataset
 - accesibility of your dataset in the major single cell data formats

<br/>
 
--->

## Required data and `anndata` structure

<p align="center">
  <img width="500" src="https://user-images.githubusercontent.com/25663501/111377611-3c8c7400-8677-11eb-8176-cf9de8f64c70.png">
</p>

<p align="center">
  <b> Figure:</b> anndata components
</p>

<br/>

The following components are required for submission to the cellxgene data portal. 

- raw count matrix (except for certain assays such as scATACseq)
- normalized expression matrix used for visulization
- cell level metadata (barcodes, cell type, tissue of origin, etc.)
- embedding (at least one required, UMAP, tSNE, spatial, PCA)

Additionally, variable and feature level metadata can be useful to include but is not required for construction of an `anndata` object.

These components should be stored in the following locations in an `anndata` object (for more information, please refer to the [anndata documentation site](https://anndata.readthedocs.io/en/latest/)):

<br/>

<div align="center">

| Component                       | `anndata` location                     | Notes                                                                  |
| ------------------------------- |:--------------------------------------:|:-----------------------------------------------------------------------|
| raw count matrix                | `adata.raw.X` or `adata.layers['raw']` | Necessary, with some exceptions (see [exceptions](#alternative-assays) |
| normalized expression matrix    | `adata.X`                              | Used for visualization in cellxgene explorer                           |
| cell level metadata             | `adata.obs`                            |                                                                        |
| variable/feature level metadata | `adata.var`                            |                                                                        |
| embedding                       | `adata.obsm`                           | Must start with the prefix 'X_' (i.e. adata.obsm['X_UMAP'])            |

</div>

<br/>

<div align="center">
  <b> Table: </b> Required data and `anndata` object structure
</div>

<br/>

<!---
**Note:** Some assays, such as scATACseq or other epigenetic assays may not have a standardized way of representing a raw count matrix. While we still require there to be a matrix in location annotated as `raw`, some alternative options exist for these assays such as putting the un-normalized gene activity matrix in the `raw` slot. Since `adata.raw.X` can take matrices that are of a different dimensionality than `adata.X`, you could potentially put in a peak x cell matrix in this location instead. We leave this up to the author's discretion, but are happy to chat about options.
--->

**Note:** In addition to these data, other representations of the expression matrix (alternative normalizations, SCTransform, corrected counts from SCTransform or background corrected counts) can all be stored as `layers` in your anndata object (as long as they maintain the same dimensionality of the main expression matrix used for visualization).

**Note:** Information which pertains to the cellxgene schema will be stored in the `adata.uns` and `adata.obs` slots of the `anndata` object and will be discussed in the [next section](#schema-definition).
 
### Format conversion
 
There are a handful of tools that can be used to convert different single cell formats (`seurat`, `loom`, `sce`, `cds`) to `anndata`. We have had the most success with the R package [sceasy](https://github.com/cellgeni/sceasy). For relatively straightforward conversions, it is suitable to use the `sceasy::convertFormat()` function as specified in the sceasy documentation. For scenarios where you have multiple assays in the same object (as you may encounter in the `seurat` toolchain) it is recommended to check out the function definitions of the methods that are called by `convertFormat()`. For instance, if you look at the source code for `sceasy::seurat2anndata()`, you will see that the extra parameters `assay`, `main_layer`, and `transfer_layers` can be specified for control over which elements of the `seurat` object are carried into the `anndata` object. Check out the following function definitions and potentially use as a starting point for more customized conversions:

- [`seurat2anndata()`](https://github.com/cellgeni/sceasy/blob/f8f0628a280e0880ea94b00100b463e1f6ba1994/R/functions.R#L14)
- [`sce2anndata()`](https://github.com/cellgeni/sceasy/blob/f8f0628a280e0880ea94b00100b463e1f6ba1994/R/functions.R#L64)
- [`loom2anndata()`](https://github.com/cellgeni/sceasy/blob/f8f0628a280e0880ea94b00100b463e1f6ba1994/R/functions.R#L116)


**Note:** While `anndata` is able to accomodate multiple representations of an expression matrix in the object, matrices that are stored in `adata.layers` are required to be the same dimensions as `adata.X` (`adata.raw.X` may be of a different dimensionality though). In some scenarios, you may need to construct different `anndata` objects to accomodate different `assays` in the same experiment (for example, spliced vs unspliced counts in a sNuc-seq experiment).

<br/>

### Alternative assays

The cellxgene data portal and explorer are able to handle a wide range of single cell data types. Due to requirements of the cellxgene schema and limitations of the `anndata` object structure, there are some assay specific considerations that are outlined below

- ATAC/mC
  - raw data: Since there is not a standard way to generate a counts matrix or gene activity matrix for single cell epigenetic data, it is suitable to put a copy of `adata.X` in place of a raw data matrix. You can specify this assignment in `adata.uns`
  - [Link to example curated dataset](www.example.com)
- RNAseq
  -  In some single cell toolchains, outputs of different computational methods may be of a different dimensionality than the input matrix (i.e SCTransform, or the scale.data slot in Seurat objects). Since `anndata` objects do not allow for matrices of different dimensions in `adata.layers`, it is suitable to pad the missing rows/features with zeros to ensure that indices (feature names) across different the different layers are matched. You can specify that these matrices were padded in `adata.uns['layer_descriptions']` as a part of the cellxgene schema. This allows users who may reuse your data to remove the padding from the matrix before working with the data further.
  -  [Link to example curated dataset](www.example.com)
- CITEseq
  - TODO: considerations for CITE-seq...
  - [Link to example curated dataset](www.example.com)
- Spatial/Visium
  - Since spatial sequencing spots are not at the single cell level, you may wish to include prediction matrices that provide deconvolution information for each spot. These can generally be stored in `adata.uns` although it can be interesting to add prediction scores for each spot as new columns to the `adata.obs` dataframe
  - [Link to example curated dataset](www.example.com)

<br/>

---

## Schema definition

<br/>

The purpose of the cellxgene schema is to support the construction of a data corpus that facilitates data integration across multiple tissues and experiments. This goal requires that we collect a standardized set of metadata about single cell datasets that are uploaded to the cellxgene data portal. To make this process easy to adhere to, we only require a few fields (detailed below) to support easy search and integration across datasets. These metadata are stored within the `anndata` object (in `adata.uns` and `adata.obs`). To access a more comprehensive decription about our schema requirements, you can refer to the [official schema definition](corpora_schema.md).

In this section, we are covering 1) what metadata are required to adhere to the cellxgene schema and 2) the locations of these metadata in the `anndata` object. While it is possible to build an object that is acceptable by the cellxgene schema manually, in the [next section](#cellxgene-curation-tools), we will show you how to perform this curation with tools provided by the CZI curation team.

#### Basic requirements ([expanded version](corpora_schema.md#basic-requirements))
- **Unique observation identifiers:** Each observation (usually a cell) must have an id that is unique within the dataset.
- **Unique feature identifiers:** Every feature (usually a gene or transcript) also needs a unique identifier.
- **No PII:** No metadata can be personally identifiable: no names, dates of birth, specific locations, etc. There's a [list](https://docs.google.com/document/d/1nlUuRiQ8Awo_QsiVFi8OWdhGONjeaXE32z5GblVCfcI/edit?usp=sharing).

#### `uns`

`adata.uns` is a dictionary of key values pairs that we use to store dataset level metadata. This is also where you store information on the different data representations that you are sharing (i.e you are sharing a raw counts matrix, normalized expression matrix, scaled and centered expression matrix)

<br/>

<div align="center">

| Key                         | Value      | Description                                                                  | Example                        |
| --------------------------- |:----------:| :----------------------------------------------------------------------------| -------------------------------|
| 'version'                   | string     | schema version (current version is `1.1.0`)                                  | `1.1.0`                        |
| 'title'                     | string     | title of publication (and title of dataset if more than one is being submitted) | `An Atlas of Gene Regulatory Elements in Adult Mouse Cerebrum: GABAergic neurons`|
| 'publication_doi'           | string     | DOI of preprint or official publication| `https://doi.org/10.1101/2020.05.10.087585` |
| 'organism'                  | string     | name of organism (first letter capitalized) | `Mouse` |
| 'organism_ontology_term_id' | string     | NCBI Taxon ID| `NCBI:txid10090` |
| 'layer_descriptions'        | dictionary | a set of key value pairs defining the locations of different representations in the `anndata` object and an description of that representation. One of these key-value pairs must be on of `raw: raw` or `raw.X: raw`    | `{'X': 'log1p' (this value can be free text)}` |

</div>

<div align="center">
  <b> Table: </b> Required and suggested dataset level metadata
</div>

<br/>

In the above table we see what type on information is necessary at the dataset level. Two fields to pay extra attention to are the `title` and the `layer_descriptions` keys. Specifically, `title` should contain the name of the publication that the dataset is coming from. If there is more than one dataset associated with the publication, then the name of the dataset should be appended to the title (i.e. `'title': 'publication name: dataset title'`). This field will be used to name and distinguish different datasets in the same collection in the portal.

At least one of the key value pairs in `layer_descriptions` needs to indicate the presence of a raw counts layers. This can either be specified as `raw: raw` or `raw.X: raw`. Additional layers may be specified and the values for these keys may be as descriptive as necessary (i.e. `scale.data: scaled and centered data`) 

<br/>
<!---
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
  <b> Table: </b> Required and suggested data layers
</div>

<br/>
--->

#### `obs`

`adata.obs` is a dataframe that is used to store cell level metadata. In addition to experiemental metadata, the following additional field are required:

| `obs` column name                    | Type   | Description                                                       | Example                                   |
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
  <b> Table: </b> Required cell level metadata
</div>

<br/>

In the above table, the types of cell level metadata that is collected for the schema. In particular, it is imporant to exapnd on how to annotate `cell_type` field. In general, we try to find most specific ontology term that accurately represents a given cell type in your dataset. In the event that your cell type is not described accurately by any term in the `CL` ontology, then it is sufficient to report the most accurate `UBERON` term and id instead. In some cases, when neither the `CL` nor the `UBERON` term are descriptive enough, then it is possible to submit an empty string (`''`) for that entry.

<br/>

---

## Cellxgene curation tools

<br/>

The cellxgene curation tools include functions that can make the curation process easier for you. Essentially the workflow looks like this:
1. create `config.yaml` that will specify how the schema is to be applied to your anndata object (see example [here](example_config.yaml))
2. run `cellxgene schema apply`, which takes your config and source anndata file to produce a curated anndata object.
3. use `cellxgene schema validate` to ensure that your curated object meets the schema requirements

You can check out the full documenatation for this tooling [here](#schema_guide.md), but a quick introduction to the tools is below:

<br/>

### Installation

`pip install cellxgene-schema`

<br/>

### `config.yaml`

<br/>

Your [`config.yaml`](example_config.yaml) file is used to update values and columns in `adata.uns` and `adata.obs` (repectively) with required schema information. Here is an example/mock of a config file:

#### `uns` section

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

In the config file snippet above, our 0 level indendentation specifies that we are modifying the `uns` slot (remember that `adata.uns` is a dictionary of key-value pairs). The next level of indentation specifies a key name to add to `uns`. The key's corresponding value will the be string following the colon, or if the key is follow by more lines that are further indented, then the corresponding value will be a dictonary containing the key-value pairs specified in the subsequent lines. More concretely, `adata.uns['layer_descriptions']`, will return a dictionary with the key value pairs `{'raw.X': 'raw', 'X': 'log1p'}` after the schema config has been applied.


**Note:** at least one of the keys in `layer_descriptions` must return the value 'raw'

<br/>

#### `obs` section
```
obs:
    tissue_ontology_term_id: UBERON:0000006              #UBERON tissue term
    assay_ontology_term_id: EFO:0010961                  #EFO assay term
    disease_ontology_term_id: PATO:0000461               #MONDO disease term or PATO:0000461
    sex: male                                            #male, female, mixed, unknown, or other
    ethnicity_ontology_term_id: na                       #HANCESTRO term
    development_stage_ontology_term_id: HsapDv:0000160   #HsapDv term
    cell_type_ontology_term_id:
      cell_label:                                        #this is column name in your obs dataframe (the field that specifies cell type)
        acinar: CL:0002064                               #mapping between cell types (specified in anndata object) and cl ontology id
        alpha: CL:0000171
        delta: CL:0000173
        endothelial: CL:0000115
        epsilon: CL:0005019

```


<br/>

In the config file snippet above, the 0 level indentation specifies that we are modifying `obs`. At the first level of indentation, we start adding schema fields to `obs`. If the first level indentation fields ar follow by an ontology value (i.e. `tissue_ontology_term_id: UBERON:0000006`), then a column called `tissue_ontology_term_id` will created in the obs dataframe, and all of its values will be `UBERON:0000006`. Additionally, another column named `tissue` will be created with the corresponding UBERON term (`islet of Langerhans`).

On other hand, if the schema field is followed by lines with have further indentation, then the next line in the config file specify a column whose values should be mapped to ontology terms (with these mappings specified by the subsequent lines). Two new columns will be added to the `obs` data frame `cell_type_ontology_term_id` and `cell_type`.

<br/>

#### `fixup_gene_symbols` section (gene symbol harmonization)

The last section describes how gene symbol conversion should be applied to each of the layers. This is similar to the
`layer_descriptions` field above, but there are only three permitted values: `raw`, `log1p`, and `sqrt`:

```
fixup_gene_symbols:
  X: log1p
  raw.X: raw
```

This tells the script how each each layer was transformed from raw values that can be directly summed. `raw` means that
the layer contains raw counts or some linear tranformation of raw counts. `log1p` means that the layer has `log(X + 1)`
for each the raw `X` values. `sqrt` means `sqrt(X)` (this is not common). For layers produced by Seurat's normalization
or SCTransform functions, the correct choice is usually `log1p`.

**Note**: If a layer is not specified in the `fixup_gene_symbols` section, then it will not be carried over into the curated object. If no layers are specified in `fixup_gene_symbols` or if `fixup_gene_symbols` is not included the `config.yaml` file at all, then all layers will be carried over into the curated object.

<br/>


### `cellxgene schema apply`
 
<br/>

In order to use the `cellxgene schema apply` command, you will need to pass the following required arguments:
- `--source-h5ad` your original anndata file
- `--remix-config` the `config.yaml` file that we specified above 
- `--output-filename` the name of the resulting anndata that is consistent with the cellxgene schema

<br/>

The next step will be to validate the resulting object.

<br/>

### `cellxgene schema validate`

<br/>

In order to validate the remixed object, needs to simply run `cellxgene schema validate remixed_anndata.h5ad`. If there has been no terminal output from the function, then your object has been validated successfuly and is ready for upload!

<br/>

### Testing locally (optional) and upload to dropbox

You can test you `anndata` object after schema application by running with a local installation of cellxgene (`cellxgene launch example.h5ad`). This allows you to preview how your dataset will appear in the cellxgene explorer view within the data portal. Following this optional testing, you need upload to dropbox which is required since datasets cannot be uploaded directly to the portal.

---

## Uploading data to the cellxgene data portal

In general, the cellxgene data portal is oriented around grouping datasets by their association with a particular publication. We refer to all of the datasets associated with a publication as a collection. When you go to the homepage of the [cellxgene data portal](#https://cellxgene.cziscience.com/), you will see a list of collections (publications). To explore any of the data associated with these publications, you would click on a particular collection and be greeted by all of the datasets associated with that publication. In this next section, we will cover how to contribute your newly curated data to the cellxgene portal. Briefly, we will cover:

- Login options
- Discovering 'My Collections'
- Creating a **private** collection
- Adding datasets to a newly created collection
- Sharing **private** collections with collaborators via private links
- Publishing your collection for public access

### Portal Sign in

sign in options for the portal (gmail, github, etc)

<!--- ![image](https://user-images.githubusercontent.com/25663501/113528770-9f46a080-958f-11eb-83a1-36620e1543d2.png) --->

<p align="center">
  <img width="250" src="https://user-images.githubusercontent.com/25663501/113528770-9f46a080-958f-11eb-83a1-36620e1543d2.png">
</p>

<p align="center">
  <b> Figure:</b> Login options
</p>

<br/>

The cellxgene data portal offers multiple login options including sign-in with your gsuite and github accounts, or via email address.

<br/>


<!--- ![image](https://user-images.githubusercontent.com/25663501/113530303-a96a9e00-9593-11eb-91b7-0898d521314a.png) --->

<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/113530303-a96a9e00-9593-11eb-91b7-0898d521314a.png">
</p>

<p align="center">
  <b> Figure:</b> Navigate to "My Collections"
</p>

<br/>

After you have logged in,  you will be able to navigate to the 'My Collections' page where you can view all of the datasets that you have uploaded to the portal.



<!--- ![image](https://user-images.githubusercontent.com/25663501/113528943-18de8e80-9590-11eb-92ab-2a606879b70c.png) --->

<!--- ![image](https://user-images.githubusercontent.com/25663501/113530607-6d840880-9594-11eb-8af7-17f75d3a8a40.png) --->


<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/113530607-6d840880-9594-11eb-8af7-17f75d3a8a40.png">
</p>

<p align="center">
  <b> Figure:</b> My Collections
</p>

<br/>

In the image above, we see an example of how the 'My Collections' page will look after you have uploaded a few datasets to the dataportal. Note that some of the collections in this example are published (publicly available) and one is private. The private collection is only available to people who you share the private link with. This allows for control over who views your datasets and for making amendmendents to datasets that have been uploaded to the portal previously. To add a new collection, we can simply click on the button highlighted above. 


### Create a collection

<!--- ![image](https://user-images.githubusercontent.com/25663501/113528984-36135d00-9590-11eb-82a5-68ae4e5ddb1b.png) --->

<p align="center">
  <img width="500" src="https://user-images.githubusercontent.com/25663501/113528984-36135d00-9590-11eb-82a5-68ae4e5ddb1b.png">
</p>

<p align="center">
  <b> Figure:</b> Creating a collection
</p>

<br/>

Once you have clicked the 'Create Collection' button, you will be prompted to enter the following information:

- Collection Name (name of the publication)
- Description (abstract of publication)
- Contact Name (PI name) and contact email (PI email)
- Add Link (link to relevant resources such as the publication itself, lab website, consortia/department website, entry in GEO or other raw data sources, etc.)

### Add a dataset to a collection

Highlight add dataset button

<!--- ![image](https://user-images.githubusercontent.com/25663501/113529554-b5556080-9591-11eb-8e0e-274fcb1ff39a.png) --->

<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/113529554-b5556080-9591-11eb-8e0e-274fcb1ff39a.png">
</p>

<p align="center">
  <b> Figure:</b> Add dataset to a collection
</p>

<br/>

Once you have created a collection, you will have the ability to add a dataset. Once you click this button, you will be guided to dropbox, where you can select a dataset for upload.


Choose from Dropbox:

<!--- ![image](https://user-images.githubusercontent.com/25663501/113529612-db7b0080-9591-11eb-9720-ae47155d62f8.png) --->

<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/113529612-db7b0080-9591-11eb-9720-ae47155d62f8.png">
</p>

<p align="center">
  <b> Figure:</b> Choose dataset from dropbox
</p>

<br/>


### Remove dataset from a collection

`screenshot removing dataset`

<!--- ![image](https://user-images.githubusercontent.com/25663501/113529632-edf53a00-9591-11eb-8ea4-575ca8c643f7.png) --->

<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/113529632-edf53a00-9591-11eb-8ea4-575ca8c643f7.png">
</p>

<p align="center">
  <b> Figure:</b> Remove dataset from a collection
</p>

<br/>

In case you have uploaded the wrong dataset or want to update a dataset within a collection, you can delete a dataset associated with a collection (without deleting the collection itself).

### Share uploaded datasets with private links
Upon successful dataset upload....


### Publish Collection to the portal

`Screenshot of Explorer mode`
