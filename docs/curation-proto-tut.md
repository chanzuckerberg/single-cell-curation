# Cellxgene Curation Tutorial

<img align="right" src="https://user-images.githubusercontent.com/25663501/111377133-b07a4c80-8676-11eb-8eb8-07ca4d7a77e9.png">

<!--- ## Overview --->

[Cellxgene's publishing platform](https://cellxgene.cziscience.com/) and [interactive single cell data explorer](https://github.com/chanzuckerberg/cellxgene) are optimized for access, exploration and reuse of single cell data. To enable this, all data published on cellxgene adhere to a standardized data schema, which allows for efficient navigation and integration of the growing number of single cell datasets. Publishing your data on the cellxgene portal offers the following benefits:

 - Link permanence (you can reference in your publication without ever worrying about dead links)
 - Sharing of private datasets with collaborators (keep the data private until it is ready for publication)
 - Enables visual exploration that works at scale (and no need for you to build your own single cell data explorer)
 - Availability of your dataset in major single cell data formats (including `AnnData` and `seurat`)

To publish your data on the platform, cellxgene requires that you create [AnnData](https://anndata.readthedocs.io/en/latest/#) objects adhering to a [succinct data schema](corpora_schema.md). 
   
This tutorial will consist of an explanation and demonstration of: 
 - how to create and structure an `AnnData` object with your single cell data
 - how to augment the `AnnData` object with metadata required by the cellxgene schema (this information will be used to make your data findable, reusable, and interoperable)
 - how to upload this object to the cellxgene data portal. 
   
If you run into any issues during this tutorial, or have any suggestions on how to improve the portal and curation experience, you can contact us via [cellxgene@chanzuckerberg.com](cellxgene@chanzuckerberg.com).

### Table of Contents

- [Quick start](#quick-start)
- [Data requirements and AnnData structure](#required-data-and-anndata-structure)
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
- [Uploading data to the cellxgene data portal](#uploading-data-to-the-cellxgene-data-portal)
  - [Portal sign in](#portal-sign-in)
  - [Create a private collection](#create-a-private-collection)
  - [Add a dataset](#add-a-dataset-to-a-private-collection)
  - [Remove a dataset](#remove-dataset-from-a-private-collection)
  - [Dataset actions](#dataset-actions)
  - [Publish collection](publish-collection-to-the-portal)
- [Gene sets (under construction)](#gene-sets)

<br/>

---

## Quick start

<br/>

If you are already familiar with cellxgene, AnnData, and the cellxgene data schema, then give curation a shot with this [quick start guide](https://github.com/chanzuckerberg/single-cell-curation#quick-start). Otherwise, continue reading!

---

## Use Case

The cellxgene curation process consists of only a few steps:

- Transform your data into AnnData format.
- Apply schema to your single cell data (curation).
- Validate that the curation process has succeeded.
- Upload curated data to the cellxgene data portal.

## Creating a curation environment (optional)

To start this tutorial, we will be creating a conda environment where we can install the cellxgene curation tools and their dependencies. If you already have conda installed ([reference for the uninitiated](https://conda.io/projects/conda/en/latest/index.html)), this is a simple process and you can run the following in a new terminal:

```bash
conda create -n "curation_env" python=3.8.0
conda activate curation_env
```
After you run the `activate` command, you should be inside of your conda environment. We can install some required packages by running the appropriate pip commands inside of our activated environment:

```bash
pip install pandas
pip install scanpy
```

## `AnnData` Preparation

The main prequisite for submitting your data to the cellxgene data portal is that it is structured in `AnnData` format. Below we see the rough steps for downloading our demo dataset: PBMC 3K from 10X. You can read more about the origin of this dataset from the [10X website](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k). We will start preparing the `AnnData` object by downloading the necessary files below (note: you can also access and download the necessary files from [this dropbox link](https://www.dropbox.com/sh/pgby9xn0v7tx3vi/AAAB8dTsXcjjFgNiIOGw7yYNa?dl=0)). 

```bash
curl https://cellxgene-curation-tutorial.s3.us-east-2.amazonaws.com/counts.mtx >> counts.mtx  # count matrix
curl https://cellxgene-curation-tutorial.s3.us-east-2.amazonaws.com/normalized_expression.mtx >> normalized_expression.mtx  # x matrix
curl https://cellxgene-curation-tutorial.s3.us-east-2.amazonaws.com/cell_metadata.tsv >> cell_metadata.tsv  # cell metadata
curl https://cellxgene-curation-tutorial.s3.us-east-2.amazonaws.com/umap_embedding.csv >> umap_embedding.csv  # embedding
curl https://cellxgene-curation-tutorial.s3.us-east-2.amazonaws.com/var_metadata.tsv >> var_metadata.tsv  # var metadata
```

Next, we can load each of these components into a python environment with ScanPy and Pandas installed ([ScanPy docs](https://anndata.readthedocs.io/en/latest/) and [Pandas docs](https://pandas.pydata.org/docs/getting_started/index.html)). 

Enter a python environment by typing the following into your terminal:

```bash
python
```

Once we see our python prompt (indicated by `>>>`), we can read in our different components and construct our Anndata object like so:

```python
import pandas as pd
import scanpy as sc

# Reading in data - you can also check out more scanpy functions for reading in different file formats...
raw_adata = sc.read_mtx('counts.mtx')  # Read in counts matrix - returns an AnnData object - count matrix is stored in raw_adata.X
adata = sc.read_mtx('normalized_expression.mtx')  # Read in normalized expression matrix - normalized expression matrix is stored in adata.X
embeddings = pd.read_csv('umap_embedding.csv', header = None).to_numpy()  # Read in embeddings and convert into numpy array; AnnData requires obsm objects be formatted as numpy arrays
obs_metadata = pd.read_csv('cell_metadata.tsv', sep='\t')  # Read in cell metadata (cluster assignment, qc metadata, etc...)
var_metadata = pd.read_csv('var_metadata.tsv', sep='\t')  # Read in variable metadata (dispersion, function annotation, etc)

# Assign components to appropriate locations in AnnData
adata.raw = raw_adata  # stash our counts/raw adata in the raw slot of the normalized adata object - counts will be accessible via adata.raw.X
adata.obsm['X_umap'] = embeddings  # embedding name must be prefixed with 'X_'
adata.obs = obs_metadata
adata.var = var_metadata
```

There may be scenarios where you want to perform conversion of the values in your `obs` data frame to change them into a representation that is human readable. In our example, our observational metadata (`adata.obs`) is not quite amenable to the curation process (specifically, the `leiden`). This is because our cluster assignments are currently integers and are not biologically interpretable. To fix this, we can update the cluster metadata in our `AnnData` object to map the cluster IDs to the cell types they label:

```python
# Create a list to map current cluster ids (stored in adata.obs['leiden']) to cell type names
# these map to the ordering of the current cluster ids in the object, (which in this case is 0, 1, 2, 3, ...)

# First convert integer values in our leiden column to categorical values instead of integer values
adata.obs['leiden'] = adata.obs.leiden.astype('category')

# Define the appropriate cell names/types to map the column to
cell_names = ['CD4 T', 'CD14 Monocytes', 'B', 'CD8 T', 'NK', 'FCGR3A Monocytes', 'Dendritic', 'Megakaryocytes']

# Perform the mapping
adata.rename_categories('leiden', cell_names)

# You can check out the resulting values in your obs dataframe like so
adata.obs.head()
```

Above we have kept the column `leiden`, but renamed all the categories in that column to labels that we will be able to link to an ontology. After this step, we are ready to write our `AnnData` file and use the cellxgene curation command line tools. Please note that the dataset also has an identical column called louvain where this step had already been performed previously. You can use either of these columns to perform the demo curation process. Going forward, we will continue using the `leiden` column for annotating cell types in the AnnData object. 

```python
adata.write_h5ad('pbmc3k.h5ad')
```

Now that we have created and written our AnnData file, we are ready to pursue the curation process. It is important to note that the curatin process is additive and does not change any of the information that was stored in the original AnnData object. To get a better sense of what the AnnData object should look like after schema injection you can refer to [this encoding document](https://github.com/chanzuckerberg/single-cell-curation/blob/ambrosecarr/schema-v1.1.1/schema/1.1.1/anndata_encoding.md)

## Cellxgene curation tools

<br/>

The cellxgene curation tools provide functions that can make the curation process easier. They allow you to specify the appopriate ontologies that capture important metadata about your object and do so in an additive way (i.e. no original metadata in your dataset are removed during schema application). The workflow to apply the schema looks like this:
1. Create a `config.yaml` file that will specify how the schema is to be applied to your `AnnData` object (see a worked example matching the pbmc3k dataset [here](example_config.yaml)).
2. Run `cellxgene schema apply`, which takes the `config.yaml` and source `AnnData` file to produce a curated `Anndata` object.
3. Use `cellxgene schema validate` to ensure that your curated object meets the schema requirements.

You can check out the full documentation for this tooling [here](#schema_guide.md), but a quick introduction to the tools is provided below.

<br/>

### Installation

You can install the cellxgene curation tools either inside of a conda environment/virtual environment (recommended) or in your system python install like so: 

```bash
pip install cellxgene-schema
```

<br/>

### `config.yaml`

<br/>

The [`config.yaml`](example_config.yaml) file is used to update values and columns in `adata.uns` and `adata.obs` (respectively) with required schema information. All original columns in the dataset will be preserved. For fields where you need to enter an ontology ID, you can use the [EBI Ontology Lookup Service](https://www.ebi.ac.uk/ols/index) to search for the appropriate terms.  This file will also standardize gene symbols in your dataset via the `fixup_gene_symbols` section. You can download the example config to your local machine (and redirect it to a local file named `config.yaml`) using the following command (and view/edit using any text editor):

```bash
curl https://raw.githubusercontent.com/chanzuckerberg/single-cell-curation/tutorial-prototype/docs/example_config.yaml >> config.yaml
```

The below sections walk through the downloaded config file which matches the pbmc3k dataset.

#### `uns` section

```yaml
uns:
    version:
        corpora_schema_version: 1.1.0                          # ex: find current schema version number in the corpora schema definition
        corpora_encoding_version: 0.1.1                        # ex: encoding version
    title: 10X PBMC Demo                                       # free text field
    publication_doi: https://doi.org/00.0000/2021.01.01.000000
    layer_descriptions:
        raw: raw                                               # it is essential for at least one the layer_descriptions to be set to 'raw'
        X: log1p                                               # free text description of normalization method
    organism: Human                                            # Specify organism name
    organism_ontology_term_id: NCBITaxon:9606                  # Specfiy organism NCBI Taxon ID

```
<br/>

In the config file snippet above, the 0 level indentation specifies that the `uns` slot is being modified (remember that `adata.uns` is a dictionary of key-value pairs). The next level of indentation specifies a key name to add to `uns`. The key's corresponding value is the string following the colon, or if the key is followed by more lines that are further indented, then the corresponding value is a dictionary containing the key-value pairs specified in the subsequent lines. More concretely, `adata.uns['layer_descriptions']` is added to `anndata.uns` as a dictionary with the key value pairs `{'raw': 'raw', 'X': 'log1p'}` when `cellxgene-schema apply` is run.

**Note:** [at least one of the keys in `layer_descriptions` must return the value 'raw'](https://github.com/chanzuckerberg/single-cell-curation/blob/tutorial-prototype/docs/schema_guide.md#unstructured-metadata)

<br/>

#### `obs` section

The cell metadata required by the cellxgene schema are structured metadata that are part of ontologies that enable powerful search and comparison capabilities. For this reason, Before describing the `obs` section of the `config` file, it's important to discuss how to look up ontologies for a particular field (i.e. cluster/cell type). In general we aim to find the highest resolution term that accurately describes a field entry. To give a more concrete example, we can consider the `louvain` column in our observational metadata. The values in this column correspond to cell type annotations. We can view the different cell types that are present in the dataset by running:

```python
>>> import scanpy as sc

>>> adata = sc.read("pbmc3k.h5ad")
>>> set(adata.obs['louvain'])
{'NK', 'Dendritic', 'CD4 T', 'FCGR3A Monocytes', 'Megakaryocytes', 'B', 'CD8 T', 'CD14 Monocytes'}
```

To find an ontology that maps to a given cell type (for this example, consider the Natural Killer cells (coded as `NK` in the `AnnData` object)), head to [EBI's Ontology Lookup Service](https://www.ebi.ac.uk/ols/index), search for 'Natural Killer Cells', filtered by the desired ontology (in this case, `CL` for cell type). A screenshot of such a lookup is presented below.

<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/118162747-1cff9680-b3ef-11eb-92ea-11defb7c88b4.png">
</p>

<p align="center">
  <b> Figure:</b> EBI Ontology Lookup Service: https://www.ebi.ac.uk/ols/index
</p>

In this image, the search presents a list of ontology terms and ids that could potentially fit "Natural Killer cells". Choose the term that most closely represents the cell type query. If the data contains very specific cell types, they may not be represented by equally specific ontology terms. In this case, cellxgene recommends using a higher level term or collapsing multiple fine grained cell types into a more general term (This might occur in brain datasets where neuronal subtypes may be denoted by combinations of specific markers and neurotransmitters).

Once this lookup has been performed for all of the cell type present in the `louvain` column, the appropriate fields can be filled in our config file. Below is a skeleton config file with the cell type ontologies filled in:

```yaml
obs:
    tissue_ontology_term_id:                             #UBERON tissue term
    assay_ontology_term_id:                              #EFO assay term
    disease_ontology_term_id:                            #MONDO disease term or PATO:0000461
    sex: unknown                                         #male, female, mixed, unknown, or other
    ethnicity_ontology_term_id: unknown                  #HANCESTRO term
    development_stage_ontology_term_id:                  #HsapDv term
    cell_type_ontology_term_id:
      louvain:                                           #this is column name in your obs dataframe (the field that specifies cell type)
        NK: CL:0000623                                   #Link an ontology ID to each value in your cell type column
        Dendritic: CL:0000451
        CD4 T: CL:0000624
        FCGR3A Monocytes: CL:0002396
        Megakaryocytes: CL:0000556
        B: CL:0000236
        CD8 T: CL:0000625
        CD14 Monocytes: CL:0001055
```

Essentially, the mapping is performed for the `cell_type_ontology_term_id` column by specifiying 1) the target column (second level of indentation) and  2) each of the unique values in the target along with their corresponding ontology term ids (third level of indentation). Now that we understand how to map ontology term ids for a column like cell type, let's find out how to fully specify the `obs` section of the config file, like in the example below:

```yaml
obs:
    tissue_ontology_term_id: UBERON:0000178              #UBERON tissue term
    assay_ontology_term_id: EFO:0008913                  #EFO assay term
    disease_ontology_term_id: PATO:0000461               #MONDO disease term or PATO:0000461
    sex: unknown                                         #male, female, mixed, unknown, or other
    ethnicity_ontology_term_id: unknown                  #HANCESTRO term
    development_stage_ontology_term_id: HsapDv:0000087   #HsapDv term
    cell_type_ontology_term_id:
      louvain:                                           #this is column name in your obs dataframe (the field that specifies cell type)
        NK: CL:0000623                                   #Link an ontology ID to each value in your cell type column
        Dendritic: CL:0000451
        CD4 T: CL:0000624
        FCGR3A Monocytes: CL:0002396
        Megakaryocytes: CL:0000556
        B: CL:0000236
        CD8 T: CL:0000625
        CD14 Monocytes: CL:0001055
```

<br/>

In the config file snippet above, the 0 level indentation specifies that we are modifying `obs`. At the first level of indentation, we start adding schema fields to `obs`. If the first level indentation field is followed by an ontology term id (i.e. `tissue_ontology_term_id: UBERON:0000178`), then a column for that field will be created in the obs dataframe (i.e. a column called `tissue_ontology_term_id` will be created), and all of its values will take on the value specified after the colon (i.e. `UBERON:0000178`). Additionally, the cellxgene curation tools will search the ontology term id and create a corresponding column with the appropriately corresponding term (i.e. a column called `tissue` will be created and all values of that column will be `blood`).

On other hand, if the schema field is followed by lines which have further indentation, then the next line in the config file will specify a column whose values should be mapped to ontology terms (with these mappings specified by the subsequent lines). This is typically more relevant for specifying cell labels. In our example above, two new columns will be added to the `obs` data frame `cell_type_ontology_term_id` and `cell_type`. The original column specifying cell identity will be retained (and will be tagged with `_original`; this is how cellxgene deals with column name collisions). The cellxgene curation tools never remove fields from the original dataset, they only add schema relevant information.

<br/>

#### `fixup_gene_symbols` section (gene symbol harmonization)

The cellxgene explorer requires unique gene symbols. The last section describes how gene symbol conversion should be applied to each of the layers, when necessary. This is similar to the
`layer_descriptions` field above, but there are only three permitted values: `raw`, `log1p`, and `sqrt`:

```yaml
fixup_gene_symbols:
  X: log1p
  raw.X: raw
```

This tells the script how each layer was transformed from raw values that can be directly summed. `raw` means that
the layer contains raw counts or some linear transformation of raw counts. `log1p` means that the layer has `log(X + 1)`
for each the raw `X` values. `sqrt` means `sqrt(X)` (this is not common). For layers produced by Seurat's normalization
or SCTransform functions, the correct choice is usually `log1p`. In our PBMC 3K example, the normalized expression matrix has been logged.

**Note**: If a layer is not specified in the `fixup_gene_symbols` section, then it will not be carried over into the curated object. If no layers are specified in `fixup_gene_symbols` or if `fixup_gene_symbols` is not included the `config.yaml` file at all, then all layers will be carried over into the curated object.

<br/>


### `cellxgene-schema apply`
 
<br/>

In order to use the `cellxgene schema apply` command, you will need to pass the following required arguments:
- `--source-h5ad` your original `AnnData` file
- `--remix-config` the `example_config.yaml` file that we specified above 
- `--output-filename` the name of the resulting `AnnData` that is consistent with the cellxgene schema

<br/>

For our pbmc3k dataset, the call to apply the schema will look like this:

```
cellxgene-schema apply --source-h5ad pbmc3k.h5ad --remix-config example_config.yaml --output-filename pbmc3k_curated.h5ad
```

It is important to note that this function is quite verbose. Even in a successful application of the schema, you will get warnings and output to the terminal:

```bash
(curation_env) user@MACOS curation-tutorial-test % cellxgene-schema apply --source-h5ad pbmc3k.h5ad --remix-config config.yaml --output-filename pbmc3k_curated.h5ad                    
/Users/user1/opt/anaconda3/envs/curation_env/lib/python3.8/site-packages/anndata/_core/anndata.py:120: ImplicitModificationWarning: Transforming to str index.
  warnings.warn("Transforming to str index.", ImplicitModificationWarning)
WARNING:root:Some symbols are simulaneously withdrawn and approved
We will treat them at approved:
{'LMO7-AS1', 'EXT2', 'HAP1', 'FHL1', 'APPL1', 'F11', 'MPP5'}
... storing 'assay_ontology_term_id' as categorical
... storing 'assay' as categorical
... storing 'ethnicity_ontology_term_id' as categorical
... storing 'ethnicity' as categorical
... storing 'sex' as categorical
... storing 'disease_ontology_term_id' as categorical
... storing 'disease' as categorical
... storing 'tissue_ontology_term_id' as categorical
... storing 'tissue' as categorical
... storing 'cell_type_ontology_term_id' as categorical
... storing 'cell_type' as categorical
... storing 'development_stage_ontology_term_id' as categorical
... storing 'development_stage' as categorical
```

In cases, where you have specified the config file incorrectly, or where the AnnData object has been constructed incorrectly, you will receive output from the function that can help you diagnose the issue. In most cases where the config file has been incorrectly specified, the function will still output a new Anndata object, however it will not pass the next validation step. In all cases where curation was unsuccessful, the `apply` subcommand will remove the schema version from the new object.

<br/>

### `cellxgene-schema validate`

<br/>

In order to validate the remixed object, needs to simply run `cellxgene schema validate path_to_remixed_anndata.h5ad`. If there has been no terminal output from the function, then your object has been validated successfully and is ready for upload!

For our pbmc3k dataset, the call to apply the schema will look like this:

```bash
cellxgene-schema validate pbmc3k_curated.h5ad
```
Contrary to the `apply` subcommand, `validate` is not very verbose and will print nothing in cases where the curation was successful:

```bash
(curation_env) user@MACOS curation-tutorial-test % cellxgene-schema validate pbmc3k_curated.h5ad
(curation_env) user@MACOS curation-tutorial-test % 
```
In cases where the curation was not a success, you will get a message that lets you that the schema version is missing. This means that you should return to the output of the `apply` subcommand to troubleshoot the issue.

<br/>

### Testing locally (optional) and upload to dropbox

You can test your curated `AnnData` object after schema application by running with a local installation of cellxgene (`cellxgene launch pbmc3k_curated.h5ad`). This allows you to preview how your dataset will appear in the cellxgene explorer view within the data portal. Following this optional testing, you need upload to dropbox which is required since datasets cannot be uploaded directly to the portal.

---

## Uploading data to the cellxgene data portal

In general, the cellxgene data portal is oriented around grouping datasets by their association with a particular publication. We refer to all of the datasets associated with a publication as a `collection`. When you go to the homepage of the [cellxgene data portal](#https://cellxgene.cziscience.com/), you will see a list of collections (publications). To explore any of the data associated with these publications, you would click on a particular collection and be greeted by all of the datasets associated with that publication. In this section, we will cover how to contribute your newly curated data to the cellxgene portal. Briefly, we will cover:

- Login options
- Discovering 'My Collections'
- Creating a **private** collection
- Adding datasets to a newly created collection
- Sharing **private** collections with collaborators via private links
- Publishing your collection for public access

### Portal Sign in

In order to upload your data go to [cellxgene.cziscience.com](cellxgene.cziscience.com), find the login section at the top right corner and select one of the sign in options for the portal (gmail, github, etc):

<p align="center">
  <img width="250" src="https://user-images.githubusercontent.com/25663501/113528770-9f46a080-958f-11eb-83a1-36620e1543d2.png">
</p>

<p align="center">
  <b> Figure:</b> Login options
</p>

<br/>

The cellxgene data portal offers multiple login options including sign-in with your gsuite and github accounts, or via email address.

<br/>

<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/113530303-a96a9e00-9593-11eb-91b7-0898d521314a.png">
</p>

<p align="center">
  <b> Figure:</b> Navigate to "My Collections"
</p>

<br/>

It should be noted that collections are displayed with specific metadata about the collections (including assays(s) used, organism(s) profiled, number of cells in the collection, the tissue profiled, and the disease state studied). After you have logged in,  you will be able to navigate to the 'My Collections' page where you can view all of the datasets that you have uploaded to the portal (if you are just starting out, the 'My Collections' view will be empty).

<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/113530607-6d840880-9594-11eb-8af7-17f75d3a8a40.png">
</p>

<p align="center">
  <b> Figure:</b> My Collections
</p>

<br/>

The image above shows an example of how the 'My Collections' page will look after a few datasets have been uploaded to the data portal. Note that some of the collections in this example are published (publicly available) and one is private. All collections start as private and are only made public when the author decides so. The private collection is only available to people who you share the private link with. This allows for control over who views your datasets and for making revisions to datasets that have been uploaded to the portal previously. To add a new collection, simply click on the button highlighted above.


### Create a private collection

<p align="center">
  <img width="500" src="https://user-images.githubusercontent.com/25663501/113528984-36135d00-9590-11eb-82a5-68ae4e5ddb1b.png">
</p>

<p align="center">
  <b> Figure:</b> Creating a collection
</p>

<br/>

Once the 'Create Collection' button has been clicked, you will be prompted to enter the following information:

- Collection Name (name of the publication)
- Description (abstract of publication)
- Contact Name (PI name) and contact email (PI email)
- Add Link (link to relevant resources such as the publication itself, lab website, consortia/department website, entry in GEO or other raw data sources, etc.)

You will also be required to agree to CZI's data submission policies which you can read in full by clicking 'Show Details'

### Add a dataset to a private collection

<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/113529554-b5556080-9591-11eb-8e0e-274fcb1ff39a.png">
</p>

<p align="center">
  <b> Figure:</b> Add dataset to a collection
</p>

<br/>

Once you have created a private collection, you will have the ability to add a dataset. Once you click this button, you will be guided to dropbox, where you can select a dataset for upload.

<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/113529612-db7b0080-9591-11eb-9720-ae47155d62f8.png">
</p>

<p align="center">
  <b> Figure:</b> Choose dataset from dropbox
</p>

<br/>


### Remove dataset from a private collection

<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/113529632-edf53a00-9591-11eb-8ea4-575ca8c643f7.png">
</p>

<p align="center">
  <b> Figure:</b> Remove dataset from a collection
</p>

<br/>

In case you have uploaded the wrong dataset or want to submit an updated dataset within a private collection, you can delete a dataset associated with a private collection (without deleting the private collection itself).

### Dataset actions

<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/113614805-f8aae000-9620-11eb-8421-359e7ec69f86.png">
</p>

<p align="center">
  <b> Figure:</b> Successful dataset upload
</p>

Once you have successfully uploaded a dataset to your private collection, you are given options to delete the dataset, download (in `anndata` (`.h5ad`) or `seurat V3` (`.rds`)), and explore the dataset using the cellxgene explorer. These options are given by the two icons on the right hand side of the entry.

### Publish Collection to the portal

<p align="center">
  <img src="https://user-images.githubusercontent.com/25663501/113615336-b930c380-9621-11eb-88ae-7855a49c5c61.png">
</p>

<p align="center">
  <b> Figure:</b> Publish your private collection
</p>

Finally, if your private collection is complete and you wish to share with the world, you simply need to hit publish (in the upper right hand side of the collection view). Congratulations! :partying_face:

## Gene sets

The cellxgene data portal will soon include the ability to upload relevant gene sets associated with the datasets in your collection. You can find a complete description of the requirements for the gene sets [here](./gene_sets.md)

In the future, you will be able to upload a file that matches the format described in the linked documentation to allow collaborators and external readers to see markers that you have deemed important for defining cell types and cell states.

## Future Work

Placeholder


## Appendices

<details>
  <summary>Appendix 1: AnnData Structure </summary>

## Required data and `AnnData` structure

<p align="center">
  <img width="500" src="https://user-images.githubusercontent.com/25663501/111377611-3c8c7400-8677-11eb-8176-cf9de8f64c70.png">
</p>

<p align="center">
  <b> Figure:</b> AnnData components
</p>

<br/>

The following components are required for submission to the cellxgene data portal:

- raw count matrix (except for certain assays such as scATACseq)
- normalized expression matrix used for visualization
- cell level metadata (barcodes, cell type, tissue of origin, etc.)
- embedding (at least one required, UMAP, tSNE, spatial, PCA)

Additionally, variable and feature level metadata can be useful to include but is not required for construction of an `AnnData` object.

These components should be stored in the following locations in an `AnnData` object (for more information, please refer to the [AnnData documentation site](https://anndata.readthedocs.io/en/latest/)):

<br/>

<div align="center">

| Component  | `AnnData` location            | Data Type                      | Notes                                                          |
| ---------- |:-----------------------------:|:------------------------------:|:---------------------------------------------------------------|
| raw count matrix | `adata.layers['raw']` | Numpy array or scipy sparse CSC matrix | Necessary, with some exceptions (see [exceptions](#alternative-assays))   |
| normalized expression matrix | `adata.X`  | Numpy array or scipy sparse CSC matrix | Used for visualization in cellxgene explorer                   |
| cell level metadata | `adata.obs` | Pandas dataframe | Categorical and continuous metadata shown in left and right cellxgene explorer sidebars respectively (can be used to color cells) |
| variable/feature level metadata | `adata.var`   | Pandas dataframe|                                                             |
| embedding                       | `adata.obsm`                           | Numpy array | Must start with the prefix 'X_' (i.e. adata.obsm['X_UMAP'])               |

</div>

<br/>

<div align="center">
  <b> Table: </b> Required data and `AnnData` object structure
</div>

<br/>

We can confirm this structure in our PBMC3K demo dataset. To get started try out the following:

```python
import scanpy as sc
adata = sc.read_h5ad('pbmc3k.h5ad')  # Read in object

adata.X.shape # dimensions of normalized expression matrix
```

**Note:** In addition to these data, other representations of the expression matrix (alternative normalizations, SCTransform, corrected counts from SCTransform or background corrected counts) can all be stored as `layers` in your `AnnData` object (as long as they maintain the same dimensionality of the main expression matrix used for visualization).
 
### Format conversion
 
There are a handful of tools that can be used to convert different single cell formats (`seurat`, `loom`, `sce`, `cds`) to `AnnData`. We have had the most success with the R package [sceasy](https://github.com/cellgeni/sceasy). For relatively straightforward conversions, it is suitable to use the `sceasy::convertFormat()` function as specified in the sceasy documentation. For scenarios where you have multiple assays in the same object (as you may encounter in the `seurat` toolchain) it is recommended to check out the function definitions of the methods that are called by `convertFormat()`. For instance, if you look at the source code for `sceasy::seurat2anndata()`, you will see that the extra parameters `assay`, `main_layer`, and `transfer_layers` can be specified for control over which elements of the `seurat` object are carried into the `AnnData` object. Check out the following function definitions and potentially use as a starting point for more customized conversions:

- [`seurat2anndata()`](https://github.com/cellgeni/sceasy/blob/f8f0628a280e0880ea94b00100b463e1f6ba1994/R/functions.R#L14)
- [`sce2anndata()`](https://github.com/cellgeni/sceasy/blob/f8f0628a280e0880ea94b00100b463e1f6ba1994/R/functions.R#L64)
- [`loom2anndata()`](https://github.com/cellgeni/sceasy/blob/f8f0628a280e0880ea94b00100b463e1f6ba1994/R/functions.R#L116)


**Note:** While `AnnData` is able to accommodate multiple representations of an expression matrix in the object, matrices that are stored in `adata.layers` are required to be the same dimensions as `adata.X` (`adata.raw.X` may be of a different dimensionality). 

<br/>

### Alternative assays

The cellxgene data portal and explorer are able to handle a wide range of single cell data types. Due to requirements of the cellxgene schema and limitations of the `AnnData` object structure, there are some assay specific considerations that are outlined below

- ATAC/mC
  - raw data: Since there is not a standard way to generate a counts matrix or gene activity matrix for single cell epigenetic data, it is suitable to put a copy of `adata.X` in place of a raw data matrix. You can specify this assignment in `adata.uns`
  - [(Not Live)Link to example curated dataset](www.example.com)
- RNAseq
  -  In some single cell toolchains, outputs of different computational methods may be of a different dimensionality than the input matrix (i.e SCTransform, or the scale.data slot in Seurat objects). Since `AnnData` objects do not allow for matrices of different dimensions in `adata.layers`, it is suitable to pad the missing rows/features with zeros to ensure that indices (feature names) across different the different layers are matched. You can specify that these matrices were padded in `adata.uns['layer_descriptions']` as a part of the cellxgene schema. This allows users who may reuse your data to remove the padding from the matrix before working with the data further.
  -  [(Not Live)Link to example curated dataset](www.example.com)
- CITEseq
  - TODO: considerations for CITE-seq...
  - [(Not Live)Link to example curated dataset](www.example.com)
- Spatial/Visium
  - Since spatial sequencing spots are not at the single cell level, you may wish to include prediction matrices that provide deconvolution information for each spot. These can generally be stored in `adata.uns` although it can be interesting to add prediction scores for each spot as new columns to the `adata.obs` dataframe (to get a cell type score for each spot)
  - [(Not Live)Link to example curated dataset](www.example.com)

<br/>

</details>
---

<details>
  <summary>Appendix 2: Schema Definition </summary>

## Schema definition

<br/>

The curation process requires that we collect metadata for a few fields. These fields are defined by the [corpora schema](corpora_schema.md) to support easy search and integration across datasets. These metadata are stored within the `AnnData` object (in `adata.uns` and `adata.obs`).In this section, we are covering:

 - what metadata are required to adhere to the cellxgene schema and
 - the locations of these metadata in the `AnnData` object. 

While it is possible to build an object that is acceptable by the cellxgene schema manually, in the [next section](#cellxgene-curation-tools), we will show you how to perform this curation with tools provided by the CZI curation team.

#### `uns`

`adata.uns` is a dictionary of key-value pairs that we use to store dataset level metadata. This is also where you store information on the different data representations that you are sharing (i.e you are sharing a raw counts matrix, normalized expression matrix, scaled and centered expression matrix). Right now, the `uns` slot in our pbmc3k dataset in not populated, but after curation it would contain the following key - value pairs:

<br/>

<div align="center">

| Key                                      | Value      | Description                                                                     | Example               |
| ---------------------------------------- |:----------:| :-------------------------------------------------------------------------------| ----------------------|
| `adata.uns['version']`                   | string     | schema version (current version is `1.1.0`)                                     | `1.1.0`               |
| `adata.uns['title']`                     | string     | title of publication (and title of dataset if more than one is being submitted) | `10X PBMC`            |
| `adata.uns['publication_doi']`           | string     | DOI of preprint or official publication                   | `https://doi.org/00.0000/2021.01.01.000000` |
| `adata.uns['organism']`                  | string     | name of organism (options are `Human` or `Mouse`)                               | `Human`               |
| `adata.uns['organism_ontology_term_id']` | string     | NCBI Taxon ID                                                                   | `NCBITaxon:9606`      |
| `adata.uns['layer_descriptions']`        | dictionary | a set of key value pairs defining the locations of different representations in the `AnnData` object and a description of that representation. One of these key-value pairs must be `raw: raw`  | `{'X': 'log1p' (this value can be free text)}` |

</div>

<div align="center">
  <b> Table: </b> Required and suggested dataset level metadata
</div>

<br/>

In the above table we see what type on information is necessary at the dataset level. You should pay extra attention to the `title` and the `layer_descriptions` keys. More specifically:
 - `title` should contain the name of the publication that the dataset is coming from. If there is more than one dataset associated with the publication, then the name of the dataset should be appended to the title (i.e. `'title': 'publication name: dataset title'`). This field will be used to name and distinguish different datasets in the same collection in the portal.
 - At least one of the key value pairs in `layer_descriptions` needs to indicate the presence of a raw counts layers. This must be specified as `raw: raw` (which implies that `adata.layers['raw']` is where your raw counts are stored). Additional layers may be specified and the values for these keys may be as descriptive as necessary (i.e. `scale.data: scaled and centered data`).

<br/>

#### `obs`

`adata.obs` is a dataframe that is used to store cell level metadata. In addition to experimental metadata, the following additional fields are required:

| `obs` column name                    | Type   | Description                                                       | Example                                   |
| :----------------------------------- |:------:| -----------------------------------------------------------------:|:------------------------------------------|
| 'tissue'                             | string | UBERON term                                                       | `blood`.                                  |
| 'tissue_ontology_term_id'            | string | UBERON term id                                                    | `UBERON:0000178`                          |
| 'assay'                              | string | EFO term                                                          | `scRNA-seq`                               |   
| 'assay_ontology_term_id'             | string | EFO term id                                                       | `EFO:0008913`                             |
| 'disease'                            | string | MONDO term or `normal`                                            | `normal`                                  |
| 'disease_ontology_term_id'           | string | MONDO term id for a disease or `PATO:0000461` for normal          | `PATO:0000461`                            |
| 'cell_type'                          | string | CL term                                                           | `megakaryocyte`                           |
| 'cell_type_ontology_term_id'         | string | CL term id                                                        | `CL:0000556`                              |
| 'sex'                                | string | `male`, `female`, `mixed`, `unknown`, or `other`                  | `unknown`                                 |
| 'ethnicity'                          | string | HANCESTRO term, `na` if non-human, `unknown` if not available     | `unknown`                                 |
| 'ethnicity_ontology_term_id	'        | string | HANCESTRO term id, `na` if non-human                              | `unknown`                                 |
| 'development_stage'                  | string | HsapDv term, `unknown` if not available                           | `Human Adult Stage`                       |
| 'development_stage_ontology_term_id	'| string | HsapDv term id if human, child of `EFO:0000399` otherwise         | `HsapDv:0000087`                          |

<div align="center">
  <b> Table: </b> Required cell level metadata
</div>

<br/>

**Note:** When annotating the `cell_type` field you should generally try to find the most specific ontology term that accurately represents a given cell type in your dataset.

**Note** The `tissue` field must be appended with " (cell culture)" or " (organoid)" if appropriate
<br/>

---

</details>
