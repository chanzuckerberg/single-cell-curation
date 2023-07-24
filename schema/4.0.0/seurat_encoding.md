
# Schema - Seurat encoding

Contact: brianraymor@chanzuckerberg.com

Document Status: _Draft_

Schema version: 4.0.0


## Overview

All data submitted to [CELLxGENE Discover](https://cellxgene.cziscience.com/) is automatically converted to a Seurat V4 object that can be loaded by the R package [Seurat](https://satijalab.org/seurat/).

This document describes the Seurat encoding for the converted data. <u><strong>Readers should be familiar with the [schema](./schema.md).</strong></u>

## Encoding

The Seurat V4 object is stored using the [native serialization RDS format](https://cran.r-project.org/doc/manuals/r-release/R-ints.html#Serialization-Formats). The [serialization interfaces](https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html) allow the object to be serialized to and restored from disk. 

A `local.rds` file downloaded from the cellxgene Data Portal can be read into an R session with the following code.

```r
seurat_object <- readRDS("local.rds")
```

The following sections describe the individual components of a dataset as encoded in the Seurat object.

### Data matrix

Matrix data is stored in the slot `assays` under the element `RNA`. Seurat allows elements of `assays` to have any arbitrary name, but all matrix data from cellxgene will be stored in `RNA`.

#### Raw data

<table><tbody>
    <tr>
      <th>Slot</th>
      <td><code>counts</code></td>
    </tr>
    <tr>
      <th>Access</th>
      <td><code>seurat_object@assays$RNA@counts</code></td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
        <code>Matrix::dgCMatrix</code> object for sparse matrices
        <br>
        <code>base::matrix</code> object for dense matrices
        </td>
    </tr>
</tbody></table>
<br>

#### Final (normalized) data

<table><tbody>
    <tr>
      <th>Slot</th>
      <td><code>data</code></td>
    </tr>
    <tr>
      <th>Access</th>
      <td><code>seurat_object@assays$RNA@data</code></td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
        <code>Matrix::dgCMatrix</code> object for sparse matrices
        <br>
        <code>base::matrix</code> object for dense matrices
        </td>
    </tr>
</tbody></table>
<br>


### Cell Metadata

It is stored as a `data.frame` in the `meta.data` slot.

<table><tbody>
    <tr>
      <th>Access</th>
      <td><code>seurat_object@meta.data</code></td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>data.frame</code></td>
    </tr>
</tbody></table>
<br>

This `data.frame` will have the following columns required by the schema. There may be additional, optional columns from the original study.

<table><tbody>
    <tr>
      <th>Column</th>
      <th>Value</th>
    </tr>
    <tr>
      <td>assay_ontology_term_id</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>assay</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>cell_type_ontology_term_id</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>cell_type</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>development_stage_ontology_term_id</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>development_stage</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>disease_ontology_term_id</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>disease</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>donor_id</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>self_reported_ethnicity_ontology_term_id</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>self_reported_ethnicity</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>is_primary_data</td>
      <td><code>logical</code></td>
    </tr>
    <tr>
      <td>organism_ontology_term_id</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>organism</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>sex_ontology_term_id</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>sex</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>suspension_type</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>tissue_ontology_term_id</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>tissue</td>
      <td><code>character</code></td>
    </tr>
</tbody></table>
<br>


### Gene Metadata

It is stored as a `data.frame` in the `meta.features` slot of the element `RNA` of the slot `assays`.

<table><tbody>
    <tr>
      <th>Access</th>
      <td><code>seurat_object@assays$RNA@meta.features</code></td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>data.frame</code></td>
    </tr>
</tbody></table>
<br>

This `data.frame` will have the following columns required by the schema. There may be additional, optional columns from the original study.

<table><tbody>
    <tr>
      <th>Column</th>
      <th>Value</th>
    </tr>
    <tr>
      <td>gene_id (<code>row.names</code>)</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>feature_biotype</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>feature_is_filtered</td>
      <td><code>logical</code></td>
    </tr>
    <tr>
      <td>feature_name</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>feature_reference</td>
      <td><code>character</code></td>
    </tr>
 </tbody></table>
<br>


### Embeddings

Each available embedding is stored as a named element of the slot `reductions`.

<table><tbody>
    <tr>
      <th>Access</th>
      <td><code>seurat_object@reductions$[embedding_name]</code> where <code>embedding name</code> is the name of embedding, usually "tsne" or "umap"   </td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code> SeuratObject::DimenReduc</code> object</td>
    </tr>
</tbody></table>
<br>

The matrix with reductions is stored in the `cell.embeddings` slot of `SeuratObject::DimenReduc`. Rows correspond to cells and are named with cell IDs. Columns correspond to reductions and are named `[key]_[number]`, where `key` is an alphanumeric string and `number` is an integer, e.g. `PC_1`.

### Dataset Metadata

Dataset metadata is stored in the `misc` slot.

<table><tbody>
    <tr>
      <th>Access</th>
      <td><code>seurat_object@misc</code></td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>list[character|list]</code></td>
    </tr>
</tbody></table>
<br>

Only fields defined in the [schema](./schema.md/#uns-dataset-metadata) are transferred to the Seurat object:


<table><tbody>
    <tr>
      <th>Name</th>
      <th>Value</th>
      <th>Optional</th>
    </tr>
    <tr>
      <td>schema_version</td>
      <td><code>character</code></td>
      <td>No</td>
    </tr>
    <tr>
      <td>title</td>
      <td><code>character</code></td>
      <td>No</td>
    </tr>
    <tr>
      <td>batch_condition</td>
      <td><code>list[character]</code></td>
      <td>Yes</td>
    </tr>
    <tr>
      <td>default_embedding</td>
      <td><code>character</code></td>
      <td>Yes</td>
    </tr>
    <tr>
      <td>X_approximate_distribution</td>
      <td><code>character</code></td>
      <td>Yes</td>
    </tr>
 </tbody></table>
<br>

