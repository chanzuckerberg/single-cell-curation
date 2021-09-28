
# Schema - Seurat encoding

Contact: acarr@chanzuckerberg.com

Document Status: _In progress_

Schema version: 2.0.0


## Overview

All data submitted to the [cellxgene Data Portal](https://cellxgene.cziscience.com/) is automatically converted to a Seurat V3 object that can be loaded by the R package [Seurat](https://satijalab.org/seurat/).

This document describes the Seurat encoding for the converted data. Readers should be familiar with the [schema](.schema.md). 

## Encoding

The Seurat V3 object is stored using the [native serialization RDS format](https://cran.r-project.org/doc/manuals/r-release/R-ints.html#Serialization-Formats). The [serialization interfaces](https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html) allow the object to be serialized to and restored from disk. 

A `local.rds` file downloaded from the cellxgene Data Portal can be read into an R session with the following code.

```r
seurat_object <- readRDS("local.rds")
```

The following sections describe the individual components of a dataset as encoded in the Seurat object.

### Data matrix

Matrix data is stored in the slot `assays` under the element `RNA`.

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
        <td>Sparse matrix of type <code>dgCMatrix</code></td>
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
        <td>Sparse matrix of type <code>dgCMatrix</code></td>
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

This `data.frame` will have at least the following columns.

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
      <td>ethnicity_ontology_term_id</td>
      <td><code>character</code></td>
    </tr>
    <tr>
      <td>ethnicity</td>
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

This `data.frame` will have at least the following columns:

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
        <td><code>DimenReduc</code> object</td>
    </tr>
</tbody></table>
<br>
