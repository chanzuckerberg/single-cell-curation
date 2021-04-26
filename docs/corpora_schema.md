# cellxgene Data Integration Schema

Contact: acarr@chanzuckerberg.com

Document Status: _Approved_

Version: 1.2.0

Date Last Modified: 2020-12-11

The key words "MUST", "MUST NOT", "REQUIRED", "SHALL", "SHALL NOT", "SHOULD", "SHOULD NOT", "RECOMMENDED", "NOT RECOMMENDED" "MAY", and "OPTIONAL" in this document are to be interpreted as described in [BCP 14](https://tools.ietf.org/html/bcp14), [RFC2119](https://www.rfc-editor.org/rfc/rfc2119.txt), and [RFC8174](https://www.rfc-editor.org/rfc/rfc8174.txt) when, and only when, they appear in all capitals, as shown here.

## Background

cellxgene aims to support the publication, sharing, and exploration of single-cell datasets.
Building on those published datasets, cellxgene seeks to create references of the phenotypes and composition of cells that make up human tissues.
Creating references from multiple datasets requires some harmonization of metadata and features in the cellxgene Data Portal. But if that harmonization is too onerous, it will burden the goal of rapid data sharing.

cellxgene balances publishing and reference creation needs by requiring datasets hosted in the cellxgene Data Portal to follow a small schema with only a few required fields.
These fields are expected to be very useful for data integration and also simply and readily available from data submitters.

Note that the requirements in the schema are just the minimum required information. Datasets often have additional metadata, which is preserved in datasets submitted to the Data Portal.

## Purpose

This document describes the schema, a type of contract, that the cellxgene application expects all files to adhere to so that it can enable searching, filtering, and integration of datasets it hosts.
It is intentionally general because cellxgene supports multiple download formats, each of which adheres to this schema.

Users interested in creating a new format that adheres to the cellxgene schema or who want to understand what metadata all cellxgene files can be expected to carry should read this document.

Users interested in the schema's implementation should review the [encoding documents](https://github.com/chanzuckerberg/single-cell-curation/tree/main/docs/encodings/), which describe how cellxgene expects data to be uploaded in AnnData, and the formats of cellxgene's downloadable matrix files.

Users interested in converting and uploading data to cellxgene should review the [schema_guide](https://github.com/chanzuckerberg/single-cell-curation/blob/main/docs/schema_guide.md), which describes a process for curating datasets that adhere to this schema.

## Overview

This schema describes data that measure the phenotypes of cells. It supports multiple assay types, but each assay takes the form of one or more two dimensional matrices whose values are quantitative measures of the genes of cells.

The schema additionally describes how the dataset, genes, and cells should be annotated to describe the biological and technical characteristics of the data.

This document is split into six main sections:
* [General requirements](#general-requirements)
* [Matrix layers](#matrix-layers), which describe the data required for different assays
* [Cell metadata](#cell-metadata), which describe each cell in the dataset
* [Gene metadata](#gene-metadata), which describe each gene in the dataset
* [Dataset metadata](#dataset-metadata), which describe the dataset as a whole
* [Presentation metadata](#presentation-metadata), which are used by the application to adjust the presentation of submitted datasets.

### General Requirements

*   **No PII**. Curators agree to this requirement as part of the data submission policy.
    However, it is not strictly enforced in our validation tooling because it is difficult for software to predict what is and is not PII.
    It is up to the submitter to ensure that no metadata can be personally identifiable: no names, dates of birth, specific locations, etc.
    There's a [list](https://docs.google.com/document/d/1sboOmbafvMh3VYjK1-3MAUt0I13UUJfkQseq8ANLPl8/edit).

### Matrix Layers

cellxgene's data requirements are tailored to optimize data reuse. Because each assay has different characteristics, our requirements differ by assay type. In general,
cellxgene requires submission of "raw" data suitable for computational reuse when a standard form exists, and strongly suggests that a "final" matrix suitable for
visualization in the explorer be included. So that cellxgene's data can be provided in download formats suitable for both R and Python, the schema imposes the following requirements:

*   All matrix layers MUST have the same shape, and have the same cell labels and gene labels.
*   Because it is impractical to retain all barcodes in raw and final matrices, any low quality cell filtering MUST be applied to both.
    By contrast, those wishing to reuse datasets require access to raw gene expression values, so genes MUST be present in both datasets.
    Summarizing, any cell barcodes that are removed from the data MUST be filtered from both raw and final matrices and genes MUST NOT be filtered from the raw matrix.
*   Any genes that publishers wish to filter from the final matrix MAY have their values replaced by a language appropriate "null" value (e.g. `np.nan` for python), which will mask them from exploration.
*   Additional layers provided at author discretion MAY be stored using author-selected keys, but MUST have the same cells and genes as other layers.

In addition to these general requirements, the following table describes the matrix data and layers requirements that are assay-specific.
If cellxgene does not support an assay you would like to publish, please post an issue on this repository to start a conversation about extending the schema.
If an entry in the table is empty, the cellxgene schema does not have any other requirements on data in those layers beyond the ones listed above.
This is usually the case when there are many ways to produce the matrix layer in question.

| Assay                                 | "raw" required? | "raw" requirements                                                                                                                            | "final" required?     | "final" requirements                                                                                 | Other layers |
|---------------------------------------|-----------------|-----------------------------------------------------------------------------------------------------------------------------------------------|-----------------------|------------------------------------------------------------------------------------------------------|--------------|
| scRNA-seq (UMI, e.g. 10x v3)          | REQUIRED        | Values MUST be de-duplicated molecule counts.                                                                                                 | STRONGLY RECOMMENDED  |                                                                                                      | OPTIONAL     |
| scRNA-seq (non-UMI, e.g. SS2)         | REQUIRED        | Values MUST be one of read counts (e.g. FeatureCounts) or  estimated fragments (e.g. output of RSEM).                                         | STRONGLY RECOMMENDED  |                                                                                                      | OPTIONAL     |
| Accessibility (e.g. ATAC-seq, mC-seq) | NOT REQUIRED    |                                                                                                                                               | REQUIRED              | Values MUST correspond to gene features (e.g. HGNC symbols if human)                                 | OPTIONAL     |
| CITE-seq                              | REQUIRED        | Values MUST be de-duplicated molecule counts. Protein features MUST suffixed with ":AB". Protein and RNA features MUST both be stored in raw. | STRONGLY RECOMMENDED  | Protein features MUST be suffixed with ":AB". Protein and RNA features MUST both be stored in final. | OPTIONAL     |

### Integration Metadata

cellxgene requires ontology terms to enable search, comparison, and integration of data.
Ontology terms MUST use [OBO-format ID](http://www.obofoundry.org/id-policy.html), meaning they are a CURIE where the prefix identifies the ontology.
For example `EFO:0000001` is a term in the `EFO` ontology.

When no appropriate ontology value is available, then the most precise accurate term MUST be used.
For example if the `cell_type` field describes a relay interneuron, but the most specific available term in the CL ontology is CL:0000099 ("Interneuron"), then the interneuron term can be used to fulfill this requirement, and ensures that users searching for "neuron" are able to find these data.
Users will still be able to access more specific cell type annotations that have been submitted with the data (but aren't required by the schema).
A dataset comprising cells of the human embryo provides a more extreme example.
In this case, the most  precise accurate term may be the root of the cell ontology `cell`, or its child term `cell in vitro`.
The Cell Ontology is expanding over time, and we hope to migrate datasets to more defined terms as they are defined.
In the mean time, using Cell Ontology terms maximizes the findability (and therefore reusability) of datasets.

#### Cell Metadata

Each cell MUST be annotated with the following ontology terms.

**Field name**|**Constraints**
:--|:--
tissue\_ontology\_term\_id|UBERON term. This field SHOULD be appended with " (cell culture)" or " (organoid)" if appropriate.
assay\_ontology\_term\_id|EFO term
disease\_ontology\_term\_id|MONDO term or [PATO:0000461](https://bioportal.bioontology.org/ontologies/PATO/?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FPATO_0000461)
cell\_type\_ontology\_term\_id|CL term
ethnicity\_ontology\_term\_id|HANCESTRO term, "na" if non-human
development\_stage\_ontology\_term\_id|HsapDv term if human, child of EFO:0000399 otherwise

With the exception of sex, which does not adhere to an ontology, the ontology label MUST also be provided and MUST match the paired ontology term:

**Field name**|**Constraints**
:--|:--
tissue|string. This field SHOULD be appended with " (cell culture)" or " (organoid)" if appropriate.
assay|string
disease|string
cell\_type|string
sex|"male", "female", "mixed", "unknown", or "other"
ethnicity|string, "na" if non-human, "unknown" if not available
development\_stage|string, "unknown" if not available

#### Gene Metadata

Cellxgene uses standard gene symbols to ensure that all datasets it stores measure the same features and can therefore be integrated.

Every gene feature MUST be assigned a unique identifier.
This is occasionally not present because of one-to-many mappings between gene symbols and other gene ids.
In cases where there are duplicated feature identifiers, they will need to be appropriately combined before submission.
This is needed for the explorer to function.
When a user requests gene expression information, the explorer needs to be able to unambiguously return a single value.

If the samples being measured are human, then the feature ids MUST be [HGNC](https://www.genenames.org/about/guidelines/#!/#tocAnchor-1-7) approved gene symbols.
If the features are mouse genes, then the feature ids SHOULD be [MGI](http://www.informatics.jax.org/mgihome/nomen/gene.shtml) gene symbols.
For other organisms, gene symbols SHOULD be the accepted standard human readable symbols for that organism.

#### Dataset Metadata

There following fields annotate the whole dataset and MUST be provided.

**Field name**|**Constraints**
:--|:--
organism|String
organism\_ontology\_term\_id|NCBITaxon term
layer\_descriptions|Dictionary\[String\]|Keys MUST be the layer names whose values are free text description of how the layer was created (e.g. "counts per million")
version|Dictionary| MUST contain key `corpora_schema_version` and its value MUST be the schema encoding version. See [here](https://github.com/chanzuckerberg/single-cell-curation/tree/main/docs/encodings/) for documentation that describes the encoding.
batch_condition|String \| List\[String\]| values MUST match cell metadata keys. Together, these keys define the "batches" that a normalization or integration algorithm should be aware of. For example if "patient" "seqBatch" are keys of vectors of cell metadata `"patient"`, `"seqBatch"`, or `["patient", "seqBatch"]` would be valid batch_condition values.


### Presentation Metadata

There are two fields that are required so that the cellxgene Data Portal and Explorer can present datasets appropriately.

* Each dataset MUST have at least one **embedding**, a mapping from each cell to a tuple of floats of length at least 2.
  These are usually generated by algorithms like umap or tsne, but can also represent `(x, y)` coordinates of cells in spatial assays.
  They are used to display the dataset in the Explorer.
* Each dataset MUST have a title.
  This is a string that describes and differentiates the dataset from others in the collection.
  It will be displayed on a page that also has the collection name.
  For example, in the collection [Cells of the adult human heart](https://cellxgene.cziscience.com/collections/b52eb423-5d0d-4645-b217-e1c6d38b2e72), the first dataset name is "All — Cells of the adult human heart".

#### Presentation Hints

The metadata fields below are optional.
They aren't needed for integration, and cellxgene can display the data fine without them, but if they are included cellxgene will do something with them.
This allows submitters to fine-tune how their datasets are presented, which is a common request.


<table>
  <tr>
   <td><strong>Field name</strong>
   </td>
   <td><strong>Description</strong>
   </td>
  </tr>
  <tr>
   <td>color_map
   </td>
   <td>Submitters can include a field called "{field}_colors" for any other categorical integer metadata field. The value must be an array of one of
       fourcolor specifications:
<ul>

<li>CSS4 color name, as supported by matplotlib
    <a href="https://matplotlib.org/3.1.0/gallery/color/named_colors.html">https://matplotlib.org/3.1.0/gallery/color/named_colors.html</a>

<li>RGB tuple/list with values ranging from 0.0 to 1.0, as in [0.5, 0.75, 1.0]

<li>RFB tuple/list with values ranging from 0 to 255, as in [128, 192, 255]

<li>Hex triplet string, as in "#08c0ff"and each string must be a hex color code.

</li>
</ul>
The color code at the nth position in the array corresponds to category n in the metadata field.

   </td>
  </tr>
  <tr>
   <td>default_embedding
   </td>
   <td>Name of the embedding that should be the default for display.
   </td>
  </tr>
</table>
