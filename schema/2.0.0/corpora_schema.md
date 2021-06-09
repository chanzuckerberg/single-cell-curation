# cellxgene Data Integration Schema

Contact: acarr@chanzuckerberg.com

Document Status: _Approved_

Version: 2.0.0

Date Last Modified: 2021-06-07

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

This schema describes data that measure the phenotypes of cells.  
It supports multiple assay types, but each assay takes the form of one or more two dimensional matrices whose values are quantitative measures of the genes of cells.

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

cellxgene's data requirements are tailored to optimize data reuse.
Because each assay has different characteristics, our requirements differ by assay type.
In general, cellxgene requires submission of "raw" data suitable for computational reuse when a standard form exists, and strongly suggests that a "final" matrix suitable for visualization in the explorer be included.
So that cellxgene's data can be provided in download formats suitable for both R and Python, the schema imposes the following requirements:

*   All matrix layers MUST describe the same cells, and have the same cell labels.
*   Because it is impractical to retain all barcodes in raw and final matrices, any low quality cell filtering MUST be applied to both.
*   Because those wishing to reuse datasets require access to raw gene expression values, genes MUST NOT be filtered from raw data, and genes present in final data MUST match or be a subset of genes in raw.
*   Additional layers provided at author discretion MAY be stored using author-selected keys, but MUST have the same cells and genes as the final layer.

In addition to these general requirements, the following table describes the matrix data and layers requirements that are assay-specific.
If cellxgene does not support an assay you would like to publish, please post an issue on this repository to start a conversation about extending the schema.
If an entry in the table is empty, the cellxgene schema does not have any other requirements on data in those layers beyond the ones listed above.
This is usually the case when there are many ways to produce the matrix layer in question.

| Assay                                 | "raw" required? | "raw" requirements                                                                                                                            | "final" required?     | "final" requirements                                                                                 | Other layers |
|---------------------------------------|-----------------|-----------------------------------------------------------------------------------------------------------------------------------------------|-----------------------|------------------------------------------------------------------------------------------------------|--------------|
| scRNA-seq (UMI, e.g. 10x v3)          | REQUIRED        | Values MUST be de-duplicated molecule counts.                                                                                                 | STRONGLY RECOMMENDED  |                                                                                                      | OPTIONAL     |
| scRNA-seq (non-UMI, e.g. SS2)         | REQUIRED        | Values MUST be one of read counts (e.g. FeatureCounts) or  estimated fragments (e.g. output of RSEM).                                         | STRONGLY RECOMMENDED  |                                                                                                      | OPTIONAL     |

### Integration Metadata

cellxgene requires ontology terms to enable search, comparison, and integration of data.
Ontology terms for cell metadata MUST use [OBO-format ID](http://www.obofoundry.org/id-policy.html), meaning they are a CURIE where the prefix identifies the ontology.
For example `EFO:0000001` is a term in the `EFO` ontology.

When no appropriate ontology value is available, then the most precise accurate term MUST be used.
For example if the `cell_type` field describes a relay interneuron, but the most specific available term in the CL ontology is CL:0000099 ("Interneuron"), then the interneuron term can be used to fulfill this requirement, and ensures that users searching for "neuron" are able to find these data.
Users will still be able to access more specific cell type annotations that have been submitted with the data (but aren't required by the schema).
A dataset comprising cells of the human embryo provides a more extreme example.
In this case, the most  precise accurate term may be the root of the cell ontology `cell`, or its child term `cell in vitro`.
The Cell Ontology is expanding over time, and we hope to migrate datasets to more defined terms as they are defined.
In the mean time, using Cell Ontology terms maximizes the findability (and therefore reusability) of datasets.


### Required Ontologies

| Ontology | Required version | Download |
|:--|:--|:--|
| [NCBI organismal classification] | [2021-02-15] | [ncbitaxon.owl] |
| [Uberon multi-species anatomy ontology] | [2021-02-12] |  [composite-vertebrate.owl]  |
| [Cell Ontology (CL)] | included with Uberon | included with Uberon |
| [Experimental Factor Ontology (EFO)] | [2021-05-17 EFO 3.30.0] | [efo-base.owl] |
| [Mondo Disease Ontology (MONDO)] | [2021-06-01] | [mondo.owl] | 
| [Human Ancestry Ontology (HANCESTRO)] | 2021-01-04 (2.5) | hancestro.owl from [OLS]  (Note: [HANCESTRO releases] do not include .owl files) |
| | | |

[NCBI organismal classification]: http://obofoundry.org/ontology/ncbitaxon.html
[2021-02-15]: https://github.com/obophenotype/ncbitaxon/releases/tag/v2021-02-15
[ncbitaxon.owl]: https://github.com/obophenotype/ncbitaxon/releases/download/v2021-02-15/ncbitaxon.owl.gz

[Uberon multi-species anatomy ontology]: http://www.obofoundry.org/ontology/uberon.html
[2021-02-12]: https://github.com/obophenotype/uberon/releases/tag/v2021-02-12
[composite-vertebrate.owl]: https://github.com/obophenotype/uberon/releases/download/v2021-02-12/composite-vertebrate.owl

[Cell Ontology (CL)]: http://obofoundry.org/ontology/cl.html

[Experimental Factor Ontology (EFO)]: http://www.ebi.ac.uk/efo
[2021-05-17 EFO 3.30.0]: https://github.com/EBISPOT/efo/releases/tag/v3.30.0
[efo-base.owl]: https://github.com/EBISPOT/efo/releases/download/v3.30.0/efo-base.owl

[Mondo Disease Ontology (MONDO)]: http://obofoundry.org/ontology/mondo.html

[2021-06-01]: https://github.com/monarch-initiative/mondo/releases/tag/v2021-06-01

[mondo.owl]: https://github.com/monarch-initiative/mondo/releases/download/v2021-06-01/mondo.owl

[Human Ancestry Ontology (HANCESTRO)]: http://www.obofoundry.org/ontology/hancestro.html
[OLS]: https://www.ebi.ac.uk/ols/ontologies/hancestro
[HANCESTRO releases]: https://github.com/EBISPOT/ancestro/releases

#### Cell Metadata

Each cell MUST be annotated with the following ontology terms by the curator.

| Field name | Constraints |  Annotator |
:--|:--|:--
| organism_ontology_term_id | NCBITaxon term | Curator |
| tissue_ontology_term_id | UBERON term. This SHOULD be appended with " (cell culture)" or " (organoid)" if appropriate. | Curator |
| assay_ontology_term_id | EFO term | Curator |
| disease_ontology_term_id | MONDO term or [PATO:0000461](http://purl.obolibrary.org/obo/PATO_0000461) | Curator |
| cell_type_ontology_term_id | CL term | Curator |
| ethnicity_ontology_term_id | `organism_ontolology_term_id` is “NCBITaxon:9606”. This MUST be either a HANCESTRO term or “unknown” if unavailable. <br><br> `organism_ontolology_term_id` is NOT “NCBITaxon:9606”. This MUST be "na". | Curator |
| development_stage_ontology_term_id | If unavailable, this MUST be "unknown". <br><br> `organism_ontolology_term_id` is “NCBITaxon:9606. This MUST be a HsapDv term.<br><br> `organism_ontolology_term_id` is NOT “NCBITaxon:9606”. This MUST be a child of "EFO:0000399". | Curator |
| sex_ontology_term_id | This MUST be a child of [PATO:0001894](http://purl.obolibrary.org/obo/PATO_0001894) | | Curator |
| | |

The ontology labels MUST be applied by curation software and MUST match the paired ontology term:

| Field name | Constraints | Annotator |
:--|:--|:--
| organism | string | Software |
| tissue | string. This MUST be appended with " (cell culture)" or " (organoid)" if set in the matching identifier. | Software |
| assay | string | Software |
| disease | string | Software |
| cell_type | string | Software |
| ethnicity | string. This MUST be "na" or "unknown" if set in the matching identifier. | Software |
| development_stage | string. This MUST be "unknown" if set in the matching identifier | Software |
| sex | string | Software |
| | |

#### Gene Metadata


Cellxgene uses ENSEMBL gene identifiers to ensure that all datasets it stores measure the same features and can therefore be integrated.

| Organism | Required version | Download |
|:--|:--|:--|
| [ENSEMBL (Human)] | Human reference GRCh38 (GENCODE v32/Ensembl 98)] matching [cellranger 2020-A (July 7, 2020) release] | [gencode.v32.primary_assembly.annotation.gtf] |
| [ENSEMBL (Mouse)] | Mouse reference mm10 (GENCODE vM23/Ensembl 98) matching [cellranger 2020-A (July 7, 2020) release]| [gencode.vM23.primary_assembly.annotation.gtf] |

[ENSEMBL (Human)]: http://www.ensembl.org/Homo_sapiens/Info/Index
[gencode.v32.primary_assembly.annotation.gtf]:http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"

[ENSEMBL (Mouse)]: http://www.ensembl.org/Mus_musculus/Info/Index
[gencode.vM23.primary_assembly.annotation.gtf]:http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz"

[cellranger 2020-A (July 7, 2020) release]: https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build

Every gene feature MUST be assigned a unique identifier.
This is occasionally not present because of one-to-many mappings between gene symbols and other gene identifiers.
In cases where there are duplicated feature identifiers, they MUST be appropriately combined before submission.
This is required for the explorer to function.
When a user requests gene expression information, the explorer must unambiguously return a single value.

| Field name | Constraints | Annotator |
:--|:--|:--
| feature_biotype | "gene" | Software |
| feature_id | ENSEMBL term from the required version corresponding to the organism for the gene feature  | Curator |
| feature_name | string. ENSEMBL gene name | Software |
|||

#### Dataset Metadata

The following fields annotate the whole dataset and MUST be provided.

**Field name**|**Constraints**
:--|:--
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

## Appendix A. Changelog

* Several editorial changes were made:
  * [#50](https://github.com/chanzuckerberg/single-cell-curation/issues/50), [#28](https://github.com/chanzuckerberg/single-cell-curation/issues/28), [#21](https://github.com/chanzuckerberg/single-cell-curation/issues/21): Removed AnnData encoding details from the schema to ensure it generalizes, and added editorial content describing why to read the schema, schema guide, and encoding documents.
  * Clarified that the meaning of "must", "should", and select other words have defined, standard meaning.
  * [#45](https://github.com/chanzuckerberg/single-cell-curation/issues/45) Updated reference to new PII content

* Empty ontology fields are no longer permitted:
  * If the source of cells is cell culture or organoid, the cell_type field **cannot** be left empty.
  * `ontology_term_id` fields may **not** be empty strings when no appropriate ontology value is available.

* `tags` and `default_field` presentation metadata are not used by the application and have been deprecated.

* Changed one metadata field:
  * [#7](https://github.com/chanzuckerberg/single-cell-curation/issues/7), [#8](https://github.com/chanzuckerberg/single-cell-curation/issues/8) Move organism from "dataset" metadata to "cell" metadata

* Added `sex_ontology_term_id` and ontology

* Updated the **Gene Metadata** section to require ENSEMBL identifiers instead of HGNC symbols. Added `gene_id` and `gene_symbol`.

* Specifies the ontology versions for validation


