
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

cellxgene balances publishing and reference creation needs by requiring datasets hosted in the cellxgene Data Portal to include a small set of metadata readily available from data submitters.

This document describes the schema, a type of contract, that cellxgene requires all datasets to adhere to so that it can enable searching, filtering, and integration of datasets it hosts.

Note that the requirements in the schema are just the minimum required information. Datasets often have additional metadata, which is preserved in datasets submitted to the cellxgene Data Portal.

## Overview

This schema supports multiple assay types. Each assay takes the form of one or more two-dimensional matrices whose values are quantitative measures of the phenotypes of cells.

The schema additionally describes how the dataset, genes, and cells are annotated to describe the biological and technical characteristics of the data.

This document is organized by:
* [General requirements](#general-requirements)
* [`X` (Matrix layers)](#x-(matrix-layers)), which describe the data required for different assays
* [`obs` (Cell metadata)](#obs-(cell-metadata)), which describe each cell in the dataset
* [`var` (Gene metadata)](#var-(gene-metadata)), which describe each gene in the dataset
* [`obsm` (Embeddings)](#obsm-(embeddings)), which describe embeddings for each dataset
* [`uns` (Dataset metadata](#uns-(dataset-metadata)), which describe the dataset as a whole

## General Requirements

* **AnnData** - The canonical data format for the cellxgene Data Portal is HDF5-backed [AnnData](https://anndata.readthedocs.io/en/latest) as written by version 0.7 of the anndata library.  Part of the rationale for selecting this format is to allow cellxgene to access both the data and metadata within a single file. The schema requirements and definitions for the AnnData `X`, `uns`, `obs`, and `obsm` attributes are described below.

*   **No PII**. Curators agree to this requirement as part of the data submission policy.
    However, it is not strictly enforced in our validation tooling because it is difficult for software to predict what is and is not PII.
    It is up to the submitter to ensure that no metadata can be personally identifiable: no names, dates of birth, specific locations, etc.
    See this [list](https://docs.google.com/document/d/1sboOmbafvMh3VYjK1-3MAUt0I13UUJfkQseq8ANLPl8/edit) for guidance.

#### *Note on types*
The types below are python3 types. Note that a python3 `str` is a sequence of Unicode code points, which is stored null-terminated and UTF-8-encoded by anndata.

## `X` (Matrix Layers)

cellxgene does not impose any additional constraints on the `X` data matrix. It may be sparse or dense and any numeric [`numpy.dtype`](https://numpy.org/doc/stable/reference/arrays.dtypes.html).

cellxgene's data requirements are tailored to optimize data reuse. Because each assay has different characteristics, the requirements differ by assay type. In general,
cellxgene requires submission of "raw" data suitable for computational reuse when a standard form exists and strongly recommends that a "final" matrix suitable for
visualization in the explorer be included. So that cellxgene's data can be provided in download formats suitable for both R and Python, the schema imposes the following requirements:

*   All matrix layers MUST have the same shape, and have the same cell labels and gene labels.
*   Because it is impractical to retain all barcodes in raw and final matrices, any low quality cell filtering MUST be applied to both.
    By contrast, those wishing to reuse datasets require access to raw gene expression values, so genes MUST be present in both datasets.
    Summarizing, any cell barcodes that are removed from the data MUST be filtered from both raw and final matrices and genes MUST NOT be filtered from the raw matrix.
*   Any genes that publishers wish to filter from the final matrix MAY have their values replaced by a language appropriate "null" value (e.g. [`np.nan`](https://numpy.org/doc/stable/reference/constants.html#numpy.nan) for python), which will mask them from exploration.
*   Additional layers provided at author discretion MAY be stored using author-selected keys, but MUST have the same cells and genes as other layers.

In addition to these general requirements, the following table describes the matrix data and layers requirements that are assay-specific. If an entry in the table is empty, the cellxgene schema does not have any other requirements on data in those layers beyond the ones listed above.
This is usually the case when there are many ways to produce the matrix layer in question.

| Assay | "raw" required? | "raw" requirements | "final" required? | "final" requirements | Other layers |
|-|-|-|-|-|-|
| scRNA-seq (UMI, e.g. 10x v3) | REQUIRED in AnnData.layers["raw"] unless no "final" is provided, then AnnData.X | Values MUST be de-duplicated molecule counts. | STRONGLY RECOMMENDED |  | OPTIONAL |
| scRNA-seq (non-UMI, e.g. SS2) | REQUIRED in AnnData.layers["raw"] unless no "final" is provided, then AnnData.X | Values MUST be one of read counts (e.g. FeatureCounts) or  estimated fragments (e.g. output of RSEM). | STRONGLY RECOMMENDED |  | OPTIONAL |
| Accessibility (e.g. ATAC-seq, mC-seq) | NOT REQUIRED |  | REQUIRED in AnnData.X | Values MUST correspond to ENSEMBL gene identifiers  | OPTIONAL |
||||

## Integration Metadata

cellxgene requires ontology terms to enable search, comparison, and integration of data.
Ontology terms for cell metadata MUST use [OBO-format identifiers](http://www.obofoundry.org/id-policy.html), meaning a CURIE (prefixed identifier) of the form **Ontology:Identifier**.
For example [EFO:0000001](https://www.ebi.ac.uk/ols/ontologies/efo/terms?short_form=EFO_0000001) is a term in the `EFO` ontology.

When no appropriate ontology value is available, then the most precise accurate term MUST be used.
For example if the `cell_type` field describes a relay interneuron, but the most specific available term in the CL ontology is [CL:0000099](https://www.ebi.ac.uk/ols/ontologies/cl/terms?obo_id=CL:0000099) for "interneuron", then the interneuron term can be used to fulfill this requirement and ensures that users searching for "neuron" are able to find these data.
Users will still be able to access more specific cell type annotations that have been submitted with the data (but aren't required by the schema).
A dataset comprising cells of the human embryo provides a more extreme example.
In this case, the most precise accurate term is [CL:0000003](https://www.ebi.ac.uk/ols/ontologies/cl/terms?obo_id=CL:0000003)
for "native cell".

~~The Cell Ontology is expanding over time, and we hope to migrate datasets to more precise terms as they are defined.
In the meantime, using Cell Ontology terms maximizes the findability (and therefore reusability) of datasets.~~

#### Gene Annotation

cellxgene requires ENSEMBL gene identifiers to ensure that all datasets it stores measure the same features and can therefore be integrated.

### Required Ontologies

The following ontology dependencies are *pinned* for this version of the schema.

**EDITOR NOTE**: *[Remove before publishing] There are more recent versions of hsapdv.owl and mmusdv.owl on the OLS site, but these downloads are missing a specific version number for pinning.*

| Ontology | OBO Prefix | Required version |
|:--|:--|:--|
| [Cell Ontology] | CL | [cl.owl] : [2021-04-22]|
| [Experimental Factor Ontology] | EFO | [efo.owl] : [2021-05-17 EFO 3.30.0]
| [Human Ancestry Ontology] | HANCESTRO |[hancestro.owl] : [2021-01-04 (2.5)] |
| [Human Developmental Stages] |  HsapDv | [hsapdv.owl] : [2016-07-06 (0.1)] |
| [Mondo Disease Ontology] | MONDO |[mondo.owl] : [2021-06-01] |
| [Mouse Developmental Stages]| MmusDv |  [mmusdv.owl] : [2016-07-06 (0.1)] |
| [NCBI organismal classification] |  NCBITaxon | [ncbitaxon.owl] : [2021-02-15] |
| [Phenotype And Trait Ontology] | PATO | [pato.owl] : [2021-05-26] |  |
| [Uberon multi-species anatomy ontology] |  UBERON | [uberon.owl] : [2021-02-12] |
| | | |

[Cell Ontology]: http://obofoundry.org/ontology/cl.html
[2021-04-22]: https://github.com/obophenotype/cell-ontology/releases/tag/v2021-04-22
[cl.owl]: https://github.com/obophenotype/cell-ontology/blob/v2021-04-22/cl.owl

[Experimental Factor Ontology]: http://www.ebi.ac.uk/efo
[2021-05-17 EFO 3.30.0]: https://github.com/EBISPOT/efo/releases/tag/v3.30.0
[efo.owl]: https://github.com/EBISPOT/efo/releases/download/v3.30.0/efo.owl

[Human Ancestry Ontology]: http://www.obofoundry.org/ontology/hancestro.html
[2021-01-04 (2.5)]: https://github.com/EBISPOT/ancestro/releases/tag/2.5
[hancestro.owl]: https://github.com/EBISPOT/ancestro/blob/2.5/hancestro.owl

[Human Developmental Stages]: http://www.obofoundry.org/ontology/hsapdv.html
[hsapdv.owl]: https://github.com/obophenotype/developmental-stage-ontologies/blob/0.1/src/hsapdv/hsapdv.owl
[2016-07-06 (0.1)]: https://github.com/obophenotype/developmental-stage-ontologies/releases/tag/0.1

[Mondo Disease Ontology]: http://obofoundry.org/ontology/mondo.html
[2021-06-01]: https://github.com/monarch-initiative/mondo/releases/tag/v2021-06-01
[mondo.owl]: https://github.com/monarch-initiative/mondo/releases/download/v2021-06-01/mondo.owl

[Mouse Developmental Stages]: http://obofoundry.org/ontology/mmusdv.html
[mmusdv.owl]: https://github.com/obophenotype/developmental-stage-ontologies/blob/0.1/src/mmusdv/mmusdv.owl

[NCBI organismal classification]: http://obofoundry.org/ontology/ncbitaxon.html
[2021-02-15]: https://github.com/obophenotype/ncbitaxon/releases/tag/v2021-02-15
[ncbitaxon.owl]: https://github.com/obophenotype/ncbitaxon/releases/download/v2021-02-15/ncbitaxon.owl.gz

[Phenotype And Trait Ontology]: http://www.obofoundry.org/ontology/pato.html
[2021-05-26]: https://github.com/pato-ontology/pato/releases/tag/v2021-05-26
[pato.owl]: https://github.com/pato-ontology/pato/blob/v2021-05-26/pato.owl

[Uberon multi-species anatomy ontology]: http://www.obofoundry.org/ontology/uberon.html
[2021-02-12]: https://github.com/obophenotype/uberon/releases/tag/v2021-02-12
[uberon.owl]: https://github.com/obophenotype/uberon/blob/v2021-02-12/uberon.owl

### Required Gene Annotations

The following gene annotation dependencies are *pinned* for this version of the schema.

| Organism | Required version | Download |
|:--|:--|:--|
| [ENSEMBL (Human)] | Human reference GRCh38 (GENCODE v32/Ensembl 98)] matching [cellranger 2020-A (July 7, 2020) release] | [gencode.v32.primary_assembly.annotation.gtf] |
| [ENSEMBL (Mouse)] | Mouse reference mm10 (GENCODE vM23/Ensembl 98) matching [cellranger 2020-A (July 7, 2020) release]| [gencode.vM23.primary_assembly.annotation.gtf] |
|||

[ENSEMBL (Human)]: http://www.ensembl.org/Homo_sapiens/Info/Index
[gencode.v32.primary_assembly.annotation.gtf]:http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"

[ENSEMBL (Mouse)]: http://www.ensembl.org/Mus_musculus/Info/Index
[gencode.vM23.primary_assembly.annotation.gtf]:http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz"

[cellranger 2020-A (July 7, 2020) release]: https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build

## `obs` (Cell Metadata)

`obs` is a [`pandas.DataFrame`](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html).

Curators MUST annotate the following columns in the `obs` dataframe:

| Key | Value | Annotator |
|-|-|-|
| assay_ontology_term_id | `str` or categorical with `str` categories. This MUST be a child of [EFO:0002694](https://www.ebi.ac.uk/ols/ontologies/efo/terms?short_form=EFO_0002694). If there is not an exact match for the assay, *clarifying text* MAY be appended in parentheses such as " (sci-plex)". | Curator |
| cell_type_ontology_term_id | `str` or categorical with `str` categories. This MUST be a CL term. | Curator |
| development_stage_ontology_term_id | `str` or categorical with `str` categories. If unavailable, this MUST be "unknown". <br><br> If `organism_ontolology_term_id` is “NCBITaxon:9606", this MUST be a HsapDv term.<br><br> If `organism_ontolology_term_id` is “NCBITaxon:10090”, this MUST a MmusDv term with the following restrictions: <br><br> **Prenatal stages** MUST be in the range beginning with "MmusDv:0000003" and ending with "MmusDv:0000035". <br><br> **Postnatal stages (1-4 weeks)** MUST be in the range beginning with "MmusDv:0000045" and ending with "MmusDv:0000048". <br><br> **Postnatal stages (after week 4)** MUST be in the range beginning with "MmusDv:0000062" and ending with "MmusDv:0000091". | Curator |
| disease_ontology_term_id | `str` or categorical with `str` categories. This MUST be a MONDO term or [PATO:0000461](http://purl.obolibrary.org/obo/PATO_0000461) | Curator |
| ethnicity_ontology_term_id | `str` or categorical with `str` categories. If `organism_ontolology_term_id` is “NCBITaxon:9606”, this MUST be either a HANCESTRO term or “unknown” if unavailable. <br><br> If `organism_ontolology_term_id` is “NCBITaxon:10090”, this MUST be "na". | Curator |
| is_primary_data | `Bool`. This MUST be `True` if this is the canonical instance of this cellular observation and `False` if not. This is commonly `False` for meta-analyses reusing data or for secondary views of data. | Curator |
| organism_ontology_term_id | `str` or categorical with `str` categories. This MUST be either "NCBITaxon:9606" or "NCBITaxon:10090". | Curator |
| sex_ontology_term_id | `str` or categorical with `str` categories. This MUST be a child of [PATO:0001894](http://purl.obolibrary.org/obo/PATO_0001894) | Curator |
| tissue_ontology_term_id | `str` or categorical with `str` categories. This MUST be an UBERON term. This SHOULD be appended with " (cell culture)" or " (organoid)" if appropriate. | Curator |
| | | |

Curators SHOULD NOT annotate the following columns in the `obs` dataframe which are the human-readable names that MUST match the corresponding ontology term identifier. These columns MUST be automatically annotated by the cellxgene Data Portal when a dataset is uploaded.

| Key | Value | Annotator |
:--|:--|:--
| assay | `str` or categorical with `str` categories. This MUST be appended with any *clarifying text* set in matching identifier. | Data Portal |
| cell_type | `str` or categorical with `str` categories. | Data Portal |
| development_stage | `str` or categorical with `str` categories. This MUST be "unknown" if set in the matching identifier | Data Portal |
| disease | `str` or categorical with `str` categories | Data Portal |
| ethnicity | `str` or categorical with `str` categories. This MUST be "na" or "unknown" if set in the matching identifier. | Data Portal |
| organism | `str` or categorical with `str` categories | Data Portal |
| sex | `str` or categorical with `str` categories | Data Portal |
| tissue | `str` or categorical with `str` categories. This MUST be appended with " (cell culture)" or " (organoid)" if set in the matching identifier. | Data Portal |
| | | |

## `var` (Gene Metadata)

`var` is a [`pandas.DataFrame`](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html).

`var.index` MUST contain unique ENSEMBL gene identifiers for features. 

Curators MUST annotate the following columns in the `var` dataframe:

| Key | Value | Annotator |
:--|:--|:--
| feature_biotype | This MUST be "gene". | Data Portal |
| feature_id (`var.index`) | `str`. This MUST be an ENSEMBL term from the gene annotation corresponding to the organism for the gene.  | Curator |
||||

Curators SHOULD NOT annotate the following column in the `var` dataframe which is the human-readable name that MUST match the corresponding gene identifier. This column MUST be automatically annotated by the cellxgene Data Portal when a dataset is uploaded.

| Key | Value | Annotator |
:--|:--|:--
| feature_name | `str`. This MUST be the ENSEMBL gene name corresponding to the `feature_id`. | Data Portal |
||||

## `obsm` (Embeddings)

For each `str` key, `obsm` stores a `numpy.ndarray` of shape `(n_obs, m)`, where `n_obs` is the number of rows in `X` and `m >= 1`.

To display a dataset in cellxgene Explorer, Curators MUST annotate **one or more** two-dimensional (`m >= 2`) embeddings (tSNE, UMAP, PCA, spatial coordinates) in `obsm`. The key for the embedding MUST be prefixed with `X_`. The text that follows this prefix is presented to users in the cellxgene Explorer *Embedding Choice* selector.

To illustrate, the [Krasnow Lab Human Lung Cell Atlas, 10X dataset](https://cellxgene.cziscience.com/e/krasnow_lab_human_lung_cell_atlas_10x-1-remixed.cxg/) in the [A molecular cell atlas of the human lung from single cell RNA sequencing collection](https://cellxgene.cziscience.com/collections/5d445965-6f1a-4b68-ba3a-b8f765155d3a) defines two embeddings in `obsm`:

* 'X_Compartment_tSNE'
* 'X_tSNE'

Users can then choose whether the 'Compartment_tSNE' or the 'tSNE' embedding  is visualized in cellxgene Explorer:

![Embeddings Illustration](images/embeddings.png
)

See also `default_embedding` in `uns`.

## `uns` (Dataset Metadata)

`uns` is a ordered dictionary with a `str` key. Curators MUST annotate the following keys and values in `uns`:

| Key | Value | Annotator |
:--|:--|:--|
| layer_descriptions | `dict` with `str` keys and values. Keys MUST be the layer names whose values are free text descriptions for how the layer was created (e.g. "counts per million"). One key MUST be "X" which describes the transformations (if any) performed to produce the `X` matrix that cellxgene Explorer displays. | Curator |
| schema_version| `str`. This MUST be "2.0.0". | Curator |
| title | `str`. This text describes and differentiates the dataset from others in the same collection. It is displayed on a page in the cellxgene Portal that also has the collection name. To illustrate, the first dataset name in the [Cells of the adult human heart collection](https://cellxgene.cziscience.com/collections/b52eb423-5d0d-4645-b217-e1c6d38b2e72) is "All — Cells of the adult human heart". | Curator |
||||

​Curators MAY also annotate the following optional keys and values in `uns`. If the key is present, then its value MUST NOT be empty.
​
| Key | Value | Annotator |
:--|:--|:--|
| batch_condition | `str` or `list[str]`. `str` values must refer to cell metadata keys in `obs`. Together, these keys define the "batches" that a normalization or integration algorithm should be aware of. For example if "patient" and "seqBatch" are keys of vectors of cell metadata, either `"patient"`, `"seqBatch"`, or `["patient", "seqBatch"]` are valid values. | Curator |
| default\_embedding|`str`. The `str` value MUST match a key to an embedding in `obsm` for the embedding to display by default. | Curator |
| <obs_column>_colors where <obs_column> MUST be a column name from `obs`. | `list` of  color values in the formats supported by [matplotlib](https://matplotlib.org/stable/tutorials/colors/colors.html). cellxgene Explorer will display [scanpy-style color information](https://github.com/chanzuckerberg/cellxgene/issues/1152#issue-564361541). | Curator |
||||

## Appendix A. Changelog

**EDITOR NOTE**: *Update prior to final commit*
