
# Schema

Contact: brianraymor@chanzuckerberg.com

Document Status: _Drafting_

Version: 7.0.0

The key words "MUST", "MUST NOT", "REQUIRED", "SHALL", "SHALL NOT", "SHOULD", "SHOULD NOT", "RECOMMENDED", "NOT RECOMMENDED" "MAY", and "OPTIONAL" in this document are to be interpreted as described in [BCP 14](https://tools.ietf.org/html/bcp14), [RFC2119](https://www.rfc-editor.org/rfc/rfc2119.txt), and [RFC8174](https://www.rfc-editor.org/rfc/rfc8174.txt) when, and only when, they appear in all capitals, as shown here.

## Schema versioning

The CELLxGENE schema version is based on [Semantic Versioning](https://semver.org/).

**Major version** is incremented when schema updates are incompatible with the AnnData data encoding or CELLxGENE API(s). Examples include:
  * Renaming metadata fields
  * Deprecating metadata fields
  * Changing the type or format of a metadata field
 
**Minor version** is incremented when schema updates may require changes only to the `cellxgene-schema` CLI or the curation process. Examples include:
  * Adding metadata fields
  * Updating pinned ontologies or gene references
  * Changing the validation requirements for a metadata field
  
**Patch version** is incremented for editorial updates and when adding organisms that do not require new metadata fields.

All changes are documented in the schema [Changelog](#appendix-a-changelog).


## Background

CELLxGENE aims to support the publication, sharing, and exploration of single-cell datasets. Building on those published datasets, CELLxGENE seeks to create references of the phenotypes and composition of cells that make up human tissues.

Creating references from multiple datasets requires some harmonization of metadata and features, but if that harmonization is too onerous, it will burden the goal of rapid data sharing. CELLxGENE balances publishing and reference creation needs by requiring datasets hosted by CELLxGENE Discover to include a small set of metadata readily available from data submitters.

This document describes the schema, a type of contract, that CELLxGENE requires all datasets to adhere to so that it can enable searching, filtering, and integration of datasets it hosts.

Note that the requirements in the schema are just the minimum required information. Datasets often have additional metadata, which is preserved in datasets submitted to CELLxGENE Discover.

## Overview

This schema supports multiple assay types. Each assay takes the form of one or more two-dimensional matrices whose values are quantitative measures of the phenotypes of cells.

The schema additionally describes how the dataset, genes, and cells are annotated to describe the biological and technical characteristics of the data.

This document is organized by:

* [General requirements](#general-requirements)
* [`X` (Matrix layers)](#x-matrix-layers), which describe the data required for different assays
* [`obs` (Cell metadata)](#obs-cell-metadata), which describe each cell in the dataset
* [`obsm` (Embeddings)](#obsm-embeddings), which describe each embedding in the dataset
* [`obsp`](#obsp), which describe pairwise annotation of observations
* [`var` and `raw.var` (Gene metadata)](#var-and-rawvar-gene-metadata), which describe each gene in the dataset
* [`varm`](#varm), which describe multi-dimensional annotation of variables/features
* [`varp`](#varp), which describe pairwise annotation of variables/features
* [`uns` (Dataset metadata)](#uns-dataset-metadata), which describe the dataset as a whole

## General Requirements

**AnnData.** The canonical data format for CELLxGENE Discover is HDF5-backed [AnnData](https://anndata.readthedocs.io/en/latest) as written by AnnData version 0.8.0 or greater. The on-disk format must be [AnnData specification (v0.1.0)](https://anndata.readthedocs.io/en/latest/fileformat-prose.html#anndata-specification-v0-1-0). Part of the rationale for selecting this format is to allow CELLxGENE to access both the data and metadata within a single file. The schema requirements and definitions for the AnnData `X`, `obs`, `var`, `raw.var`, `obsm`, and `uns` attributes are described below.

**Reserved Names**. The names of metadata fields MUST NOT start with `"__"`. The names of the metadata fields specified by the schema are **reserved** for the purposes and specifications described in the schema.

**Unique Names**. The names of schema and data submitter metadata fields in `obs` and `var` MUST be unique. For example, duplicate `"feature_biotype"` keys in AnnData `var` are not allowed.

Reserved Names from previous schema versions that have since been deprecated MUST NOT be present in datasets:

<table>
<thead>
  <tr>
    <th>Reserved Name</th>
    <th>AnnData</th>
    <th>Deprecated in</th>
  </tr>
</thead>
<tbody>
 <tr>
    <td>organism</td>
    <td>obs</td>
    <td>6.0.0</td>
  </tr>
  <tr>
    <td>organism_ontology_term_id</td>
    <td>obs</td>
    <td>6.0.0</td>
  </tr>
  <tr>
  <tr>
    <td>ethnicity</td>
    <td>obs</td>
    <td>3.0.0</td>
  </tr>
  <tr>
    <td>ethnicity_ontology_term_id</td>
    <td>obs</td>
    <td>3.0.0</td>
  </tr>
  <tr>
    <td>X_normalization</td>
    <td>uns</td>
    <td>3.0.0</td>
  </tr>
  <tr>
    <td>default_field</td>
    <td>uns</td>
    <td>2.0.0</td>
  </tr>
  <tr>
    <td>layer_descriptions</td>
    <td>uns</td>
    <td>2.0.0</td>
  </tr>
  <tr>
    <td>tags</td>
    <td>uns</td>
    <td>2.0.0</td>
  </tr>
   <tr>
    <td>version</td>
    <td>uns</td>
    <td>2.0.0</td>
  </tr>
  <tr>
    <td>contributors</td>
    <td>uns</td>
    <td>1.1.0</td>
  </tr>
  <tr>
    <td>preprint_doi</td>
    <td>uns</td>
    <td>1.1.0</td>
  </tr>
  <tr>
    <td>project_description</td>
    <td>uns</td>
    <td>1.1.0</td>
  </tr>
  <tr>
    <td>project_links</td>
    <td>uns</td>
    <td>1.1.0</td>
  </tr>
  <tr>
    <td>project_name</td>
    <td>uns</td>
    <td>1.1.0</td>
  </tr>
  <tr>
    <td>publication_doi</td>
    <td>uns</td>
    <td>1.1.0</td>
  </tr>
  <tr>
    <td></td>
    <td></td>
    <td></td>
  </tr>
</tbody>
</table>

**Redundant Metadata**. It is STRONGLY RECOMMENDED to avoid multiple metadata fields containing identical or similar information.

**No Personal Identifiable Information (PII)**.  This is not strictly enforced by validation because it is difficult for software to predict what is and is not PII; however, curators MUST agree to the data submission policies of CELLxGENE Discover on behalf of data submitters which includes this requirement:

> It is my responsibility to ensure that this data is not identifiable. In particular, I commit that I will remove any [direct personal identifiers](https://docs.google.com/document/d/1sboOmbafvMh3VYjK1-3MAUt0I13UUJfkQseq8ANLPl8/edit) in the metadata portions of the data, and that CZI may further contact me if it believes more work is needed to de-identify it.

This includes names, emails, or other PII for researchers or curators involved in the data generation and submission.

#### *Note on types*
The types below are python3 types. Note that a python3 `str` is a sequence of Unicode code points, which is stored null-terminated and UTF-8-encoded by AnnData.

## `X` (Matrix Layers)

The data stored in the `AnnData.X` data matrix is the data that is viewable in CELLxGENE Explorer. For `AnnData.X`, `AnnData.raw.X`, and all layers, if a data matrix contains 50% or more values that are zeros, it MUST be encoded as a [`scipy.sparse.csr_matrix`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html) with zero values encoded as <a href="https://docs.scipy.org/doc/scipy/tutorial/sparse.html#sparse-arrays-implicit-zeros-and-duplicates">implicit zeros</a>.


CELLxGENE's matrix layer requirements are tailored to optimize data reuse. Because each assay has different characteristics, the requirements differ by assay type. In general, CELLxGENE requires submission of "raw" data suitable for computational reuse when a standard raw matrix format exists for an assay. It is STRONGLY RECOMMENDED to also include a "normalized" matrix with processed values ready for data analysis and suitable for visualization in CELLxGENE Explorer. So that CELLxGENE's data can be provided in download formats suitable for both R and Python, the schema imposes the following requirements:

*   All matrix layers MUST have the same shape, and have the same cell labels and gene labels.
*   Because it is impractical to retain all barcodes in raw and normalized matrices, any cell filtering MUST be applied to both.
    By contrast, those wishing to reuse datasets require access to raw gene expression values, so genes SHOULD NOT be filtered from either dataset.
    Summarizing, any cell barcodes that are removed from the data MUST be filtered from both raw and normalized matrices and genes SHOULD NOT be filtered from the raw matrix.
*   Any genes that publishers wish to filter from the normalized matrix MAY have their values replaced by zeros and MUST be flagged in the column [`feature_is_filtered`](#feature_is_filtered) of [`var`](#var-and-rawvar-gene-metadata), which will mask them from exploration.
*   Additional layers provided at author discretion MAY be stored using author-selected keys, but MUST have the same cells and genes as other layers. It is STRONGLY RECOMMENDED that these layers have names that accurately summarize what the numbers in the layer represent (e.g. `"counts_per_million"`, `"SCTransform_normalized"`, or `"RNA_velocity_unspliced"`).

### Definitions for scATAC-seq assays

<b>paired assay</b>. `assay_ontology_term_id` is a descendant of both <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010891"><code>"EFO:0010891"</code></a> for <i>scATAC-seq</i> and <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008913"><code>"EFO:0008913"</code></a> for <i>single-cell RNA sequencing</i>. A gene expression matrix (RNA data) is required.

<b>unpaired assay</b>. `assay_ontology_term_id` is <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010891"><code>"EFO:0010891"</code></a> for <i>scATAC-seq</i> or a descendant and is not a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008913"><code>"EFO:0008913"</code></a> for <i>single-cell RNA sequencing</i>. A gene activity matrix and not a peak matrix is required.

Also see the requirements for [scATAC-seq assets](#scatac-seq-assets).<br><br>

The following table describes the matrix data and layers requirements that are **assay-specific**. If an entry in the table is empty, the schema does not have any other requirements on data in those layers beyond the ones listed above.

| Assay | "raw" required? | "raw" location | "normalized" required? | "normalized" location |
|-|-|-|-|-|
| scRNA-seq (UMI, e.g. 10x multiome, 10x v3, Slide-seqV2) | REQUIRED. Values MUST be de-duplicated molecule counts. Each cell MUST contain at least one non-zero value. All non-zero values MUST be positive integers stored as `numpy.float32`. Any two cells MUST NOT contain identical values for all their features. | `AnnData.raw.X` unless no "normalized" is provided, then `AnnData.X` | STRONGLY RECOMMENDED | `AnnData.X` |
| Visium Spatial (e.g. V1, CytAssist) | REQUIRED. Values MUST be de-duplicated molecule counts. All non-zero values MUST be positive integers stored as `numpy.float32`.<br><br>If <code>uns['spatial']['is_single']</code> is <code>False</code> then each cell MUST contain at least one non-zero value.<br><br>If <code>uns['spatial']['is_single']</code> is <code>True</code> then the unfiltered feature-barcode matrix (<code>raw_feature_bc_matrix</code>) MUST be used. See <a href="https://www.10xgenomics.com/support/software/space-ranger/analysis/outputs/space-ranger-feature-barcode-matrices">Space Ranger Feature-Barcode Matrices</a>.<br><br>if <code>assay_ontology_term_id</code> is <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0022860"><code>"EFO:0022860"</code></a> for <i>Visium CytAssist Spatial Gene Expression, 11mm</i>, this matrix MUST contain 14336 rows; otherwise, this matrix MUST contain 4992 rows.<br><br>If the <code>obs['in_tissue']</code> value is <code>1</code>, then the cell MUST contain at least one non-zero value and any two cells MUST NOT contain identical values for all their features.<br><br>If any <code>obs['in_tissue']</code> values are <code>0</code>, then at least one cell corresponding to a <code>obs['in_tissue']</code> with a value of <code>0</code> MUST contain a non-zero value.| `AnnData.raw.X` unless no "normalized" is provided, then `AnnData.X` | STRONGLY RECOMMENDED | `AnnData.X` |
| scRNA-seq (non-UMI, e.g. SS2) | REQUIRED. Values MUST be one of read counts (e.g. FeatureCounts) or  estimated fragments (e.g. output of RSEM). Each cell MUST contain at least one non-zero value. All non-zero values MUST be positive integers stored as `numpy.float32`. Any two cells MUST NOT contain identical values for all their features. | `AnnData.raw.X` unless no "normalized" is provided, then `AnnData.X` | STRONGLY RECOMMENDED | `AnnData.X` |
| unpaired Accessibility (e.g. ATAC-seq, mCT-seq) | NOT REQUIRED | | REQUIRED | `AnnData.X` | STRONGLY RECOMMENDED |
|||||

## Integration Metadata

CELLxGENE requires ontology terms to enable search, comparison, and integration of data. With the exception of Cellosaurus, ontology terms for cell metadata MUST use [OBO-format identifiers](http://www.obofoundry.org/id-policy.html), meaning a CURIE (prefixed identifier) of the form **Ontology:Identifier**. For example, [EFO:0000001](https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0000001) is a term in the Experimental Factor Ontology (EFO). Cellosaurus requires a prefixed identifier of the form **Ontology_Identifier** such as [CVCL_1P02](https://www.cellosaurus.org/CVCL_1P02).


The most accurate ontology term MUST always be used. If an exact or approximate ontology term is not available, a new term may be requested:

- For the [Cell Ontology], data submitters may [suggest a new term](https://github.com/obophenotype/cell-ontology/issues/new?assignees=bvarner-ebi&labels=new+term+request%2C+cellxgene&template=a_adding_term_cellxgene.md&title=%5BNTR-cxg%5D) and [notify the curation team](mailto:cellxgene@chanzuckerberg.com) of the pending term request, so that the datasets can be updated once the term is available.

  To meet CELLxGENE schema requirements, the most accurate available CL term MUST be used until the new term is available. For example if `cell_type_ontology_term_id` describes a relay interneuron, but the most accurate available term in the CL ontology is [CL:0000099](https://www.ebi.ac.uk/ols4/ontologies/cl/classes?obo_id=CL%3A0000099) for *interneuron*, then the interneuron term can be used to fulfill this requirement and ensures that users searching for "neuron" are able to find these data.  If no appropriate term can be found (e.g. the cell type is unknown), then `"unknown"` MUST be used. Users will still be able to access more specific cell type annotations that have been submitted with the dataset (but aren't required by the schema).

   
- For all other ontologies, data submitters may submit a [request to the curation team](mailto:cellxgene@chanzuckerberg.com) during the submission process.

Terms documented as obsolete in an ontology MUST NOT be used. For example, [EFO:0009310](https://www.ebi.ac.uk/ols4/ontologies/efo/classes/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0009310) for *obsolete_10x v2* was marked as obsolete in EFO version 3.31.0 and replaced by [EFO:0009899](https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0009899) for *10x 3' v2*.

### Required Ontologies

The following ontology dependencies are *pinned* for this version of the schema.

| Ontology | Prefix | Release | Download |
|:--|:--|:--|:--|
| [C. elegans Development Ontology] | WBls: |  [2025-04-01 WS297](https://github.com/obophenotype/c-elegans-development-ontology/releases/tag/v2025-04-01) | [wbls.owl](https://github.com/obophenotype/c-elegans-development-ontology/blob/v2025-04-01/wbls.owl) |
| [C. elegans Gross Anatomy Ontology] | WBbt: | [2025-03-26 WS297](https://github.com/obophenotype/c-elegans-gross-anatomy-ontology/releases/tag/v2025-03-26) | [wbbt.owl](https://github.com/obophenotype/c-elegans-gross-anatomy-ontology/blob/v2025-03-26/wbbt.owl) |
| [Cell Ontology] | CL: |  [2025-04-10](https://github.com/obophenotype/cell-ontology/releases/tag/v2025-04-10) | [cl.owl](https://github.com/obophenotype/cell-ontology/releases/download/v2025-04-10/cl.owl)|
| [Cellosaurus] | CVCL_ | 52.0 | [cellosaurus.obo ](https://ftp.expasy.org/databases/cellosaurus/cellosaurus.obo)_(Versioned releases are unavailable. Cellosaurus may replace this download with a newer release.)_ |
| [Drosophila Anatomy Ontology] | FBbt: | [2025-03-27](https://github.com/FlyBase/drosophila-anatomy-developmental-ontology/releases/tag/v2025-03-27)| [fbbt.owl](https://github.com/FlyBase/drosophila-anatomy-developmental-ontology/releases/download/v2025-03-27/fbbt.owl) |
| [Drosophila Development Ontology] | FBdv: | [2025-03-26](https://github.com/FlyBase/drosophila-developmental-ontology/releases/tag/v2025-03-26) | [fbdv.owl](https://github.com/FlyBase/drosophila-developmental-ontology/releases/download/v2025-03-26/fbdv.owl) |
| [Experimental Factor Ontology] | EFO: | [2025-05-15 EFO 3.78.0](https://github.com/EBISPOT/efo/releases/tag/v3.78.0) | [efo.owl](https://github.com/EBISPOT/efo/releases/download/v3.78.0/efo.owl) |
| [Human Ancestry Ontology] | AfPO:<br>HANCESTRO: | [2025-04-01](https://github.com/EBISPOT/hancestro/releases/tag/v2025-04-01) | [hancestro-base.owl](https://github.com/EBISPOT/hancestro/blob/v2025-04-01/hancestro-base.owl) |
| [Human Developmental Stages] |  HsapDv: | [2025-01-23](https://github.com/obophenotype/developmental-stage-ontologies/releases/tag/v2025-01-23) | [hsapdv.owl](https://github.com/obophenotype/developmental-stage-ontologies/releases/download/v2025-01-23/hsapdv.owl) |
| [Mondo Disease Ontology] | MONDO: | [2025-05-06](https://github.com/monarch-initiative/mondo/releases/tag/v2025-05-06) | [mondo.owl](https://github.com/monarch-initiative/mondo/releases/download/v2025-05-06/mondo.owl) |
| [Mouse Developmental Stages]| MmusDv: | [2025-01-23](https://github.com/obophenotype/developmental-stage-ontologies/releases/tag/v2025-01-23) | [mmusdv.owl](https://github.com/obophenotype/developmental-stage-ontologies/releases/download/v2025-01-23/mmusdv.owl) |
| [NCBI organismal classification] |  NCBITaxon: | [2025-03-13](https://github.com/obophenotype/ncbitaxon/releases/tag/v2025-03-13) | [ncbitaxon.owl](https://github.com/obophenotype/ncbitaxon/releases/download/v2025-03-13/ncbitaxon.owl.gz) |
| [Phenotype And Trait Ontology] | PATO: | [2025-05-14](https://github.com/pato-ontology/pato/releases/tag/v2025-05-14) | [pato.owl](https://github.com/pato-ontology/pato/blob/v2025-05-14/pato.owl)  |
| [Uberon multi-species anatomy ontology] |  UBERON: | [2025-05-28](https://github.com/obophenotype/uberon/releases/tag/v2025-05-28) | [uberon.owl](https://github.com/obophenotype/uberon/releases/download/v2025-05-28/uberon.owl) |
| [Zebrafish Anatomy Ontology] | ZFA:<br>ZFS: | [2025-01-28](https://github.com/ZFIN/zebrafish-anatomical-ontology/releases/tag/v2025-01-28) | [zfa.owl](https://github.com/ZFIN/zebrafish-anatomical-ontology/blob/v2025-01-28/zfa.owl) |
| | | | |

[C. elegans Development Ontology]: https://obofoundry.org/ontology/wbls.html

[C. elegans Gross Anatomy Ontology]: https://obofoundry.org/ontology/wbbt.html

[Cell Ontology]: https://obofoundry.org/ontology/cl.html

[Cellosaurus]: https://www.cellosaurus.org/description.html

[Drosophila Anatomy Ontology]: https://obofoundry.org/ontology/fbbt.html

[Drosophila Development Ontology]: https://obofoundry.org/ontology/fbdv.html

[Experimental Factor Ontology]: https://www.ebi.ac.uk/efo

[Human Ancestry Ontology]: https://www.obofoundry.org/ontology/hancestro.html

[Human Developmental Stages]: https://obofoundry.org/ontology/hsapdv.html

[Mondo Disease Ontology]: https://obofoundry.org/ontology/mondo.html

[Mouse Developmental Stages]: https://obofoundry.org/ontology/mmusdv.html

[NCBI organismal classification]: https://obofoundry.org/ontology/ncbitaxon.html

[Phenotype And Trait Ontology]: https://www.obofoundry.org/ontology/pato.html

[Uberon multi-species anatomy ontology]: https://www.obofoundry.org/ontology/uberon.html

[Zebrafish Anatomy Ontology]: https://obofoundry.org/ontology/zfa.html


### Required Gene Annotations

ENSEMBL identifiers are required for genes and [External RNA Controls Consortium (ERCC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4978944/) identifiers for [RNA Spike-In Control Mixes] to ensure that all datasets measure the same features and can therefore be integrated.

The following gene annotation dependencies are *pinned* for this version of the schema.

| NCBITaxon | Source | Required version | Download |
|:--|:--|:--|:--|
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A6239"><code>NCBITaxon:6239</code></a><br>for <i>Caenorhabditis elegans</i> | [ENSEMBL](https://www.ensembl.org/Caenorhabditis_elegans/Info/Index) | WBcel235<br>(GCA_000002985.3)<br>Ensembl 114 | [Caenorhabditis_elegans.WBcel235.114.gtf.gz](https://ftp.ensembl.org/pub/release-114/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.114.gtf.gz) |
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9483"><code>NCBITaxon:9483</code></a><br>for <i>Callithrix jacchus</i>  | [ENSEMBL](https://www.ensembl.org/Callithrix_jacchus/Info/Index) | mCalJac1.pat.X<br>(GCA_011100555.1)<br>Ensembl 114 | [Callithrix_jacchus.mCalJac1.pat.X.114.gtf.gz](https://ftp.ensembl.org/pub/release-114/gtf/callithrix_jacchus/Callithrix_jacchus.mCalJac1.pat.X.114.gtf.gz) |
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A7955"><code>NCBITaxon:7955</code></a><br>for <i>Danio rerio</i> |  [ENSEMBL](https://www.ensembl.org/Danio_rerio/Info/Index) | GRCz11<br>(GCA_000002035.4)<br>Ensembl 114 | [Danio_rerio.GRCz11.114.gtf.gz](https://ftp.ensembl.org/pub/release-114/gtf/danio_rerio/Danio_rerio.GRCz11.114.gtf.gz) |
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A7227"><code>NCBITaxon:7227</code></a><br>for <i>Drosophila melanogaster</i>| [ENSEMBL](https://www.ensembl.org/Drosophila_melanogaster/Info/Index) | BDGP6.54<br>(GCA_000001215.4)<br>Ensembl 114 | [Drosophila_melanogaster.BDGP6.54.114.gtf.gz](https://ftp.ensembl.org/pub/release-114/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.114.gtf.gz) |
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9595"><code>NCBITaxon:9595</code></a><br>for <i>Gorilla gorilla gorilla</i>  | [ENSEMBL](https://www.ensembl.org/Gorilla_gorilla/Info/Index) | gorGor4<br>(GCA_000151905.3)<br>Ensembl 114 | [Gorilla_gorilla.gorGor4.114.gtf.gz](https://ftp.ensembl.org/pub/release-114/gtf/gorilla_gorilla/Gorilla_gorilla.gorGor4.114.gtf.gz) |
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9606"><code>NCBITaxon:9606</code></a><br>for <i>Homo sapiens</i> | [GENCODE](https://www.gencodegenes.org/human/) | GENCODE v48<br>(GRCh38.p14)<br>Ensembl 114 | [gencode.v48.primary_assembly.annotation.gtf](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz) |
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9541"><code>NCBITaxon:9541</code></a><br>for <i>Macaca fascicularis</i>  | [ENSEMBL](https://www.ensembl.org/Macaca_fascicularis/Info/Index) | Macaca_fascicularis_6.0<br>(GCA_011100615.1)<br>Ensembl 114 | [Macaca_fascicularis.Macaca_fascicularis_6.0.114.gtf.gz](https://ftp.ensembl.org/pub/release-114/gtf/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_6.0.114.gtf.gz) |
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9544"><code>NCBITaxon:9544</code></a><br>for <i>Macaca mulatta</i>  | [ENSEMBL](https://www.ensembl.org/Macaca_mulatta/Info/Index) | Mmul_10<br>(GCA_003339765.3)<br>Ensembl 114 | [Macaca_mulatta.Mmul_10.114.gtf.gz](https://ftp.ensembl.org/pub/release-114/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.114.gtf.gz) |
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A30608"><code>NCBITaxon:30608</code></a><br>for <i>Microcebus murinus</i>  | [ENSEMBL](https://www.ensembl.org/Microcebus_murinus/Info/Index) | Mmur_3.0<br>(GCA_000165445.3)<br>Ensembl 114| [Microcebus_murinus.Mmur_3.0.114.gtf.gz](https://ftp.ensembl.org/pub/release-114/gtf/microcebus_murinus/Microcebus_murinus.Mmur_3.0.114.gtf.gz) |
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A10090"><code>NCBITaxon:10090</code></a><br>for <i>Mus musculus</i> | [GENCODE](https://www.gencodegenes.org/mouse/) |  GENCODE vM37<br>(GRCm39)<br>Ensembl 114 | [gencode.vM37.primary_assembly.annotation.gtf](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/gencode.vM37.primary_assembly.annotation.gtf.gz) |
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9986"><code>NCBITaxon:9986</code></a><br>for <i>Oryctolagus cuniculus</i>  | [ENSEMBL](https://www.ensembl.org/Oryctolagus_cuniculus/Info/Index) | OryCun2.0<br>(GCA_000003625.1)<br>Ensembl 114 | [Oryctolagus_cuniculus.OryCun2.0.114.gtf.gz](https://ftp.ensembl.org/pub/release-114/gtf/oryctolagus_cuniculus/Oryctolagus_cuniculus.OryCun2.0.114.gtf.gz) |
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9598"><code>NCBITaxon:9598</code></a><br>for <i>Pan troglodytes</i>  | [ENSEMBL](https://www.ensembl.org/Pan_troglodytes/Info/Index) | Pan_tro_3.0<br>(GCA_000001515.5)<br>Ensembl 114 | [Pan_troglodytes.Pan_tro_3.0.114.gtf.gz](https://ftp.ensembl.org/pub/release-114/gtf/pan_troglodytes/Pan_troglodytes.Pan_tro_3.0.114.gtf.gz) |
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A10116"><code>NCBITaxon:10116</code></a><br>for <i>Rattus norvegicus</i>  | [ENSEMBL](https://www.ensembl.org/Rattus_norvegicus/Info/Index) | GRCr8<br>(GCA_036323735.1)<br>Ensembl 114| [Rattus_norvegicus.GRCr8.114.gtf.gz](https://ftp.ensembl.org/pub/release-114/gtf/rattus_norvegicus/Rattus_norvegicus.GRCr8.114.gtf.gz)  |
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A2697049"><code>NCBITaxon:2697049</code></a><br>for <i>SARS-CoV-2</i>  | [ENSEMBL](https://covid-19.ensembl.org/index.html) | SARS-CoV-2 reference (ASM985889v3) | [Sars\_cov\_2.ASM985889v3.101.gtf](https://ftp.ensemblgenomes.org/pub/viruses/gtf/sars_cov_2/Sars_cov_2.ASM985889v3.101.gtf.gz) |
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9823"><code>NCBITaxon:9823</code></a><br>for <i>Sus scrofa</i> | [ENSEMBL](https://www.ensembl.org/Sus_scrofa/Info/Index) |  Sscrofa11.1<br>(GCA_000003025.6)<br>Ensembl 114 | [Sus_scrofa.Sscrofa11.1.114.gtf.gz](https://ftp.ensembl.org/pub/release-114/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.114.gtf.gz) |
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A32630"><code>NCBITaxon:32630</code></a><br>for <i>synthetic construct</i> | [ThermoFisher ERCC<br>Spike-Ins] | ThermoFisher ERCC RNA Spike-In Control Mixes (Cat # 4456740, 4456739) | [cms_095047.txt] |
|||||

[RNA Spike-In Control Mixes]: https://www.thermofisher.com/document-connect/document-connect.html?url=https%3A%2F%2Fassets.thermofisher.com%2FTFS-Assets%2FLSG%2Fmanuals%2Fcms_086340.pdf&title=VXNlciBHdWlkZTogRVJDQyBSTkEgU3Bpa2UtSW4gQ29udHJvbCBNaXhlcyAoRW5nbGlzaCAp

[ThermoFisher ERCC<br>Spike-Ins]: https://www.thermofisher.com/order/catalog/product/4456740#/4456740
[cms_095047.txt]: https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095047.txt


## `obs` (Cell Metadata)

`obs` is a [`pandas.DataFrame`](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html).

Curators MUST annotate the following columns in the `obs` dataframe:

### index of pandas.DataFrame

<table><tbody>
    <tr>
      <th>Key</th>
      <td>index of <code>pandas.DataFrame</code></td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>str</code>. The index of the pandas.DataFrame MUST contain unique identifiers for observations.<br><br></td>
    </tr>
</tbody></table>
<br>

### array_col

<table><tbody>
    <tr>
      <th>Key</th>
      <td>array_col</td>
    </tr>
    <tr>
      <th>Annotator</th>
         <td>Curator MUST annotate if <code>assay_ontology_term_id</code> is  a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> and <code>uns['spatial']['is_single']</code> is <code>True</code>; otherwise, this key MUST NOT be present.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>int</code>. This MUST be the value of the column coordinate for the corresponding spot from the <code>array_col</code> field in <code>tissue_positions_list.csv</code> or <code>tissue_positions.csv</code>. See <a href="https://www.10xgenomics.com/support/software/space-ranger/analysis/outputs/spatial-outputs">Space Ranger Spatial Outputs</a>.<br>
        <br>
          <table>
          <thead>
          <tr>
          <th>For Assay</th>
          <th>Value MUST be in<br>the range between</th>
          </tr>
          </thead>
          <tbody>
            <tr>
              <td><i>Visium Spatial Gene Expression V1</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0022857"><code>EFO:0022857</code></a>]</td>
              <td><code>0</code> and <code>127</code></td>
           </tr> 
            <tr>
              <td><i>Visium CytAssist Spatial Gene Expression, 6.5mm</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0022859"><code>EFO:0022859</code></a>]</td>
              <td><code>0</code> and <code>127</code></td>
           </tr>
            <tr>
              <td><i>Visium CytAssist Spatial Gene Expression, 11mm</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0022860"><code>EFO:0022860</code></a>]</td>
              <td><code>0</code> and <code>223</code></td>
           </tr>
          </tbody></table>
        </td>
    </tr>
</tbody></table>
<br>

### array_row

<table><tbody>
    <tr>
      <th>Key</th>
      <td>array_row</td>
    </tr>
    <tr>
      <th>Annotator</th>
         <td>Curator MUST annotate if <code>assay_ontology_term_id</code> is a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> and <code>uns['spatial']['is_single']</code> is <code>True</code>; otherwise, this key MUST NOT be present.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>int</code>. This MUST be value of the row coordinate for the corresponding spot from the <code>array_row</code> field in in <code>tissue_positions_list.csv</code> or <code>tissue_positions.csv</code>. See <a href="https://www.10xgenomics.com/support/software/space-ranger/analysis/outputs/spatial-outputs">Space Ranger Spatial Outputs</a>.<br>
        <br>
          <table>
          <thead>
          <tr>
          <th>For Assay</th>
          <th>Value MUST be in<br>the range between</th>
          </tr>
          </thead>
          <tbody>
            <tr>
              <td><i>Visium Spatial Gene Expression V1</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0022857"><code>EFO:0022857</code></a>]</td>
              <td><code>0</code> and <code>77</code></td>
           </tr> 
            <tr>
              <td><i>Visium CytAssist Spatial Gene Expression, 6.5mm</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0022859"><code>EFO:0022859</code></a>]</td>
              <td><code>0</code> and <code>77</code></td>
           </tr>
            <tr>
              <td><i>Visium CytAssist Spatial Gene Expression, 11mm</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0022860"><code>EFO:0022860</code></a>]</td>
              <td><code>0</code> and <code>127</code></td>
           </tr>
          </tbody></table>
        </td>
    </tr>
</tbody></table>
<br>

### assay_ontology_term_id

<table><tbody>
    <tr>
      <th>Key</th>
      <td>assay_ontology_term_id</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be an EFO term and either:<br><br>
          <ul><li>
            the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0002772"><code>"EFO:0002772"</code></a> for <i>assay by molecule</i> excluding <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> while allowing its descendants
          </li>
          <li>
            the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010183"><code>"EFO:0010183"</code></a>  for <i>single cell library construction</i> excluding <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> while allowing its descendants
          </li></ul>
        If <code>assay_ontology_term_id</code> is either a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> or <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030062"><code>"EFO:0030062"</code></a> for <i>Slide-seqV2</i> then all observations MUST contain the same value.<br><br>
        If <code>assay_ontology_term_id</code> is either <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010891"><code>"EFO:0010891"</code></a> for <i>scATAC-seq</i> or its descendants, there are additional requirements for separate fragments file assets documented in <a href="#scatac-seq-assets">scATAC-seq assets</a>.<br><br>
        An assay based on 10X Genomics products SHOULD be the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008995"><code>"EFO:0008995"</code></a> for <i>10x technology</i>. An assay based on <i>SMART (Switching Mechanism at the 5' end of the RNA Template) or SMARTer technology</i> SHOULD either be <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010184"><code>"EFO:0010184"</code></a> for <i>Smart-like</i> or preferably its most accurate descendant.<br><br>
       <br>Recommended values for specific assays:
          <br><br>
          <table>
          <thead>
          <tr>
          <th>For</th>
          <th>Use</th>
          </tr>
          </thead>
          <tbody>
            <tr>
              <td><i>10x 3' v2</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0009899"><code>"EFO:0009899"</code></a></td>
            </tr>
            <tr>
              <td><i>10x 3' v3</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0009922"><code>"EFO:0009922"</code></a></td>
            </tr>
            <tr>
              <td><i>10x 3' v4</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0022604"><code>"EFO:0022604"</code></a></td>
            </tr>
            <tr>
              <td><i>10x 5' v1</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0011025"><code>"EFO:0011025"</code></a></td>
            </tr>
            <tr>
              <td><i>10x 5' v2</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0009900"><code>"EFO:0009900"</code></a></td>
            </tr>            <tr>
              <td><i>10x 5' v3</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0022605"><code>"EFO:0022605"</code></a></td>
            </tr>            </tr>            <tr>
              <td><i>10x multiome</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030059"><code>"EFO:0030059"</code></a></td>
            </tr>
            <tr>
              <td><i>Smart-seq2</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008931"><code>"EFO:0008931"</code></a></td>
            </tr>
            <tr>
              <td><i>Visium Spatial Gene Expression V1</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0022857"><code>"EFO:0022857"</code></a></td>
            </tr>
            <tr>
              <td><i>Visium CytAssist Spatial Gene Expression, 6.5mm</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0022859"><code>"EFO:0022859"</code></a></td>
            </tr>
            <tr>
              <td><i>Visium CytAssist Spatial Gene Expression, 11mm</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0022860"><code>"EFO:0022860"</code></a></td>
            </tr>
          </tbody></table>
        </td>
    </tr>
</tbody></table>
<br>

### cell_type_ontology_term_id

<table><tbody>
  <tr>
    <th>Key</th>
    <td>cell_type_ontology_term_id</td>
  </tr>
  <tr>
    <th>Annotator</th>
    <td>Curator MUST annotate.</td>
  </tr>
  <tr>
    <th>Value</th>
    <td>
      categorical with <code>str</code> categories.<br><br>This MUST be <code>"unknown"</code> when:
      <ul>
        <li>
          no appropriate term can be found (e.g. the cell type is unknown)
        </li>
        <li>
          <code>assay_ontology_term_id</code> is a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i>,<br><code>uns['spatial']['is_single']</code> is <code>True</code>,<br>and the corresponding value of <code>in_tissue</code> is <code>0</code>
        </li>
      </ul>
        <br>If <code>tissue_type</code> is <code>"cell line"</code>, this MAY be <code>"na"</code>, but then all observations where <code>tissue_type</code> is <code>"cell line"</code> MUST be <code>"na"</code>.<br><br>The following CL terms MUST NOT be used:
        <ul><li>
          <a href="https://www.ebi.ac.uk/ols4/ontologies/cl/terms?obo_id=CL:0000255"><code>"CL:0000255"</code></a> for <i>eukaryotic cell</i>
        </li>
        <li>
          <a href="https://www.ebi.ac.uk/ols4/ontologies/cl/terms?obo_id=CL:0000257"><code>"CL:0000257"</code></a> for <i>Eumycetozoan cell</i>
        </li>
        <li>
            <a href="https://www.ebi.ac.uk/ols4/ontologies/cl/terms?obo_id=CL:0000548"><code>"CL:0000548"</code></a> for <i>animal cell</i>
         </li></ul><br>
      <table>
        <thead><tr>
          <th>For <code>organism_ontology_term_id</code></th>
          <th>Value</th>
        </tr></thead>
        <tbody>
          <tr>
            <td>
              <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A6239"><code>"NCBITaxon:6239"</code></a><br>for <i>Caenorhabditis elegans</i>
            </td>
            <td>
              MUST be either a CL term or the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/wbbt/classes?obo_id=WBbt%3A0004017"><code>WBbt:0004017</code></a><br>for <i>Cell</i> excluding <a href="https://www.ebi.ac.uk/ols4/ontologies/wbbt/classes?obo_id=WBbt%3A0006803"><code>WBbt:0006803</code></a> for <i>Nucleus</i> and its descendants
            </td>
          </tr>
          <tr>
            <td>
              <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A7955"><code>"NCBITaxon:7955"</code></a><br>for <i>Danio rerio</i>
            </td>
            <td>
              MUST be either a CL term or the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/zfa/classes?obo_id=ZFA%3A0009000"><code>ZFA:0009000</code></a> <br>for <i>cell</i>
            </td>
          </tr>
          <tr>
            <td>
              <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A7227"><code>"NCBITaxon:7227"</code></a><br>for <i>Drosophila melanogaster</i>
            </td>
            <td>
              MUST be either a CL term or the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/fbbt/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FFBbt_00007002?lang=en"><code>FBbt:00007002</code></a><br>for <i>cell</i>
            </td>
          </tr>
          <tr>
            <td>
              For all other organisms
            </td>
            <td>
              MUST be a CL term
            </td>
          </tr>
        </tbody>
      </table>
    </td>
  </tr>
</tbody></table>
<br>

### development_stage_ontology_term_id

<table><tbody>
    <tr>
      <th>Key</th>
      <td>development_stage_ontology_term_id</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
      <td>
        categorical with <code>str</code> categories.<br><br>If <code>tissue_type</code> is <code>"cell line"</code>, this MUST be <code>"na"</code>.<br><br>If unavailable, this MUST be <code>"unknown"</code>.<br><br>
        <table>
          <thead>
            <tr>
              <th>For <code>organism_ontology_term_id</code></th>
              <th>Value</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td>
                <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A6239"><code>"NCBITaxon:6239"</code></a><br>for <i>Caenorhabditis elegans</i>
              </td>
              <td>
                MUST be <a href="https://www.ebi.ac.uk/ols4/ontologies/wbls/classes?obo_id=WBls%3A0000669"><code>WBls:0000669</code></a> for <i>unfertilized egg Ce</i>,<br>the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/wbls/classes?obo_id=WBls%3A0000803"><code>WBls:0000803</code></a><br>for <i>C. elegans life stage occurring during embryogenesis</i>, or<br>the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/wbls/classes?obo_id=WBls%3A0000804"><code>WBls:0000804</code></a><br>for <i>C. elegans life stage occurring post embryogenesis</i> 
              </td>
            </tr>
            <tr>
              <td>
                <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A7955"><code>"NCBITaxon:7955"</code></a><br>for <i>Danio rerio</i>
              </td>
              <td>
                MUST be the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/zfs/classes?obo_id=ZFS%3A0100000"><code>ZFS:0100000</code></a> for <i>zebrafish stage</i><br>excluding <a href="https://www.ebi.ac.uk/ols4/ontologies/zfs/classes?obo_id=ZFS%3A0000000"><code>ZFS:0000000</code></a> for <i>Unknown</i>
              </td>
            </tr>
            <tr>
              <td>
                <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A7227"><code>"NCBITaxon:7227"</code></a><br>for <i>Drosophila melanogaster</i>
              </td>
              <td>
                MUST be either the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/fbdv/classes?obo_id=FBdv%3A00007014"><code>FBdv:00007014</code></a> for<br><i>adult age in days</i> or the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/fbdv/classes?obo_id=FBdv%3A00005259"><code>FBdv:00005259</code></a> for<br><i>developmental stage</i> excluding <a href="https://www.ebi.ac.uk/ols4/ontologies/fbdv/classes?obo_id=FBdv%3A00007012"><code>FBdv:00007012</code></a> for <i>life stage</i>
              </td>
            </tr>
            <tr>
              <td>
                <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9606"><code>"NCBITaxon:9606"</code></a><br>for <i>Homo sapiens</i>
              </td>
              <td>
                MUST be the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/hsapdv/classes?obo_id=HsapDv%3A0000001"><code>HsapDv:0000001</code></a> for <i>life cycle</i>
              </td>
            </tr>
            <tr>
              <td>
                <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A10090"><code>"NCBITaxon:10090"</code></a><br>for <i>Mus musculus</i> or one of its descendants
              </td>
              <td>
                MUST be the accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/mmusdv/classes?obo_id=MmusDv%3A0000001"><code>MmusDv:0000001</code></a> for <i>life cycle</i>
              </td>
            </tr>
            <tr>
              <td>
                For all other organisms
              </td>
              <td>
                MUST be the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/uberon/classes?obo_id=UBERON%3A0000105"><code>UBERON:0000105</code></a> for <i>life cycle stage</i>, excluding <a href="https://www.ebi.ac.uk/ols4/ontologies/uberon/classes?obo_id=UBERON%3A0000071"><code>UBERON:0000071</code></a> for <i>death stage</i>.
              </td>
            </tr>
          </tbody>
        </table>
      </td>
  </tr>
</tbody></table>
<br>

### disease_ontology_term_id

<table><tbody>
    <tr>
      <th>Key</th>
      <td>disease_ontology_term_id</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be one of:<br><br>
        <ul>
          <li><a href="https://www.ebi.ac.uk/ols4/ontologies/pato/classes?obo_id=PATO%3A0000461"><code>"PATO:0000461"</code></a> for <i>normal</i> or <i>healthy</i>.</li>
          <li>one or more MONDO terms in ascending lexical order separated by the delimiter <code>" || "</code> with no duplication of terms. For example, if the terms are <code>"MONDO:1030008"</code>, <code>"MONDO:0800349"</code>, <code>"MONDO:0004604"</code>, and <code>"MONDO:0043004"</code> then the value MUST be <code>"MONDO:0004604 || MONDO:0043004 || MONDO:0800349 || MONDO:1030008"</code>.</li>
       </ul>
       MONDO terms MUST be either:
       <ul>
          <li>a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/mondo/classes?obo_id=MONDO%3A0000001"><code>"MONDO:0000001"</code></a> for <i>disease</i></li>
          <li><a href="https://www.ebi.ac.uk/ols4/ontologies/mondo/classes?obo_id=MONDO%3A0021178"><code>"MONDO:0021178"</code></a> for <i>injury</i> or <b>preferably</b> its most accurate descendant</li>       </ul>
        </td>
    </tr>
</tbody></table>
<br>

### donor_id

<table><tbody>
    <tr>
      <th>Key</th>
      <td>donor_id</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories.<br><br>If <code>tissue_type</code> is <code>"cell line"</code>, this MUST be <code>"na"</code>; otherwise, this MUST NOT be <code>"na"</code>, but MUST be free-text that identifies a unique individual that data were derived from. It is STRONGLY RECOMMENDED that this identifier be designed so that it is unique to:<br><br>
          <ul><li>a given individual within the collection of datasets that includes this dataset</li>
          <li>a given individual across all collections in CELLxGENE Discover</li></ul><br>
          It is STRONGLY RECOMMENDED that <code>"pooled"</code> be used  for observations from a sample of multiple individuals that were not confidently assigned to a single individual through demultiplexing.<br><br>It is STRONGLY RECOMMENDED that <code>"unknown"</code> ONLY be used for observations in a dataset when it is not known which observations are from the same individual.<br><br>
        </td>
    </tr>
</tbody></table>
<br>


### in_tissue

<table><tbody>
    <tr>
      <th>Key</th>
      <td>in_tissue</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate if <code>assay_ontology_term_id</code> is a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> and <code>uns['spatial']['is_single']</code> is <code>True</code>; otherwise, this key MUST NOT be present.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>int</code>. This MUST be the value for the corresponding spot from the <code>in_tissue</code> field in <code>tissue_positions_list.csv</code> or <code>tissue_positions.csv</code> which is either <code>0</code> if the spot falls outside tissue or <code>1</code> if the spot falls inside tissue. See <a href="https://www.10xgenomics.com/support/software/space-ranger/analysis/outputs/spatial-outputs">Space Ranger Spatial Outputs</a>.
        </td>
    </tr>
</tbody></table>
<br>


### is_primary_data

<table><tbody>
    <tr>
      <th>Key</th>
      <td>is_primary_data</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>bool</code>. This MUST be <code>False</code> if <code>uns['spatial']['is_single']</code> is <code>False</code>. This MUST be <code>True</code> if this is the canonical instance of this cellular observation and <code>False</code> if not. This is commonly <code>False</code> for meta-analyses reusing data or for secondary views of data.
        </td>
    </tr>
</tbody></table>
<br>

### self_reported_ethnicity_ontology_term_id

<table>
  <tbody>
    <tr>
      <th>Key</th>
      <td>self_reported_ethnicity_ontology_term_id</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
      <td>
        categorical with <code>str</code> categories.<br><br>If <code>tissue_type</code> is <code>"cell line"</code>, this MUST be <code>"na"</code><br><br>If <code>organism_ontolology_term_id</code> is NOT
        <code>"NCBITaxon:9606"</code> for <i>Homo sapiens</i>, this MUST be <code>"na"</code>.<br><br>Otherwise, if
        <code>organism_ontolology_term_id</code> is
        <code>"NCBITaxon:9606"</code> for <i>Homo sapiens</i>, this MUST be <code>"unknown"</code> if unavailable; otherwise, this MUST meet the following requirements:<br /><br />
        <ul>
          <li>
            The value MUST be formatted as one or more AfPO or HANCESTRO terms in ascending lexical order separated by the delimiter <code>" || "</code> with no duplication of terms.
          </li>
          <li>
            Each AfPO or HANCESTRO term MUST be a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0601"><code>"HANCESTRO:0601"</code></a> for <i>ethnicity category</i> or <a href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0602"><code>"HANCESTRO:0602"</code></a> for <i>geography-based population category</i>.<br><br>
          </li>
            For example, if the terms are <code>"HANCESTRO:0590</code> and <code>HANCESTRO:0580"</code> then the value of <code>self_reported_ethnicity_ontology_term_id</code> MUST be <code>"HANCESTRO:0580 || HANCESTRO:0590"</code>.<br><br>
        </ul>
      </td>
    </tr>
  </tbody>
</table>
<br />   

### sex_ontology_term_id

<table><tbody>
    <tr>
      <th>Key</th>
      <td>sex_ontology_term_id</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories.<br><br>If <code>tissue_type</code> is <code>"cell line"</code>, this MUST be <code>"na"</code>.<br><br>If unavailable, this MUST be <code>"unknown"</code>.<br><br>If <code>organism_ontology_term_id</code> is <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A6239"><code>"NCBITaxon:6239"</code></a> for <i>Caenorhabditis elegans</i>, this MUST be <a href="https://www.ebi.ac.uk/ols4/ontologies/pato/classes?obo_id=PATO%3A0000384"><code>"PATO:0000384"</code></a> for <i>male</i> or <a href="https://www.ebi.ac.uk/ols4/ontologies/pato/classes?obo_id=PATO%3A0001340"><code>"PATO:0001340"</code></a> for <i>hermaphrodite</i>; otherwise, this MUST be one of:<br><br>
        <ul>
        <li><a href="https://www.ebi.ac.uk/ols4/ontologies/pato/classes?obo_id=PATO%3A0000383"><code>"PATO:0000383"</code></a> for  <i>female</i></li>
        <li><a href="https://www.ebi.ac.uk/ols4/ontologies/pato/classes?obo_id=PATO%3A0000384"><code>"PATO:0000384"</code></a> for  <i>male</i></li>
        <li><a href="https://www.ebi.ac.uk/ols4/ontologies/pato/classes?obo_id=PATO%3A0001340"><code>"PATO:0001340"</code></a> for  <i>hermaphrodite</i></li>
        </ul>
        </td>
    </tr>
</tbody></table>
<br>


### suspension_type

<table><tbody>
    <tr>
      <th>Key</th>
      <td>suspension_type</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be <code>"cell"</code>, <code>"nucleus"</code>, or <code>"na"</code>.<br>
        <br>This MUST be the correct type for the corresponding assay:
          <br><br>
          <table>
          <thead>
          <tr>
          <th>For Assay</th>
          <th>MUST Use</th>
          </tr>
          </thead>
          <tbody>
            <tr>
              <td><i>10x transcription profiling</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030080"><code>EFO:0030080</code></a>] and its descendants</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr> 
            <tr>
              <td><i>ATAC-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0007045"><code>EFO:0007045</code></a>] and its descendants</td>
              <td><code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>BD Rhapsody Targeted mRNA</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0700004"><code>EFO:0700004</code></a>]</td>
              <td><code>"cell"</code></td>
           </tr>
            <tr>
              <td><i>BD Rhapsody Whole Transcriptome Analysis</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0700003"><code>EFO:0700003</code></a>]</td>
              <td><code>"cell"</code></td>
           </tr>
            <tr>
              <td><i>CEL-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008679"><code>EFO:0008679</code></a>]</td>
              <td><code>"cell"</code></td>
           </tr>
            <tr>
              <td><i>CEL-seq2</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010010"><code>EFO:0010010</code></a>] and its descendants</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>DroNc-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008720"><code>EFO:0008720</code></a>]</td>
              <td><code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>Drop-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008722"><code>EFO:0008722</code></a>]</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>GEXSCOPE technology</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0700011"><code>EFO:0700011</code></a>]</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr> 
            <tr>
              <td><i>inDrop</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008780"><code>EFO:0008780</code></a>]</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>MARS-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008796"><code>EFO:0008796</code></a>]</td>
              <td><code>"cell"</code></td>
           </tr>
           <tr>
             <td><i>mCT-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030060"><code>EFO:0030060</code></a>]</td>
             <td><code>"cell"</code> or <code>"nucleus"</code></td>
          </tr>
          <tr>
            <td><i>MERFISH</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008992"><code>EFO:0008992</code></a>]</td>
            <td><code>"na"</code></code></td>
          </tr>
          <tr>
           <td><i>methylation profiling by high throughput sequencing</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0002761"><code>EFO:0002761</code></a>] and its descendants</td>
          <td><code>"nucleus"</code></td>
         </tr>
          <tr>
              <td><i>microwell-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030002"><code>EFO:0030002</code></a>]</td>
              <td><code>"cell"</code></td>
           </tr>    
            <tr>
              <td><i>Patch-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008853"><code>EFO:0008853</code></a>]</td>
              <td><code>"cell"</code></td>
           </tr>
            <tr>
              <td><i>Quartz-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008877"><code>EFO:0008877</code></a>]</td>
              <td><code>"cell"</code></td>
           </tr>
          <tr>
            <td><i>ScaleBio single cell RNA sequencing</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0022490"><code>EFO:0022490</code></a>]</td>
           <td><code>"cell"</code> or <code>"nucleus"</code></td>
          </tr>
            <tr>
              <td><i>sci-Plex</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030026"><code>EFO:0030026</code></a>]</td>
              <td><code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>sci-RNA-seq3</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030028"><code>EFO:0030028</code></a>]</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>Seq-Well</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008919"><code>EFO:0008919</code></a>] and its descendants</td>
              <td><code>"cell"</code></td>
           </tr>
            <tr>
              <td><i>Smart-like</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010184"><code>EFO:0010184</code></a>] and its descendants</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>spatial transcriptomics</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008994"><code>EFO:0008994</code></a>] and its descendants</td>
              <td><code>"na"</code></td>
           </tr> 
            <tr>
              <td><i>SPLiT-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0009919"><code>EFO:0009919</code></a>] and its descendants</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr> 
            <tr>
              <td><i>STRT-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008953"><code>EFO:0008953</code></a>] and its descendants</td>
              <td><code>"cell"</code></td>
           </tr>
            <tr>
              <td><i>TruDrop</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0700010"><code>EFO:0700010</code></a>]</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr> 
          </tbody></table>
          <br>If the assay does not appear in this table, the most appropriate value MUST be selected and <a href="mailto:cellxgene@chanzuckerberg.com">the curation team informed</a> during submission so that the assay can be added to the table.<br>
        </td>
    </tr>
</tbody></table>
<br>

### tissue_ontology_term_id

<table><tbody>
    <tr>
      <th>Key</th>
      <td>tissue_ontology_term_id</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
      <td>
        categorical with <code>str</code> categories.<br><br>If <code>tissue_type</code> is <code>"cell line"</code>, this MUST be a Cellosaurus term.<br><br>If <code>tissue_type</code> is <code>"primary cell culture"</code>, this MUST follow the requirements for <code>cell_type_ontology_term_id</code>.<br><br>If <code>tissue_type</code> is <code>"organoid"</code>, this MUST NOT be <a href="https://www.ebi.ac.uk/ols4/ontologies/uberon/classes?obo_id=UBERON%3A0000922"><code>UBERON:0000922</code></a> for <i>embryo</i>. If the organoid is an embryoid, it is STRONGLY RECOMMENDED that the value is <a href="https://www.ebi.ac.uk/ols4/ontologies/uberon/classes?obo_id=UBERON%3A0014374"><code>UBERON:0014374</code></a> for <i>embryoid body</i>. If the organoid is a gastruloid, it is STRONGLY RECOMMENDED that the value is <a href="https://www.ebi.ac.uk/ols4/ontologies/uberon/classes?obo_id=UBERON%3A0004734"><code>UBERON:0004734</code></a> for <i>gastrula</i>.<br><br>Otherwise, if <code>tissue_type</code> is <code>"organoid"</code> or <code>"tissue"</code> then:<br><br>
        <table>
          <thead>
            <tr>
              <th>For <code>organism_ontology_term_id</code></th>
              <th>Value</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td>
                <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A6239"><code>"NCBITaxon:6239"</code></a><br>for <i>Caenorhabditis elegans</i>
              </td>
              <td>
                MUST be either the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/uberon/classes?obo_id=UBERON%3A0001062"><code>UBERON:0001062</code></a> for <i>anatomical entity</i> or the most accurate descendant<br>of <a href="https://www.ebi.ac.uk/ols4/ontologies/wbbt/classes?obo_id=WBBT%3A0005766"><code>WBbt:0005766</code></a> for <i>Anatomy</i> excluding <a href="https://www.ebi.ac.uk/ols4/ontologies/wbbt/classes?obo_id=WBBT%3A0007849"><code>WBbt:0007849</code></a> for <i>hermaphrodite</i>,<br><a href="https://www.ebi.ac.uk/ols4/ontologies/wbbt/classes?obo_id=WBBT%3A0007850"><code>WBbt:0007850</code></a> for <i>male</i>, <a href="https://www.ebi.ac.uk/ols4/ontologies/wbbt/classes?obo_id=WBBT%3A0008595"><code>WBbt:0008595</code></a> for <i>female</i>, <a href="https://www.ebi.ac.uk/ols4/ontologies/wbbt/classes?obo_id=WBBT%3A0004017"><code>WBbt:0004017</code></a> for <i>Cell</i><br>and its descendants, and <a href="https://www.ebi.ac.uk/ols4/ontologies/wbbt/classes?obo_id=WBBT%3A0006803"><code>WBbt:00006803</code></a> for <i>Nucleus</i> and its descendants
              </td>
            </tr>
            <tr>
              <td>
                <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A7955"><code>"NCBITaxon:7955"</code></a><br>for <i>Danio rerio</i>
              </td>
              <td>
                MUST be either the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/uberon/classes?obo_id=UBERON%3A0001062"><code>UBERON:0001062</code></a> for <i>anatomical entity</i> or the most accurate descendant of<br><a href="https://www.ebi.ac.uk/ols4/ontologies/zfa/classes?obo_id=ZFA%3A0100000"><code>ZFA:0100000</code></a> for <i>zebrafish anatomical entity</i> excluding <a href="https://www.ebi.ac.uk/ols4/ontologies/zfa/classes?obo_id=ZFA%3A0001093"><code>ZFA:0001093</code></a> for<br><i>unspecified</i> and <a href="https://www.ebi.ac.uk/ols4/ontologies/zfa/classes?obo_id=ZFA%3A0009000"><code>ZFA:0009000</code></a> for <i>cell</i> and its descendants
              </td>
            </tr>
            <tr>
              <td>
                <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A7227"><code>"NCBITaxon:7227"</code></a><br>for <i>Drosophila melanogaster</i>
              </td>
              <td>
                MUST be either the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/uberon/classes?obo_id=UBERON%3A0001062"><code>UBERON:0001062</code></a> for <i>anatomical entity</i> or the most accurate descendant of<br><a href="https://www.ebi.ac.uk/ols4/ontologies/fbbt/classes?obo_id=FBBT%3A10000000"><code>FBbt:10000000</code></a> for <i>anatomical entity</i> excluding <a href="https://www.ebi.ac.uk/ols4/ontologies/fbbt/classes?obo_id=FBbt%3A00007002"><code>FBbt:00007002</code></a> for <i>cell</i><br>and its descendants
              </td>
              </tr>    
              <tr>
              <td>
                For all other organisms
              </td>
              <td>
              MUST be the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/uberon/classes?obo_id=UBERON%3A0001062"><code>UBERON:0001062</code></a> for <i>anatomical entity</i>
              </td>
              </tr>
          </tbody>
        </table>
      </td>
  </tr>
</tbody></table>
<br>

### tissue_type

<table><tbody>
    <tr>
      <th>Key</th>
      <td>tissue_type</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be one of:
          <ul>
            <li><code>"cell line"</code></li>
            <li><code>"organoid"</code></li>
            <li><code>"primary cell culture"</code></li>
            <li><code>"tissue"</code></li>
         </ul>
    </tr>
</tbody></table>
<br>

When a dataset is uploaded, CELLxGENE Discover MUST automatically add the matching human-readable name for the corresponding ontology term to the `obs` dataframe. Curators MUST NOT annotate the following columns.

### assay

<table><tbody>
    <tr>
      <th>Key</th>
      <td>assay</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be the human-readable name assigned to the value of <code>assay_ontology_term_id</code>. 
        </td>
    </tr>
</tbody></table>
<br>

### cell_type

<table><tbody>
    <tr>
      <th>Key</th>
      <td>cell_type</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories.<br><br>This MUST be <code>"na"</code> if the value of  <code>cell_type_ontology_term_id</code> is <code>"na"</code>.<br><br>This MUST be <code>"unknown"</code> if the value of  <code>cell_type_ontology_term_id</code> is <code>"unknown"</code>.<br><br>Otherwise, this MUST be the human-readable name assigned to the value of <code>cell_type_ontology_term_id</code>.
        </td>
    </tr>
</tbody></table>
<br>


### development_stage

<table><tbody>
    <tr>
      <th>Key</th>
      <td>development_stage</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories.<br><br>This MUST be <code>"na"</code> if the value of <code>development_stage_ontology_term_id</code> is <code>"na"</code>.<br><br>This MUST be <code>"unknown"</code> if the value of <code>development_stage_ontology_term_id</code> is <code>"unknown"</code>.<br><br>Otherwise, this MUST be the human-readable name assigned to the value of <code>development_stage_ontology_term_id</code>.
        </td>
    </tr>
</tbody></table>
<br>

When a dataset is uploaded, CELLxGENE Discover MUST annotate a unique observation identifier for each cell. Curators MUST NOT annotate the following column.

### observation_joinid

<table><tbody>
    <tr>
      <th>Key</th>
      <td>observation_joinid</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>str</code>
    </tr>
</tbody></table>
<br>

### self_reported_ethnicity

<table><tbody>
    <tr>
      <th>Key</th>
      <td>self_reported_ethnicity</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be <code>"na"</code> if the value of <code>self_reported_ethnicity_ontology_term_id</code> is <code>"na"</code>. This MUST be <code>"unknown"</code> if the value of <code>self_reported_ethnicity_ontology_term_id</code> is <code>"unknown"</code>. Otherwise, this MUST be one or more human-readable names for the terms in <code>self_reported_ethnicity_ontology_term_id</code> in the same order separated by the delimiter <code>" || "</code>.<br><br> For example, if the value of <code>self_reported_ethnicity_ontology_term_id</code> is <code>"HANCESTRO:0005 || HANCESTRO:0014"</code> then the value of <code>self_reported_ethnicity</code> MUST be <code>"European || Hispanic or Latin American"</code>.
        </td>
    </tr>
</tbody></table>
<br>

### sex

<table><tbody>
    <tr>
      <th>Key</th>
      <td>sex</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories.<br><br>This MUST be <code>"na"</code> if the value of  <code>sex_ontology_term_id</code> is <code>"na"</code>.<br><br>This MUST be <code>"unknown"</code> if the value of  <code>sex_ontology_term_id</code> is <code>"unknown"</code>.<br><br>Otherwise, this MUST be the human-readable name assigned to the value of <code>sex_ontology_term_id</code>.
        </td>
    </tr>
</tbody></table>
<br>

## `obsm` (Embeddings)


The value for each `str` key MUST be a  `numpy.ndarray` of shape `(n_obs, m)`, where `n_obs` is the number of rows in `X` and `m >= 1`. 

To display a dataset in CELLxGENE Explorer, Curators MUST annotate **one or more** embeddings of at least two-dimensions (e.g. tSNE, UMAP, PCA, spatial coordinates) as `numpy.ndarrays` in `obsm`.<br><br>

### spatial

<table><tbody>
    <tr>
      <th>Key</th>
      <td>spatial</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate if <code>uns['spatial']['is_single']</code> is <code>True</code>.<br><br>Curator MAY annotate if <code>uns['spatial']['is_single']</code> is <code>False</code>.
      <br><br>Otherwise, this key MUST NOT be present.</td>
    </tr>
        <tr>
      <th>Value</th>
        <td><code>numpy.ndarray</code> with the following requirements<br><br>
          <ul>
          <li>MUST have the same number of rows as <code>X</code> and MUST include at least two columns</li>
          <li>MUST be a <a href="https://numpy.org/doc/stable/reference/generated/numpy.dtype.kind.html"><code>numpy.dtype.kind</code></a> of <code>"f"</code>, <code>"i"</code>, or "<code>u"</code></li>
          <li>MUST NOT contain any <a href="https://numpy.org/devdocs/reference/constants.html#numpy.inf">positive infinity (<code>numpy.inf</code>)</a> or <a href="https://numpy.org/devdocs/reference/constants.html#numpy.NINF">negative infinity (<code>numpy.NINF</code>)</a> values </li>
          <li>MUST NOT contain <a href="https://numpy.org/devdocs/reference/constants.html#numpy.nan">Not a Number (<code>numpy.nan</code>)</a> values</li></ul><br>If <code>assay_ontology_term_id</code>is a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> and <code>uns['spatial']['is_single']</code> is <code>True</code>, the array MUST be created from the corresponding <code>pxl_row_in_fullres</code> and <code>pxl_col_in_fullres</code> fields from <code>tissue_positions_list.csv</code> or <code>tissue_positions.csv</code>. See <a href="https://www.10xgenomics.com/support/software/space-ranger/analysis/outputs/spatial-outputs">Space Ranger Spatial Outputs</a>.
        </td>
    </tr>
</tbody></table>
<br>

### X_{suffix}

<table><tbody>
    <tr>
      <th>Key</th>
      <td>X_{suffix} with the following requirements:<br><br>
      <ul>
        <li>{suffix} MUST be at least one character in length.</li>
        <li>The first character of {suffix} MUST be a letter of the alphabet and the remaining characters MUST be alphanumeric characters, <code>'_'</code>, <code>'-'</code>, or <code>'.'</code> (This is equivalent to the regular expression pattern <code>"^[a-zA-Z][a-zA-Z0-9_.-]*$"</code>.)</li>
         <li>{suffix} MUST NOT be <code>"spatial"</code>.
      </ul><br>
      {suffix} is presented as text to users in the <b>Embedding Choice</b> selector in CELLxGENE Explorer so it is STRONGLY RECOMMENDED that it be descriptive.<br><br>See also <code>default_embedding</code> in <code>uns</code>.</td>
    </tr>
    <tr>
      <th>Annotator</th>
         <td>Curator MUST annotate if <code>assay_ontology_term_id</code> is neither a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> nor <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030062"><code>"EFO:0030062"</code></a> for <i>Slide-seqV2</i>.<br><br>Curator MAY annotate if <code>assay_ontology_term_id</code> is either a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> or <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030062"><code>"EFO:0030062"</code></a> for <i>Slide-seqV2</i>.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>numpy.ndarray</code> with the following requirements<br><br>
          <ul>
          <li>MUST have the same number of rows as <code>X</code> and MUST include at least two columns</li>
          <li>MUST be a <a href="https://numpy.org/doc/stable/reference/generated/numpy.dtype.kind.html"><code>numpy.dtype.kind</code></a> of <code>"f"</code>, <code>"i"</code>, or "<code>u"</code></li>
          <li>MUST NOT contain any <a href="https://numpy.org/devdocs/reference/constants.html#numpy.inf">positive infinity (<code>numpy.inf</code>)</a> or <a href="https://numpy.org/devdocs/reference/constants.html#numpy.NINF">negative infinity (<code>numpy.NINF</code>)</a> values </li>
          <li>MUST NOT contain all <a href="https://numpy.org/devdocs/reference/constants.html#numpy.nan">Not a Number (<code>numpy.nan</code>) </a>
values</li></ul>
        </td>
    </tr>
</tbody></table>
<br>

## `obsp`

The size of the ndarray stored for a key in `obsp` MUST NOT be zero.
<br>

## `uns` (Dataset Metadata)

`uns` is a ordered dictionary with a `str` key. The data stored as a value for a key in `uns` MUST be `True`, `False`, `None`, or its size MUST NOT be zero.

Curators MUST annotate the following keys and values in `uns`:

### organism_ontology_term_id

<table><tbody>
    <tr>
      <th>Key</th>
      <td>organism_ontology_term_id</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>str</code>. One of the following terms MUST be used:
          <br><br>
          <table>
          <thead>
          <tr>
          <th>For Organism</th>
          <th>MUST Use</th>
          </tr>
          </thead>
          <tbody>
         <tr>
            <td><i>Caenorhabditis elegans</i></td>
            <td>
              <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A6239"><code>"NCBITaxon:6293"</code></a>
            </td>
          </tr>
         <tr>
            <td><i>Callithrix jacchus</i></td>
            <td>
              <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9483"><code>"NCBITaxon:9483"</code></a>
            </td>
          </tr>
          <tr>
            <td><i>Danio rerio</i></td>
            <td>
              <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A7955"><code>"NCBITaxon:7955"</code></a></td>
          </tr>
          <tr>
            <td><i>Drosophila melanogaster</i></td>
            <td>
              <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A7227"><code>"NCBITaxon:7227"</code></a>
            </td>
          </tr>
          <tr>
            <td><i>Gorilla gorilla gorilla</i></td>
            <td>
              <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9595"><code>"NCBITaxon:9595"</code></a>
            </td>
           </tr>
            <tr>
              <td><i>Homo sapiens</i></td>
              <td>
                <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9606"><code>"NCBITaxon:9606"</code></a>
              </td>
            </tr>
            <tr>
              <td><i>Macaca fascicularis<br>and its descendants</i></td>
              <td>
                <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9541"><code>"NCBITaxon:9541"</code></a>or one of its descendants
              </td>
            </tr>
          <tr>
              <td><i>Macaca mulatta<br>and its descendants</i></td>
              <td>
                <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9544"><code>"NCBITaxon:9544"</code></a>or one of its descendants
              </td>
            </tr>
            <tr>
              <td><i>Microcebus murinus</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A30608"><code>"NCBITaxon:30608"</code></a></td>
            </tr>
            <tr>
              <td><i>Mus musculus<br>and its descendants</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A10090"><code>"NCBITaxon:10090"</code></a>or one of its descendants</td>
            </tr>
            <tr>
              <td><i>Oryctolagus cuniculus<br>and its descendants</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9986"><code>"NCBITaxon:9986"</code></a>or one of its descendants</td>
            </tr>
            <tr>
              <td><i>Pan troglodytes<br>and its descendants</i></td>
              <td>
                <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9598"><code>"NCBITaxon:9598"</code></a>or one of its descendants
             </td>
            </tr>
            <tr>
              <td><i>Rattus norvegicus<br>and its descendants</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A10116"><code>"NCBITaxon:10116"</code></a>or one of its descendants</td>
            </tr>
            <tr>
              <td><i>Sus scrofa<br>and its descendants</i></td>
              <td>
               <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9823"><code>"NCBITaxon:9823"</code></a>or one of its descendants
            </td>
          </tr>
          </tbody></table>
        </td>
    </tr>
</tbody></table>
<br>


### spatial

<table><tbody>
    <tr>
      <th>Key</th>
      <td>spatial</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate if <code>assay_ontology_term_id</code> is descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> or <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030062"><code>"EFO:0030062"</code></a> for <i>Slide-seqV2</i>; otherwise, this key MUST NOT be present.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>dict</code>. The requirements for the key-value pairs are documented in the following sections:
          <ul>
          <li>spatial['is_single']</li>         
          <li>spatial[<i>library_id</i>]</li>
          <li>spatial[<i>library_id</i>]['images']</li>
          <li>spatial[<i>library_id</i>]['images']['fullres']</li>
          <li>spatial[<i>library_id</i>]['images']['hires']</li>
          <li>spatial[<i>library_id</i>]['scalefactors']</li>
          <li>spatial[<i>library_id</i>]['scalefactors']['spot_diameter_fullres']</li>
          <li>spatial[<i>library_id</i>]['scalefactors']['tissue_hires_scalef']</li>
         </ul><br>Additional key-value pairs MUST NOT be present.
        </td>
    </tr>
</tbody></table>
<br>

#### is_single

<table><tbody>
    <tr>
      <th>Key</th>
      <td>is_single</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate if <code>assay_ontology_term_id</code> is a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> or <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030062"><code>"EFO:0030062"</code></a> for <i>Slide-seqV2</i>; otherwise, this key MUST NOT be present.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>bool</code>. This MUST be <code>True</code>:
        <ul>
        <li>if <code>assay_ontology_term_id</code> is a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> and the dataset represents one Space Ranger output for a single tissue section
      </li>
      <li> if <code>assay_ontology_term_id</code> is <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030062"><code>"EFO:0030062"</code></a> for <i>Slide-seqV2</i> and the dataset represents the output for a single array on a puck </li>
        </ul>
        Otherwise, this MUST be <code>False</code>.
        </td>
    </tr>
</tbody></table>
<br>

#### spatial[_library_id_]
<table><tbody>
    <tr>
      <th>Key</th>
      <td>Identifier for the Visium library</td>
    </tr>
    <tr>
      <th>Annotation</th>
      <td>Curator MUST annotate if <code>assay_ontology_term_id</code> is a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> and <code>uns['spatial']['is_single']</code> is <code>True</code>; otherwise, this key MUST NOT be present.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>dict</code>. There MUST be only one <code><i>library_id</i></code>.</td>
    </tr>
</tbody></table>
<br>

#### spatial[_library_id_]['images']
<table><tbody>
    <tr>
      <th>Key</th>
      <td>images</td>
    </tr>
    <tr>
      <th>Annotation</th>
      <td>Curator MUST annotate if <code>assay_ontology_term_id</code> is a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> and <code>uns['spatial']['is_single']</code> is <code>True</code>; otherwise, this key MUST NOT be present.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          <code>dict</code>
        </td>
    </tr>
</tbody></table>
<br>

#### spatial[_library_id_]['images']['fullres']
<table><tbody>
    <tr>
      <th>Key</th>
      <td>fullres</td>
    </tr>
    <tr>
      <th>Annotation</th>
      <td>Curator MAY annotate if <code>assay_ontology_term_id</code> is a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> and <code>uns['spatial']['is_single']</code> is <code>True</code>; otherwise, this key MUST NOT be present.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          The full resolution image MUST be converted to a<code>numpy.ndarray</code> with the following requirements:<br><br>
          <ul>
          <li>The length of <code>numpy.ndarray.shape</code> MUST be <code>3</code></li>
          <li>The <code>numpy.ndarray.dtype</code> MUST be <code>numpy.uint8</code></li>
          <li>The <code>numpy.ndarray.shape[2]</code> MUST be either <code>3</code> (RGB color model for example) or <code>4</code> (RGBA color model for example)</li>
          </ul>
        </td>
    </tr>
</tbody></table>
<br>

#### spatial[_library_id_]['images']['hires']
<table><tbody>
    <tr>
      <th>Key</th>
      <td>hires</td>
    </tr>
    <tr>
      <th>Annotation</th>
      <td>Curator MUST annotate if <code>assay_ontology_term_id</code> is a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> and <code>uns['spatial']['is_single']</code> is <code>True</code>; otherwise, this key MUST NOT be present.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          <code>tissue_hires_image.png</code> MUST be converted to a<code>numpy.ndarray</code> with the following requirements:<br><br>
          <ul>
          <li>The length of <code>numpy.ndarray.shape</code> MUST be <code>3</code></li>
          <li>The <code>numpy.ndarray.dtype</code> MUST be <code>numpy.uint8</code></li>
          <li>If <code>assay_ontology_term_id</code> is <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0022860"><code>"EFO:0022860"</code></a> for <i>Visium CytAssist Spatial Gene Expression, 11mm</i>, the largest dimension in <code>numpy.ndarray.shape[:2]</code> MUST be <code>4000</code>pixels; otherwise, the largest dimension in <code>numpy.ndarray.shape[:2]</code> MUST be <code>2000</code>pixels. See <a href="https://www.10xgenomics.com/support/software/space-ranger/analysis/outputs/spatial-outputs">Space Ranger Spatial Outputs</a></li>
          <li>The <code>numpy.ndarray.shape[2]</code> MUST be either <code>3</code> (RGB color model for example) for <code>4</code> (RGBA color model for example)</li>
          </ul>
        </td>
    </tr>
</tbody></table>
<br>

#### spatial[_library_id_]['scalefactors']
<table><tbody>
    <tr>
      <th>Key</th>
      <td>scalefactors</td>
    </tr>
    <tr>
      <th>Annotation</th>
      <td>Curator MUST annotate if <code>assay_ontology_term_id</code> is a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> and <code>uns['spatial']['is_single']</code> is <code>True</code>; otherwise, this key MUST NOT be present.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          <code>dict</code>
        </td>
    </tr>
</tbody></table>
<br>

#### spatial[_library_id_]['scalefactors']['spot_diameter_fullres']
<table><tbody>
    <tr>
      <th>Key</th>
      <td>spot_diameter_fullres</td>
    </tr>
    <tr>
      <th>Annotation</th>
      <td>Curator MUST annotate if <code>assay_ontology_term_id</code> is a descendant <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> and <code>uns['spatial']['is_single']</code> is <code>True</code>; otherwise, this key MUST NOT be present.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>float</code>. This must be the value of the <code>spot_diameter_fullres</code> field from <code>scalefactors_json.json</code>. See <a href="https://www.10xgenomics.com/support/software/space-ranger/analysis/outputs/spatial-outputs">Space Ranger Spatial Outputs</a>.</td>
    </tr>
</tbody></table>
<br>

#### spatial[_library_id_]['scalefactors']['tissue_hires_scalef']
<table><tbody>
    <tr>
      <th>Key</th>
      <td>tissue_hires_scalef</td>
    </tr>
    <tr>
      <th>Annotation</th>
      <td>Curator MUST annotate if <code>assay_ontology_term_id</code> is a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a> for <i>Visium Spatial Gene Expression</i> and <code>uns['spatial']['is_single']</code> is <code>True</code>; otherwise, this key MUST NOT be present.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          <code>float</code>. This must be the value of the <code>tissue_hires_scalef</code> field from <code>scalefactors_json.json</code>. See <a href="https://www.10xgenomics.com/support/software/space-ranger/analysis/outputs/spatial-outputs">Space Ranger Spatial Outputs</a>.
        </td>
    </tr>
</tbody></table>
<br>

### title

<table><tbody>
    <tr>
      <th>Key</th>
      <td>title</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          <code>str</code>. This text describes and differentiates the dataset from other datasets in the same collection. It is displayed on a page in CELLxGENE Discover that also has the collection name. To illustrate, the first dataset name in the <a href="https://cellxgene.cziscience.com/collections/b52eb423-5d0d-4645-b217-e1c6d38b2e72">Cells of the adult human heart collection</a> is "All — Cells of the adult human heart".<br><br>It is STRONGLY RECOMMENDED that each dataset <code>title</code> in a collection is unique and does not depend on other metadata such as a different  <code>assay</code> to disambiguate it from other datasets in the collection.
        </td>
    </tr>
</tbody></table>
<br>

​Curators MAY also annotate the following optional keys and values in `uns`. If the key is present, then its value MUST NOT be empty.
​
### batch_condition

<table><tbody>
    <tr>
      <th>Key</th>
      <td>batch_condition</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MAY annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          <code>list[str]</code>. <code>str</code> values MUST refer to cell metadata keys in <code>obs</code>. Together, these keys define the <i>batches</i> that a normalization or integration algorithm should be aware of. For example if <code>"patient"</code> and <code>"seqBatch"</code> are keys of vectors of cell metadata, either <code>["patient"]</code>, <code>["seqBatch"]</code>, or <code>["patient", "seqBatch"]</code> are valid values.
        </td>
    </tr>
</tbody></table>
<br>


### {column}_colors

<table><tbody>
  <tr>
    <th>Key</th>
      <td>
        {column}_colors where {column} MUST be the name of a <code>category</code> data type column in <code>obs</code> that<br> is annotated by the data submitter or curator. The following columns that are annotated by CELLxGENE<br> Discover MUST NOT be specified as {column}:<br><br>
      <ul>
        <li>assay</li>
        <li>cell_type</li>
        <li>development_stage</li>
        <li>disease</li>
        <li>self_reported_ethnicity</li>
        <li>sex</li>
        <li>tissue</li>       
      </ul><br>
      Instead annotate {column}_ontology_term_id_colors for these columns such as <code>assay_ontology_term_id</code>.<br><br>
    </td>
  </tr>
  <tr>
    <th>Annotator</th>
    <td>Curator MAY annotate.</td>
  </tr>
  <tr>
    <th>Value</th>
      <td>
        <code>numpy.ndarray</code>. This MUST be a 1-D array of shape <code>(, c)</code>, where <code>c</code> is greater than or equal to the<br> number of categories in the {column} as calculated by:<br><br>
           <samp>anndata.obs.{column}.cat.categories.size</samp><br><br>
        The color code at the Nth position in the <code>ndarray</code> corresponds to the Nth category of <samp>anndata.obs.{column}.cat.categories</samp>.<br><br>For example, if <code>cell_type_ontology_term_id</code> includes two categories:<br><br>
        <samp>anndata.obs.cell_type_ontology_term_id.cat.categories.values</samp><br><br>
        <samp>array(['CL:0000057', 'CL:0000115'],
      dtype='object')</samp><br><br>then <code>cell-type_ontology_term_id_colors</code> MUST contain two or more colors such as:<br><br>
        <samp>['aqua' 'blueviolet']</samp><br><br>where <code>'aqua'</code> is the color assigned to <code>'CL:0000057'</code> and <code>'blueviolet'</code> is the color assigned to<br> <code>'CL:0000115'</code>.<br><br>All elements in the <code>ndarray</code> MUST use the same color model, limited to:<br><br>
          <table>
          <thead>
            <tr>
              <th>Color Model</th>
              <th>Element Format</th>
            </tr>
          </thead><tbody>
            <tr>
              <td>
              <a
              href="https://www.w3.org/TR/css-color-4/#named-colors"
              ><i>Named Colors </i>
              </a>
            </td>
              <td><code>str</code>. MUST be a case-insensitive CSS4 color name with no spaces such as<br> <code>"aliceblue"</code>
            </td>
            </tr>
            <tr>
             <td>
              <a
              href="https://www.w3.org/TR/css-color-4/#hex-notation"
              ><i>Hex Triplet</i>
              </a>
            </td>
              <td><code>str</code>. MUST start with <code>"#"</code> immediately followed by six case-insensitive hexadecimal<br> characters as in <code>"#08c0ff"</code></td>
            </tr>
          </tbody></table>
        </td>
    </tr>
</tbody></table>
<br>

### default_embedding

<table><tbody>
    <tr>
      <th>Key</th>
      <td>default_embedding</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MAY annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          <code>str</code>. The value MUST match a key to an embedding in <code>obsm</code> for the embedding to display by default in CELLxGENE Explorer.
        </td>
    </tr>
</tbody></table>
<br>

### X_approximate_distribution

<table><tbody>
    <tr>
      <th>Key</th>
      <td>X_approximate_distribution</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MAY annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          <code>str</code>. CELLxGENE Discover runs a heuristic to detect the approximate distribution of the data in X so that it can accurately calculate statistical properties of the data. This field enables the curator to override this heuristic and specify the data distribution explicitly. The value MUST be <code>"count"</code> (for data whose distributions are best approximated by counting distributions like Poisson, Binomial, or Negative Binomial) or <code>"normal"</code> (for data whose distributions are best approximated by the Gaussian distribution.)
        </td>
    </tr>
</tbody></table>
<br>

Curators MUST NOT annotate the following keys and values in `uns`.

When a dataset is uploaded, CELLxGENE Discover MUST automatically add the `citation` key and set its value.

### citation

<table><tbody>
    <tr>
      <th>Key</th>
      <td>citation</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>str</code>. Its format MUST use the following template:
          <br><br>
          <table>
          <thead>
          <tr>
          <th>Citation Element</th>
          <th>Value</th>
          </tr>
          </thead>
          <tbody>
            <tr>
              <td><i><code>"Publication: "</code></i></td>
              <td>Publication DOI url for the collection<br><br> This element MUST only be present if a<br> Publication DOI is defined for the collection;<br> otherwise, it MUST NOT be present.</td>
            </tr>
            <tr>
              <td><i><code>"Dataset Version: "</code></i></td>
              <td>Permanent url to this version of the dataset</td>
            </tr>
              <td><i><code>" curated and distributed by<br> CZ CELLxGENE Discover in Collection: "</code> </i></td>
              <td>Permanent url to the collection</td>
            </tr>
          </tbody></table>
          A citation for a H5AD dataset with a Publication DOI:<br><br>"<code><b>Publication:</b> https://doi.org/10.1126/science.abl4896 <b>Dataset Version:</b> https://datasets.cellxgene.cziscience.com/dbd8b789-3efa-4a63-9243-90cff64f2045.h5ad <b>curated and distributed by CZ CELLxGENE Discover in Collection:</b> https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5"</code><br><br>
          A citation for a RDS dataset without a Publication DOI:<br><br><code>"<b>Dataset Version:</b> https://datasets.cellxgene.cziscience.com/08ea16dc-3f4e-4c84-8692-74d70be22d12.rds <b>curated and distributed by CZ CELLxGENE Discover in Collection:</b> https://cellxgene.cziscience.com/collections/10bf5c50-8d85-4c5f-94b4-22c1363d9f31"</code><br><br>
        </td>
    </tr>
</tbody></table>
<br>

When a dataset is uploaded, CELLxGENE Discover MUST automatically add the `organism` key and set its value to the matching human-readable name for the corresponding ontology term.

### organism

<table><tbody>
    <tr>
      <th>Key</th>
      <td>organism</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be the human-readable name assigned to the value of <code>organism_ontology_term_id</code>.
        </td>
    </tr>
</tbody></table>
<br>

When a dataset is uploaded, CELLxGENE Discover MUST automatically add the `schema_reference` key and set its value to the permanent URL of this document. 

### schema_reference

<table><tbody>
    <tr>
      <th>Key</th>
      <td>schema_reference</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          This MUST be <code>"https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/7.0.0/schema.md"</code>.
        </td>
    </tr>
</tbody></table>
<br>

---

When a dataset is uploaded, CELLxGENE Discover MUST automatically add the `schema_version` key and its value to `uns`.

### schema_version

<table><tbody>
    <tr>
      <th>Key</th>
      <td>schema_version</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          This MUST be <code>"7.0.0"</code>.
        </td>
    </tr>
</tbody></table>
<br>

## `var` and `raw.var` (Gene Metadata)

`var` and `raw.var` are both of type [`pandas.DataFrame`](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html).

Curators MUST annotate the following columns in the `var` dataframe and if present, the `raw.var` dataframe.

### index of pandas.DataFrame

<table><tbody>
    <tr>
      <th>Key</th>
      <td>index of <code>pandas.DataFrame</code></td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>str</code>. The index of the <code>pandas.DataFrame</code> MUST contain unique identifiers for features. If present, the index of <code>raw.var</code> MUST be identical to the index of <code>var</code>.<br><br>If the feature is a RNA Spike-In Control Mix then the value MUST be an ERCC Spike-In identifier (e.g. <code>"ERCC-0003"</code>).<br><br>If the feature is a gene then the value MUST be the <code>gene_id</code> attribute from the corresponding gene reference documented in <a href="#required-gene-annotations">Required Gene Annotations</a> for either the <code>organism_ontology_term_id</code> or <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A2697049"><code>NCBITaxon:2697049</code></a> for <i>SARS-CoV-2</i>.<br><br>Version numbers MUST be removed from the <code>gene_id</code> if it is prefixed with <code>"ENS"</code> for <i>Ensembl stable identifier</i>. See <a href="https://ensembl.org/Help/Faq?id=488">I have an Ensembl ID, what can I tell about it from the ID?</a> For example, if the <code>gene_id</code> is <code>“ENSG00000186092.7”</code>, then the value MUST be <code>“ENSG00000186092”</code>.
        <br><br>
        </td>
    </tr>
</tbody></table>
<br>

Curators MUST annotate the following column only in the `var` dataframe. This column MUST NOT be present in `raw.var`:

### feature_is_filtered

<table><tbody>
    <tr>
      <th>Key</th>
      <td>feature_is_filtered</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          <code>bool</code>. When a raw matrix is not present, the value for all features MUST be <code>False</code>.<br><br>
          When both a raw and normalized matrix are present, this MUST be <code>True</code> if the feature was filtered out in the normalized matrix (<code>X</code>) but is present in the raw matrix (<code>raw.X</code>). The value for all cells of the given feature in the normalized matrix MUST be <code>0</code>. If a feature contains all zeroes in the normalized matrix, then either the corresponding feature in the raw matrix MUST be all zeroes or the value MUST be <code>True</code>.
        <td>
    </tr>
</tbody></table>
<br>

Curators MUST NOT annotate the following columns in the `var` dataframe and if present, the `raw.var` dataframe.

When a dataset is uploaded, CELLxGENE Discover MUST automatically add the matching human-readable name for the corresponding feature biotype, identifier, and the NCBITaxon term for the reference organism to the `var` and `raw.var` dataframes. In addition, it MUST
add the feature length and type.

### feature_biotype

<table><tbody>
    <tr>
      <th>Key</th>
      <td>feature_biotype</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>This MUST be <code>"gene"</code> or <code>"spike-in"</code>.  
        </td>
    </tr>
</tbody></table>
<br>

### feature_length

<table><tbody>
    <tr>
      <th>Key</th>
      <td>feature_length</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
        <code>uint</code> number of base-pairs (bps). The value is the median of the lengths of isoforms, reusing the median calculation from <a href="https://doi.org/10.1093/bioinformatics/btac561">GTFtools: a software package for analyzing various features of gene models.</a>
      </td>
    </tr>
</tbody></table>
<br>

### feature_name

<table><tbody>
    <tr>
      <th>Key</th>
      <td>feature_name</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>str</code>. If the <code>feature_biotype</code> is <code>"spike-in"</code> then this MUST be the ERCC Spike-In identifier appended with <code>" (spike-in control)"</code>.<br><br>If the <code>feature_biotype</code> is <code>"gene"</code> and a <code>gene_name</code> attribute is assigned to the <code>var.index</code> feature identifier in its corresponding gene reference, this MUST be the value of the <code>gene_name</code>. If a <code>gene_name</code> attribute is not assigned, then this MUST default to the <code>var.index</code> feature identifier. 
        </td>
    </tr>
</tbody></table>
<br>

### feature_reference

<table><tbody>
    <tr>
      <th>Key</th>
      <td>feature_reference</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>str</code>. This MUST be the reference organism for a feature:
          <br><br>
          <table>
          <thead>
          <tr>
          <th>Reference Organism</th>
          <th>MUST Use</th>
          </tr>
          </thead>
          <tbody>
         <tr>
            <td><i>Caenorhabditis elegans</i></td>
            <td>
              <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A6239"><code>"NCBITaxon:6293"</code></a>
            </td>
          </tr>
         <tr>
            <td><i>Callithrix jacchus</i></td>
            <td>
              <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9483"><code>"NCBITaxon:9483"</code></a>
            </td>
          </tr>
          <tr>
            <td><i>Danio rerio</i></td>
            <td>
              <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A7955"><code>"NCBITaxon:7955"</code></a></td>
          </tr>
          <tr>
            <td><i>Drosophila melanogaster</i></td>
            <td>
              <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A7227"><code>"NCBITaxon:7227"</code></a>
            </td>
          </tr>
          <tr>
            <td><i>Gorilla gorilla gorilla</i></td>
            <td>
              <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9595"><code>"NCBITaxon:9595"</code></a>
            </td>
           </tr>
            <tr>
              <td><i>Homo sapiens</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9606"><code>"NCBITaxon:9606"</code></a></td>
            </tr>
            <tr>
              <td><i>Macaca fascicularis</i></td>
              <td>
                <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9541"><code>"NCBITaxon:9541"</code></a>
              </td>
            </tr>
           <tr>
            <td><i>Macaca mulatta</i></td>
            <td>
              <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9544"><code>"NCBITaxon:9544"</code></a>
            </td>
            </tr>
            <tr>
              <td><i>Microcebus murinus</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A30608"><code>"NCBITaxon:30608"</code></a></td>
            </tr>
            <tr>
              <td><i>Mus musculus</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A10090"><code>"NCBITaxon:10090"</code></a></td>
            </tr>
            <tr>
              <td><i>Oryctolagus cuniculus</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9986"><code>"NCBITaxon:9986"</code></a></td>
            </tr>
            <tr>
              <td><i>Pan troglodytes</i></td>
              <td>
                <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9598"><code>"NCBITaxon:9598"</code></a>
             </td>
            </tr>
            <tr>
              <td><i>Rattus norvegicus</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A10116"><code>"NCBITaxon:10116"</code></a></td>
            </tr>
            <tr>
              <td><i>SARS-CoV-2</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A2697049"><code>"NCBITaxon:2697049"</code></a></td>
            </tr>
            <tr>
              <td><i>Sus scrofa</i></td>
              <td>
               <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9823"><code>"NCBITaxon:9823"</code></a>
            </td>
          </tr>
            <tr>
              <td><i>ERCC Spike-Ins</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A32630"><code>"NCBITaxon:32630"</code></a></td>
            </tr>
          </tbody></table>
        </td>
    </tr>
</tbody></table>
<br>

### feature_type

<table><tbody>
    <tr>
      <th>Key</th>
      <td>feature_type</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>str</code>. If the <code>feature_biotype</code> is <code>"gene"</code> then this MUST be the gene type assigned to the feature identifier in <code>var.index</code>. If the <code>feature_biotype</code> is <code>"spike-in"</code> then this MUST be <code>"synthetic"</code>.<br><br>See  <a href="https://www.gencodegenes.org/pages/biotypes.html ">GENCODE</a> and <a href="https://useast.ensembl.org/info/genome/genebuild/biotypes.html ">Ensembl</a> references.
        </td>
    </tr>
</tbody></table>
<br>


## `varm`

The size of the ndarray stored for a key in `varm` MUST NOT be zero.
<br>

## `varp`
The size of the ndarray stored for a key in `varp` MUST NOT be zero.
<br>

---

## scATAC-seq assets

### Requirements

A Dataset MUST meet all of the following requirements to be eligible for scATAC-seq assets:
* <code>assay_ontology_term_id</code> values MUST be either all <i>paired assays</i> or <i>unpaired assays</i>
* <code>is_primary_data</code> values MUST be all <code>True</code>
* <code>organism_ontology_term_id</code> value MUST be either <code>"NCBITaxon:9606"</code> for <i>Homo sapiens</i> or <code>"NCBITaxon:10090"</code> for <i>Mus musculus</i> or one of its descendants. The value determines the required Chromosome Table.

If the <code>assay_ontology_term_id</code> values are all <i>paired assays</i> then the Dataset MAY have a fragments file asset.

If the <code>assay_ontology_term_id</code> values are all <i>unpaired assays</i> then the Dataset MUST have a fragments file asset.

## scATAC-seq Asset: Submitted Fragment File

This MUST be a gzipped tab-separated values (TSV) file.

The curator MUST annotate the following header-less columns. Additional columns and header lines beginning with `#` MUST NOT be included. Each row MUST represent a unique fragment.

### first column

<table><tbody>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          <code>str</code>. This MUST be the reference genome chromosome the fragment is located on.<br><br>If the value of organism_ontology_term_id</code> in the associated Dataset is <code>"NCBITaxon:9606"</code> for <i>Homo sapiens</i> then the first column value MUST be a value from the <code>Chromosome</code> column in the <a href="#human-grch38p14">Human Chromosome Table</a>.<br><br>
          If the value of <code>organism_ontology_term_id</code> in the associated Dataset is <code>"NCBITaxon:10090"</code> for <i>Mus musculus</i> then the first column value MUST be a value from the <code>Chromosome</code> column in the <a href="#mouse-grcm39">Mouse Chromosome Table</a>.
        </td>
    </tr>
</tbody></table>
<br>


### second column

<table><tbody>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>int</code>. This MUST be the 0-based start coordinate of the fragment.
        </td>
    </tr>
</tbody></table>
<br>

### third column

<table><tbody>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>int</code>. This MUST be the 0-based end coordinate of the fragment. The end position is exclusive, representing the position immediately following the fragment interval. The value MUST be greater than the start coordinate specified in the second column and less than or equal to the <code>Length</code> of the <code>Chromosome</code> specified in the first column, as specified in the appropriate Chromosome Table.
        </td>
    </tr>
</tbody></table>
<br>

### fourth column

<table><tbody>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>str</code>. This MUST be an observation identifier from the <a href="#index-of-pandasdataframe"><code>obs</code> index</a> of the associated Dataset. Every <code>obs</code> index value of the associated Dataset MUST appear at least once in this column.
        </td>
    </tr>
</tbody></table>
<br>

### fifth column

<table><tbody>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>int</code>. This MUST be the total number of read pairs associated with this fragment. The value MUST be <code>1</code> or greater.
        </td>
    </tr>
</tbody></table>
<br>

## scATAC-seq Asset: Processed Fragments File

From every submitted fragments file asset, CELLxGENE Discover MUST generate <code>{artifact_id}-fragments.tsv.gz</code>, a tab-separated values (TSV) file position-sorted and compressed by bgzip.

## scATAC-seq Asset: Fragments File index

From every processed fragments file asset, CELLxGENE Discover MUST generate <code>{artifact_id}-fragments.tsv.gz.tbi</code>, a <a href="https://www.htslib.org/doc/tabix.html">tabix</a> index of the fragment intervals from the fragments file.

## Chromosome Tables

Chromosome Tables are determined by the reference assembly for the gene annotation versions pinned in this version of the schema. Only chromosomes or scaffolds that have at least one gene feature present are included.

### <a href="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz">Human (GRCh38.p14)</a>

<table>
  <thead>
  <tr>
  <th>Chromosome</th>
  <th>Length</th>
  </tr>
  </thead>
  <tbody>
    <tr>
        <td>chr1</td>
        <td>248956422</td>
    </tr>
    <tr>
        <td>chr2</td>
        <td>242193529</td>
    </tr>
    <tr>
        <td>chr3</td>
        <td>198295559</td>
    </tr>
    <tr>
        <td>chr4</td>
        <td>190214555</td>
    </tr>
    <tr>
        <td>chr5</td>
        <td>181538259</td>
    </tr>
    <tr>
        <td>chr6</td>
        <td>170805979</td>
    </tr>
    <tr>
        <td>chr7</td>
        <td>159345973</td>
    </tr>
    <tr>
        <td>chr8</td>
        <td>145138636</td>
    </tr>
    <tr>
        <td>chr9</td>
        <td>138394717</td>
    </tr>
    <tr>
        <td>chr10</td>
        <td>133797422</td>
    </tr>
    <tr>
        <td>chr11</td>
        <td>135086622</td>
    </tr>
    <tr>
        <td>chr12</td>
        <td>133275309</td>
    </tr>
    <tr>
        <td>chr13</td>
        <td>114364328</td>
    </tr>
    <tr>
        <td>chr14</td>
        <td>107043718</td>
    </tr>
    <tr>
        <td>chr15</td>
        <td>101991189</td>
    </tr>
    <tr>
        <td>chr16</td>
        <td>90338345</td>
    </tr>
    <tr>
        <td>chr17</td>
        <td>83257441</td>
    </tr>
    <tr>
        <td>chr18</td>
        <td>80373285</td>
    </tr>
    <tr>
        <td>chr19</td>
        <td>58617616</td>
    </tr>
    <tr>
        <td>chr20</td>
        <td>64444167</td>
    </tr>
    <tr>
        <td>chr21</td>
        <td>46709983</td>
    </tr>
    <tr>
        <td>chr22</td>
        <td>50818468</td>
    </tr>
    <tr>
        <td>chrX</td>
        <td>156040895</td>
    </tr>
    <tr>
        <td>chrY</td>
        <td>57227415</td>
    </tr>
    <tr>
        <td>chrM</td>
        <td>16569</td>
    </tr>
    <tr>
        <td>GL000009.2</td>
        <td>201709</td>
    </tr>
    <tr>
        <td>GL000194.1</td>
        <td>191469</td>
    </tr>
    <tr>
        <td>GL000195.1</td>
        <td>182896</td>
    </tr>
    <tr>
        <td>GL000205.2</td>
        <td>185591</td>
    </tr>
    <tr>
        <td>GL000213.1</td>
        <td>164239</td>
    </tr>
    <tr>
        <td>GL000216.2</td>
        <td>176608</td>
    </tr>
    <tr>
        <td>GL000218.1</td>
        <td>161147</td>
    </tr>
    <tr>
        <td>GL000219.1</td>
        <td>179198</td>
    </tr>
    <tr>
        <td>GL000220.1</td>
        <td>161802</td>
    </tr>
    <tr>
        <td>GL000225.1</td>
        <td>211173</td>
    </tr>
    <tr>
        <td>KI270442.1</td>
        <td>392061</td>
    </tr>
    <tr>
        <td>KI270711.1</td>
        <td>42210</td>
    </tr>
    <tr>
        <td>KI270713.1</td>
        <td>40745</td>
    </tr>
    <tr>
        <td>KI270721.1</td>
        <td>100316</td>
    </tr>
    <tr>
        <td>KI270726.1</td>
        <td>43739</td>
    </tr>
    <tr>
        <td>KI270727.1</td>
        <td>448248</td>
    </tr>
    <tr>
        <td>KI270728.1</td>
        <td>1872759</td>
    </tr>
    <tr>
        <td>KI270731.1</td>
        <td>150754</td>
    </tr>
    <tr>
        <td>KI270733.1</td>
        <td>179772</td>
    </tr>
    <tr>
        <td>KI270734.1</td>
        <td>165050</td>
    </tr>
    <tr>
        <td>KI270744.1</td>
        <td>168472</td>
    </tr>
    <tr>
        <td>KI270750.1</td>
        <td>148850</td>
    </tr>
</tbody></table>

### <a href="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/GRCm39.primary_assembly.genome.fa.gz">Mouse (GRCm39)</a>

<table>
  <thead>
  <tr>
  <th>Chromosome</th>
  <th>Length</th>
  </tr>
  </thead>
  <tbody>
    <tr>
        <td>chr1</td>
        <td>195154279</td>
    </tr>
    <tr>
        <td>chr2</td>
        <td>181755017</td>
    </tr>
    <tr>
        <td>chr3</td>
        <td>159745316</td>
    </tr>
    <tr>
        <td>chr4</td>
        <td>156860686</td>
    </tr>
    <tr>
        <td>chr5</td>
        <td>151758149</td>
    </tr>
    <tr>
        <td>chr6</td>
        <td>149588044</td>
    </tr>
    <tr>
        <td>chr7</td>
        <td>144995196</td>
    </tr>
    <tr>
        <td>chr8</td>
        <td>130127694</td>
    </tr>
    <tr>
        <td>chr9</td>
        <td>124359700</td>
    </tr>
    <tr>
        <td>chr10</td>
        <td>130530862</td>
    </tr>
    <tr>
        <td>chr11</td>
        <td>121973369</td>
    </tr>
    <tr>
        <td>chr12</td>
        <td>120092757</td>
    </tr>
    <tr>
        <td>chr13</td>
        <td>120883175</td>
    </tr>
    <tr>
        <td>chr14</td>
        <td>125139656</td>
    </tr>
    <tr>
        <td>chr15</td>
        <td>104073951</td>
    </tr>
    <tr>
        <td>chr16</td>
        <td>98008968</td>
    </tr>
    <tr>
        <td>chr17</td>
        <td>95294699</td>
    </tr>
    <tr>
        <td>chr18</td>
        <td>90720763</td>
    </tr>
    <tr>
        <td>chr19</td>
        <td>61420004</td>
    </tr>
    <tr>
        <td>chrX</td>
        <td>169476592</td>
    </tr>
    <tr>
        <td>chrY</td>
        <td>91455967</td>
    </tr>
    <tr>
        <td>chrM</td>
        <td>16299</td>
    </tr>
    <tr>
        <td>GL456210.1</td>
        <td>169725</td>
    </tr>
    <tr>
        <td>GL456211.1</td>
        <td>241735</td>
    </tr>
    <tr>
        <td>GL456212.1</td>
        <td>153618</td>
    </tr>
    <tr>
        <td>GL456219.1</td>
        <td>175968</td>
    </tr>
    <tr>
        <td>GL456221.1</td>
        <td>206961</td>
    </tr>
    <tr>
        <td>GL456239.1</td>
        <td>40056</td>
    </tr>
    <tr>
        <td>GL456354.1</td>
        <td>195993</td>
    </tr>
    <tr>
        <td>GL456372.1</td>
        <td>28664</td>
    </tr>
    <tr>
        <td>GL456381.1</td>
        <td>25871</td>
    </tr>
    <tr>
        <td>GL456385.1</td>
        <td>35240</td>
    </tr>
    <tr>
        <td>JH584295.1</td>
        <td>1976</td>
    </tr>
    <tr>
        <td>JH584296.1</td>
        <td>199368</td>
    </tr>
    <tr>
        <td>JH584297.1</td>
        <td>205776</td>
    </tr>
    <tr>
        <td>JH584298.1</td>
        <td>184189</td>
    </tr>
    <tr>
        <td>JH584299.1</td>
        <td>953012</td>
    </tr>
    <tr>
        <td>JH584303.1</td>
        <td>158099</td>
    </tr>
    <tr>
        <td>JH584304.1</td>
        <td>114452</td>
    </tr>
</tbody></table>

---

## Appendix A. Changelog

### schema v7.0.0
* General Requirements
  * Integration Metadata
    * Updated the requirements for prefixed ontology identifiers to address the Cellosaurus exception
* Required Ontologies
  * Added Cellosaurus release 52.0
  * Updated HANCESTRO to release 2025-04-01
* Required Gene Annotations
  * Updated *Caenorhabditis elegans* to WBcel235 (GCA_000002985.3) Ensembl 114
  * Updated *Callithrix jacchus* to mCalJac1.pat.X (GCA_011100555.1) Ensembl 114
  * Updated *Danio rerio* to GRCz11 (GCA_000002035.4) Ensembl 114
  * Updated *Drosophila melanogaster* to BDGP6.54 (GCA_000001215.4) Ensembl 114
  * Updated *Gorilla gorilla gorilla* to gorGor4 (GCA_000151905.3) Ensembl 114
  * Updated *Homo sapiens* to GENCODE v48 (GRCh38.p14) Ensembl 114
  * Updated *Macaca fascicularis* to Macaca_fascicularis_6. (GCA_011100615.1) Ensembl 114
  * Updated *Macaca mulatta* to Mmul_10 (GCA_003339765.3) Ensembl 114
  * Updated *Microcebus murinus* to Mmur_3.0 (GCA_000165445.3) Ensembl 114
  * Updated *Mus musculus* to GENCODE vM37 (GRCm39) Ensembl 114
  * Updated *Oryctolagus cuniculus* to OryCun2.0 (GCA_000003625.1) Ensembl 114
  * Updated *Pan troglodytes* to Pan_tro_3.0 (GCA_000001515.5) Ensembl 114
  * Updated *Rattus norvegicus* to GRCr8 (GCA_036323735.1) Ensembl 114
  * Updated *Sus scrofa* to Sscrofa11.1 (GCA_000003025.6) Ensembl 114
* obs (Cell metadata)
  * Updated the requirements for <code>cell_type</code> to require <code>"na"</code> when the <code>cell_type_ontology_term_id</code> is <code>"na"</code>
  * Updated the requirements for <code>cell_type_ontology_term_id</code> to allow <code>"na"</code> when the <code>tissue_type</code> is <code>"cell line"</code>
  * Updated the requirements for <code>development_stage</code> to require <code>"na"</code> when the <code>development_stage__ontology_term_id</code> is <code>"na"</code>
  * Updated the requirements for <code>development_stage_ontology_term_id</code> to require <code>"na"</code> when the <code>tissue_type</code> is <code>"cell line"</code>
  * Updated the requirements for <code>donor_id</code> to require <code>"na"</code> when the <code>tissue_type</code> is <code>"cell line"</code>
  * Updated the requirements for <code>self_reported_ethnicity_ontology_term_id</code> to require <code>"na"</code> when the <code>tissue_type</code> is <code>"cell line"</code>
  * Updated the requirements for <code>self_reported_ethnicity_ontology_term_id</code> to require HANCESTRO or AfPO terms that are descendants of <code>"HANCESTRO:0601"</code> for <i>ethnicity category</i> or <code>"HANCESTRO:0602"</code> for <i>geography-based population category</i>
  * Updated the requirements for <code>sex</code> to require <code>"na"</code> when the <code>sex_ontology_term_id</code> is <code>"na"</code>
  * Updated the requirements for <code>sex_ontology_term_id</code> to require <code>"na"</code> when the <code>tissue_type</code> is <code>"cell line"</code>
  * Updated the requirements for <code>suspension_type</code>:
    * Added Cel-seq
    * Added Quartz-seq
    * Deleted sci-RNA-seq
    * Updated CEL-seq2 to CEL-seq2 and its descendants
    * Updated SPLiT-seq to SPLiT-seq and its descendants
    * Updated STRT-seq to STRT-seq and its descendants
  * **Breaking change**. Updated the requirements for <code>tissue_ontology_term_id</code> to rename the <code>tissue_type</code> of <code>"cell culture"</code> to <code>"primary cell culture"</code>
  * Updated the requirements for <code>tissue_ontology_term_id</code> to add the <code>tissue_type</code> of <code>"cell line"</code>
  * Updated the requirements for <code>tissue_ontology_term_id</code> when the <code>tissue_type</code> is <code>"organoid"</code> 
  * Updated the requirements for <code>tissue_ontology_term_id</code> for species with taxon specific ontologies to require the most accurate descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/uberon/classes?obo_id=UBERON%3A0001062"><code>UBERON:0001062</code></a> for <i>anatomical entity</i> when the <code>tissue_type</code> is <code>"tissue"</code> or <code>"organoid"</code>
  * **Breaking change**. Updated the requirements for <code>tissue_type</code> to rename <code>"cell culture"</code> to <code>"primary cell culture"</code>
  * Added <code>"cell line"</code> to <code>tissue_type</code>
* uns (Dataset Metadata)
  * Updated `schema_reference` to <code>"https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/7.0.0/schema.md"</code>
  * Updated `schema_version` to <code>"7.0.0"</code>

### schema v6.0.0

* General Requirements
  * Added <code>organism</code> in <code>obs</code> to Reserved Names
  * Added <code>organism_ontology_term_id</code> in <code>obs</code> to Reserved Names
* Required Ontologies
  * Updated C. elegans Development Ontology (WBls) to release 2025-04-01 WS297
  * Updated C. elegans Gross Anatomy Ontology (WBbt) to release 2025-03-26 WS297
  * Updated Cell Ontology (CL) to release 2025-04-10
  * Updated Drosophila Anatomy Ontology (FBbt) to release 2025-03-27
  * Updated Drosophila Development Ontology (FBdv) to release 2025-03-26
  * Updated Expermental Factor Ontology (EFO) to release 2025-05-15 EFO 3.78.0
  * Updated Mondo Disease Ontology (MONDO) to release 2025-05-06
  * Updated NCBI organismal classification (NCBITaxon) to release 2025-03-13
  * Updated Phenotype And Trait Ontology (PATO) to release 2025-05-14
  * Updated  Uberon multi-species anatomy ontology (UBERON) to release 2025-05-28
* X (Matrix Layers)
  * Updated requirements in table to not allow duplicate `obs`by raw counts
* obs (Cell metadata)
  * Updated <code>disease</code> to address multiple labels
  * Updated <code>disease_ontology_term_id</code> to allow multiple terms
  * Updated <code>self_reported_ethnicity</code> delimiter requirements
  * Updated <code>self_reported_ethnicity_ontology_term_id</code> delimiter requirements
  * Deprecated <code>organism</code>
  * Deprecated <code>organism_ontology_term_id</code>
* uns (Dataset Metadata)
  * Updated `{column}_colors` requirements to not reference <code>organism</code> 
  * Added <code>organism</code>
  * Added <code>organism_ontology_term_id</code>
  * Updated `schema_reference` to <code>"https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/6.0.0/schema.md"</code>
  * Updated `schema_version` to <code>"6.0.0"</code>
* var and raw.var (Gene Metadata)
  * Updated <code>index of pandas.DataFrame</code> requirements to remove version suffixes from Ensembl stable identifiers
  * Clarified <code>feature_is_filtered</code> requirements
  * Updated <code>feature_name</code> requirements to use feature identifier if a <code>gene_name</code> is unavailable
* scATAC-seq assets
  * Updated requirements for *one organism per dataset*

### schema v5.3.0
* Schema Versioning
  * **Patch versions** may be used to add organisms that do not require new metadata fields. 
* General Requirements
* Required Ontologies
  * Added C. elegans Development Ontology (WBls) release 2025-01-04 WS296
  * Added C. elegans Gross Anatomy Ontology (WBbt) release 2025-01-02 WS296
  * Updated Cell Ontology (CL) to the 2025-02-13 release
  * Added Drosophila Anatomy Ontology (FBbt) release 2025-02-13
  * Added Drosophila Development Ontology (FBdv) release 2025-02-12
  * Updated Expermental Factor Ontology (EFO) to release 2025-02-17 EFO 3.75.0
  * Updated Human Developmental Stages (HsapDv) to release 2025-01-23
  * Updated Mondo Disease Ontology (MONDO) to release 2025-02-04
  * Updated Mouse Developmental Stages (MmusDv) to release 2025-01-2
  * Updated NCBI organismal classification (NCBITaxon) to release 2024-11-25
  * Updated Phenotype And Trait Ontology (PATO) to release 2025-02-01
  * Updated Uberon multi-species anatomy ontology (UBERON) to release 2025-01-15
  * Added Zebrafish Anatomy Ontology (ZFA+ZFS) release 2025-01-28
* Required Gene Annotations
  * Refactored table to include NCBI Taxon for supported organisms
  * Added *Caenorhabditis elegans* WBcel235 (GCA_000002985.3) Ensembl 113
  * Added *Callithrix jacchus* mCalJac1.pat.X (GCA_011100555.1) Ensembl 113
  * Added *Danio rerio* GRCz11 (GCA_000002035.4) Ensembl 113
  * Added *Drosophila melanogaster* BDGP6.46 (GCA_000001215.4) Ensembl 113
  * Added *Gorilla gorilla gorilla* gorGor4 (GCA_000151905.3) Ensembl 113
  * Added *Macaca fascicularis* Macaca_fascicularis_6.0 (GCA_011100615.1) Ensembl 113
  * Added *Macaca mulatta* Mmul_10 (GCA_003339765.3) Ensembl 113
  * Added *Microcebus murinus* Mmur_3.0 (GCA_000165445.3) Ensembl 113
  * Added *Oryctolagus cuniculus* OryCun2.0 (GCA_000003625.1) Ensembl 113
  * Added *Pan troglodytes* Pan_tro_3.0 (GCA_000001515.5) Ensembl 113
  * Added *Rattus norvegicus* mRatBN7.2 (GCA_015227675.2) Ensembl 113
  * Added *Sus scrofa* Sscrofa11.1 (GCA_000003025.6) Ensembl 113
* X (Matrix Layers)
  * Updated _Visium Spatial Gene Expression_ table row to _Descendants of Visium Spatial Gene Expression_
  * Added matrix requirements for _Visium CytAssist Spatial Gene Expression, 11mm_.
  * Updated the STRONGLY RECOMMENDED requirement to a MUST. A matrix with 50% or more values that are zeros MUST be encoded as `scipy.sparse.csr_matrix`.
* obs (Cell metadata)
  * Updated the requirements for `array_col`:
    * MUST be annotated if the `assay_ontology_term_id` is a descendant of _Visium Spatial Gene Expression_
    * Added ranges for _Visium CytAssist Spatial Gene Expression, 6.5mm_ and _Visium CytAssist Spatial Gene Expression, 11mm_ 
  * Updated the requirements for `array_row`:
     * MUST be annotated if the `assay_ontology_term_id` is a descendant of _Visium Spatial Gene Expression_
    * Added ranges for _Visium CytAssist Spatial Gene Expression, 6.5mm_ and _Visium CytAssist Spatial Gene Expression, 11mm_ 
  * Updated the requirements for `assay_ontology_term_id`:
    * For _Visium Spatial Gene Expression_, only its descendants are allowed. All observations must contain the same value.
    * For _scATAC-seq and its descendants_, added reference to scATAC-seq fragment asset requirements 
    * Added more recommended terms for assays
  * Updated the requirements for `cell_type_ontology_term_id` to include descendants of  _Visium Spatial Gene Expression_.
  * Updated the ontology requirements for `cell_type_ontology_term_id` to include:
    * WBbt for *Caenorhabditis elegans*
    * ZFA for *Danio rerio*
    * FBbt for *Drosophila melanogaster* 
  * Updated the ontology requirements for `development_stage_ontology_term_id` to include:
    * WBls for *Caenorhabditis elegans*
    * ZFS for *Danio rerio*
    * FBdv for *Drosophila melanogaster*
  * Updated `development_stage_ontology_term_id` to include descendants of <i>Mus musclus</i>
  * Updated the requirements for `in_tissue` to include descendants of  _Visium Spatial Gene Expression_.
  * Updated the requirements for `organism_ontology_term_id` to limit its values to a list of NCBITaxon terms
  * Updated the requirements for `sex_ontology_term_id`:
    * For *Caenorhabditis elegans*, values are limited to <i>hermaphrodite</i>, <i>male</i>, or `"unknown"`
    * For all other organisms, values are limited to <i>female</i>, <i>hermaphrodite</i>, <i>male</i>, or `"unknown"`
  * Updated the ontology requirements for `tissue_ontology_term_id` to include:
    * WBbt for *Caenorhabditis elegans*
    * ZFA for *Danio rerio*
    * FBbt for *Drosophila melanogaster* 
* obsm (Embeddings)
  * Updated the requirements for `spatial` to include descendants of  _Visium Spatial Gene Expression_ and to prohibit 'Not a Number' values. 
  * Updated the requirements for `X_{suffix}` to include descendants of  _Visium Spatial Gene Expression_.
* uns (Dataset Metadata)
  * Updated `schema_reference` to <code>"https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.3.0/schema.md"</code>
  * Updated `schema_version` to <code>"5.3.0"</code>
  * Updated the requirements for `spatial` to include descendants of  _Visium Spatial Gene Expression_.
  * Updated the requirements for  <code>spatial['is_single']</code> to include descendants of  _Visium Spatial Gene Expression_.
  * Updated the requirements for  <code>spatial[<i>library_id</i>]</code> to include descendants of  _Visium Spatial Gene Expression_.
  * Updated the requirements for  <code>spatial[<i>library_id</i>]['images']</code> to include descendants of  _Visium Spatial Gene Expression_.
  * Updated the requirements for  <code>spatial[<i>library_id</i>]['images']['fullres']</code> to include descendants of  _Visium Spatial Gene Expression_.
  * Updated the requirements for  <code>spatial[<i>library_id</i>]['images']['hires']</code> to include descendants of  _Visium Spatial Gene Expression_. Added requirements for <i>Visium CytAssist Spatial Gene Expression, 11mm</i>.
  * Updated the requirements for  <code>spatial[<i>library_id</i>]['scalefactors']</code> to include descendants of  _Visium Spatial Gene Expression_.
  * Updated the requirements for  <code>spatial[<i>library_id</i>]['scalefactors']['spot_diameter_fullres']</code> to include descendants of  _Visium Spatial Gene Expression_.
  * Updated the requirements for  <code>spatial[<i>library_id</i>]['scalefactors']['tissue_hires_scalef']</code> to include descendants of  _Visium Spatial Gene Expression_.
* var and raw.var (Gene Metadata)
  * Updated `feature_reference` to include MUST requirements for:
    * *Caenorhabditis elegans*
    * *Danio rerio*
    * *Drosophila melanogaster*
    * *Macaca fascicularis*
    * *Microcebus murinus*
    * *Oryctolagus cuniculus*
    * *Rattus norvegicus*
  * Updated `feature_reference` to include STRONGLY RECOMMENDED requirements and warning for current organisms with orthologous gene references:
    * *Callithrix jacchus*
    * *Gorilla gorilla gorilla*
    * *Macaca mulatta*
    * *Pan troglodytes*
    * *Sus scrofa and its descendants*
* Added **scATAC-seq assets** section


### schema v5.2.0

* General Requirements
  * Updated AnnData from version 0.8.0 to version 0.8.0 or greater
* Required Ontologies
  * Updated CL to the 2024-08-16 release
  * Updated EFO to the 2024-08-15 EFO 3.69.0 release
  * Updated HsapDv to the 2024-05-28 release
  * Updated MONDO to the 2024-08-06 release
  * Updated MmusDv to the 2024-05-28 release
  * Updated UBERON to the 2024-08-07 release
* obs (Cell metadata)
  * Updated requirements for `development_stage_ontology_term_id` to require the most accurate descendant of _life cycle_. 
    * If <code>organism_ontology_term_id</code> is <code>"NCBITaxon:9606"</code> for <i>Homo sapiens</i>, this MUST be the most accurate descendant of `HsapDv:0000001` for _life cycle_
    * If <code>organism_ontology_term_id</code> is <code>"NCBITaxon:10090"</code> for <i>Mus musculus</i>, this MUST be the most accurate descendant of `MmusDv:0000001` for _life cycle_
  * Updated requirements for `suspension_type`
    * Added mCT-seq
    * Added MERFISH
    * Added ScaleBio single cell RNA sequencing
    * Added sci-RNA-seq3
    * Removed CITE-seq and its descendants
    * Removed smFISH and its descendants
    * Removed snmC-seq
    * Removed spatial proteomics and its descendants
    * Replaced snmC-seq2 with methylation profiling by high throughput sequencing and its descendants
* uns (Dataset metadata)
  * Updated `schema_reference` to <code>"https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.2.0/schema.md"</code>
  * Updated `schema_version` to <code>"5.2.0"</code>
* var and raw.var (Gene metadata)
  * Updated the requirements for `feature_length`. All `feature_biotypes` are now included. The calculation of the value changed from the merged length of isoforms to the median of the lengths of isoforms.
  * Added `feature_type`

### schema v5.1.0

* All references to "child" and "children" have been changed to "descendant" and "descendants" for accuracy.
* Required Ontologies
  * Updated CL to the 2024-04-05 release
  * Updated EFO to the 2024-04-15 EFO 3.65.0 release
  * Updated MONDO to the 2024-05-08 release
  * Updated UBERON to the 2024-03-22 release
* X (Matrix Layers)
  * Added _Visium Spatial Gene Expression_ to the table of assays
* obs (Cell metadata)
  * Added `array_col` for _Visium Spatial Gene Expression_ when <code>uns['spatial']['is_single']</code> is <code>True</code>
  * Added `array_row` for _Visium Spatial Gene Expression_ when <code>uns['spatial']['is_single']</code> is <code>True</code>
  * Updated the requirements for `assay_ontology_term_id` for _Visium Spatial Gene Expression_ and _Slide-seqV2_. All observations must contain the same value.
  * Updated the requirements for `cell_type_ontology_term_id` for _Visium Spatial Gene Expression_ when <code>uns['spatial']['is_single']</code> is <code>True</code>. The value must be `"unknown"` if the corresponding value of `in_tissue` is `0`.
  * Added `in_tissue` for _Visium Spatial Gene Expression_ when <code>uns['spatial']['is_single']</code> is <code>True</code>
  * Updated the requirements for `is_primary_data` for _Visium Spatial Gene Expression_. The value must be <code>False</code>when <code>uns['spatial']['is_single']</code> is <code>False</code>.
  * Updated the requirements for `self_reported_ethnicity_ontology_term_id`. There must be no duplication of terms.
* obsm (Embeddings)
  * Restored v3.1.0 requirement allowing only `numpy.ndarray` values with specific shapes due to Seurat conversion failures
  * Added `spatial` for _Visium Spatial Gene Expression_ and _Slide-seqV2_
  * Updated requirements for `X_{suffix}`. {suffix} MUST NOT be `"spatial"`.
* uns (Dataset metadata)
  * Updated `{column}_colors` instructions
  * Updated `schema_reference` to <code>"https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md"</code>
  * Updated `schema_version` to <code>"5.1.0"</code>
  * Added `spatial` for _Visium Spatial Gene Expression_ and _Slide-seqV2_, including scale factors and underlay images for _Visium Spatial Gene Expression_.

### schema v5.0.0

* General Requirements
  * Updated requirements to prohibit duplicate data submitter metadata field names in `obs` and `var`
* Required Ontologies
  * Updated CL to the 2024-01-04 release
  * Updated EFO to the 2024-01-15 EFO 3.62.0 release
  * Updated MONDO to the 2024-01-03 release
  * Updated UBERON to the 2024-01-18 release
* Required Gene Annotations
  * Updated GENCODE (Human) to Human Reference GRCh38.p14 (GENCODE v44/Ensembl 110)
  * Updated GENCODE (Mouse) to Mouse reference GRCm39 (GENCODE vM33/Ensembl 110)
* obs (Cell metadata)
  * Updated the requirements for `assay_ontology_term_id` to not allow the parent terms `EFO:0002772` for _assay by molecule_ and `EFO:0010183` for _single cell library construction_. Their most accurate children are still valid. 
  * **Breaking change**. Updated the requirements for `cell_type` to annotate `"unknown"` as the label when the `cell_type_ontology_term_id` value is  `"unknown"`. 
  * **Breaking change**. Updated the requirements for `cell_type_ontology_term_id` to replace `"CL:0000003"` for *native cell* with `"unknown"` to indicate that the cell type is unknown. 
  * Updated the requirements for `disease_ontology_term_id` to restrict MONDO terms to the most accurate child of `"MONDO:0000001"` for _disease_ or `"MONDO:0021178"` for _injury_ or preferably its most accurate child.
* obsm (Embeddings)
  * Updated requirements for `X_{suffix}` to change the regular expression pattern from `"^[a-zA-Z][a-zA-Z0-9]*$"` to `"^[a-zA-Z][a-zA-Z0-9_.-]*$"`
* uns (Dataset metadata)
  * Updated requirements. The data stored as a value for a key in `uns` MUST be `True`, `False`, `None`, or its size MUST NOT be zero.
  * Updated schema_reference to <code>"https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.0.0/schema.md"</code>
  * Updated schema_version to <code>"5.0.0"</code>

### schema v4.0.0

* Required Ontologies
  * Updated CL to the 2023-08-24 release
  * Updated EFO to the 2023-08-15 EFO 3.57.0 release
  * Updated HANCESTRO to the 3.0 release
  * Updated MONDO to the 2023-08-02 release
  * Updated UBERON to the 2023-09-05 release
* obs (Cell metadata)
  * Updated the requirements for `cell_type_ontology_term_id`
  * Added `index`
  * Added `observation_joinid`
  * Updated the requirements for `self_reported_ethnicity`
  * Updated the requirements for `self_reported_ethnicity_ontology_term_id`
  * Added `tissue_type`
  * Updated the requirements for `tissue`
  * Updated the requirements for `tissue_ontology_term_id`
* obsm (Embeddings)
  * Prohibited ndarrays with a size of zero
  * Updated requirements for `X_{suffix}`
* obsp
  * Added section and prohibited ndarrays with a size of zero
* uns (Dataset metadata)
  * Added `citation`
  * Added `{column}_colors`
  * Added `schema_reference`
  * Updated the requirements for `schema_version` to be annotated by CELLxGENE Discover and not the curator
* var and raw.var (Gene metadata)
  * Added `feature_length`
* varm
  * Added section and prohibited ndarrays with a size of zero
* varp
  * Added section and prohibited ndarrays with a size of zero
* X (Matrix Layers)
  * Updated requirements for raw matrices

### schema v3.1.0

* Added section for Schema versioning
* Required Ontologies
  * Updated CL to the 2023-07-20 release
  * Updated EFO to the 2023-07-17 EFO 3.56.0 release
  * Updated MONDO to the 2023-07-03 release
  * Updated NCBITaxon to the 2023-06-20 release
  * Updated PATO to the 2023-05-18 release
  * Updated UBERON to the 2023-06-28 release
* obs (Cell metadata)
  * `assay_ontology_term_id`
    * Added Visium Spatial Gene Expression to recommended values
    * Removed Smart-seq from recommended values
  * `suspension_type`
    * Added MARS-seq
    * Added BD Rhapsody Whole Transcriptome Analysis
    * Added BD Rhapsody Targeted mRNA
    * Added inDrop
    * Added STRT-seq
    * Added TruDrop
    * Added GEXSCOPE technology
    * Added SPLiT-seq
    * Changed spatial transcriptomics by high-throughput sequencing [EFO:0030005] and its children to spatial transcriptomics [EFO:0008994] and its children
    * Updated Seq-Well [EFO:0008919] to Seq-Well [EFO:0008919] and its children
* uns (Dataset metadata)
  * `schema_version`
    * Must be annotated by CELLxGENE Discover and not the Curator.


### schema v3.0.0

* Updated AnnData version 0.7 to version 0.8.0
* All references to the "final" matrix has been replaced with "normalized" for clarity.
* General Requirements
  * Reserved Names from previous schema versions that have since been deprecated MUST NOT be present.
  * Updated *pinned* ontologies to require the most recent version
* obs (Cell metadata)
  * Removed guidance in `assay_ontology_term_id` that allowed clarifying text enclosed in parentheses if there was not an exact match for an assay.
  * Added `donor_id`
  * Renamed `ethnicity_ontology_term_id` to `self_reported_ethnicity_ontology_term_id`. Added `"multiethnic"` value.
  * Renamed `ethnicity` to `self_reported_ethnicity`. Added `"multiethnic"` value.
  * Added `suspension_type`
* var and raw.var (Gene metadata)
  * `feature_biotype` must be annotated by CELLxGENE Discover and not the Curator.
* uns (Dataset metadata)
  * Updated `schema_version`
  * Deprecated `X_normalization`

### schema v2.0.0

schema v2.0.0 substantially *remodeled* schema v1.1.0:

* "must", "should", and select other words have a defined, standard meaning.

* Curators are responsible for annotating ontology and gene identifiers. CELLxGENE Discover adds the assigned human-readable names for all identifiers.

* Documented and *pinned* the required versions of ontologies and gene annotations used in schema validation.

* General Requirements
  * AnnData is now the canonical data format. The schema outline and descriptions are AnnData-centric.

  * Metazoan multi-organism data is accepted by CELLxGene Discover. For data that is neither Human, Mouse, nor SARS-COV-2, features MUST be translated into orthologous genes from the Human and Mouse gene annotations. 

  * Policies for reserved names and redundant metadata are documented.

  * [#45](https://github.com/chanzuckerberg/single-cell-curation/issues/45) Updated reference to new PII content

* X (matrix layers)
  * Added guidance for sparse matrices
  * Clarified matrix requirements by assay

* obs (cell metadata)
  * Empty ontology fields are no longer permitted.
  * Moved organism from uns to obs
  * Clarified requirements and added detailed guidance for assays, tissue, and development stages
  * Added ontology for mouse development stages
  * Added ontology for sex
  * Added `is_primary_data`

* var
  * Replaced HGNC gene **symbols** as `var.index` with ENSEMBL or ERCC spike-in **identifiers** 
  * Added `feature_name`, `index`, and `feature_reference`
  * Added `feature_is_filtered`
  * Added requirements for `raw.var` which must be identical to `var`

* uns
  * Added `batch_condition`
  * Added `X_approximate_distribution`
  * Replaced `layer_descriptions` with `X_normalization`
  * Replaced `version` which included `corpora_schema_version` and `corpora_encoding_version` with `schema_version`
  * Deprecated `tags` and `default_field` presentation metadata
  * Removed <code><i>obs_column</i>_colors</code>
