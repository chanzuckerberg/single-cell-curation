
# Schema

Contact: brianraymor@chanzuckerberg.com

Document Status: _Approved_

Version: 5.0.0

The key words "MUST", "MUST NOT", "REQUIRED", "SHALL", "SHALL NOT", "SHOULD", "SHOULD NOT", "RECOMMENDED", "NOT RECOMMENDED" "MAY", and "OPTIONAL" in this document are to be interpreted as described in [BCP 14](https://tools.ietf.org/html/bcp14), [RFC2119](https://www.rfc-editor.org/rfc/rfc2119.txt), and [RFC8174](https://www.rfc-editor.org/rfc/rfc8174.txt) when, and only when, they appear in all capitals, as shown here.

## Schema versioning

The CELLxGENE schema version is based on [Semantic Versioning](https://semver.org/).

**Major version** is incremented when schema updates are incompatible with the AnnData and Seurat data encodings or CELLxGENE API(s). Examples include:
  * Renaming metadata fields
  * Deprecating metadata fields
  * Changing the type or format of a metadata field
 
**Minor version** is incremented when schema updates may require changes only to the `cellxgene-schema` CLI or the curation process. Examples include:
  * Adding metadata fields
  * Updating pinned ontologies or gene references
  * Changing the validation requirements for a metadata field
  
**Patch version** is incremented for editorial updates to the schema.

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

**AnnData.** The canonical data format for CELLxGENE Discover is HDF5-backed [AnnData](https://anndata.readthedocs.io/en/latest) as written by version 0.8 of the anndata library.  Part of the rationale for selecting this format is to allow CELLxGENE to access both the data and metadata within a single file. The schema requirements and definitions for the AnnData `X`, `obs`, `var`, `raw.var`, `obsm`, and `uns` attributes are described below.

All data submitted to CELLxGENE Discover is automatically converted to a Seurat V5 object that can be loaded by the R package Seurat. See the [Seurat encoding](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.0.0/seurat_encoding.md) for further information.

**Organisms**. Data MUST be from a Metazoan organism or SARS-COV-2 and defined in the NCBI organismal classification. For data that is neither Human, Mouse, nor SARS-COV-2, features MUST be translated into orthologous genes from the pinned Human and Mouse gene annotations.

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
The types below are python3 types. Note that a python3 `str` is a sequence of Unicode code points, which is stored null-terminated and UTF-8-encoded by anndata.

## `X` (Matrix Layers)

The data stored in the `X` data matrix is the data that is viewable in CELLxGENE Explorer. CELLxGENE does not impose any additional constraints on the `X` data matrix.

In any layer, if a matrix has 50% or more values that are zeros, it is STRONGLY RECOMMENDED that the matrix be encoded as a [`scipy.sparse.csr_matrix`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html).

CELLxGENE's matrix layer requirements are tailored to optimize data reuse. Because each assay has different characteristics, the requirements differ by assay type. In general, CELLxGENE requires submission of "raw" data suitable for computational reuse when a standard raw matrix format exists for an assay. It is STRONGLY RECOMMENDED to also include a "normalized" matrix with processed values ready for data analysis and suitable for visualization in CELLxGENE Explorer. So that CELLxGENE's data can be provided in download formats suitable for both R and Python, the schema imposes the following requirements:

*   All matrix layers MUST have the same shape, and have the same cell labels and gene labels.
*   Because it is impractical to retain all barcodes in raw and normalized matrices, any cell filtering MUST be applied to both.
    By contrast, those wishing to reuse datasets require access to raw gene expression values, so genes SHOULD NOT be filtered from either dataset.
    Summarizing, any cell barcodes that are removed from the data MUST be filtered from both raw and normalized matrices and genes SHOULD NOT be filtered from the raw matrix.
*   Any genes that publishers wish to filter from the normalized matrix MAY have their values replaced by zeros and MUST be flagged in the column [`feature_is_filtered`](#feature_is_filtered) of [`var`](#var-and-rawvar-gene-metadata), which will mask them from exploration.
*   Additional layers provided at author discretion MAY be stored using author-selected keys, but MUST have the same cells and genes as other layers. It is STRONGLY RECOMMENDED that these layers have names that accurately summarize what the numbers in the layer represent (e.g. `"counts_per_million"`, `"SCTransform_normalized"`, or `"RNA_velocity_unspliced"`).

The following table describes the matrix data and layers requirements that are **assay-specific**. If an entry in the table is empty, the schema does not have any other requirements on data in those layers beyond the ones listed above.

| Assay | "raw" required? | "raw" location | "normalized" required? | "normalized" location |
|-|-|-|-|-|
| scRNA-seq (UMI, e.g. 10x v3) | REQUIRED. Values MUST be de-duplicated molecule counts. Each cell MUST contain at least one non-zero value. All non-zero values MUST be positive integers stored as `numpy.float32`.| `AnnData.raw.X` unless no "normalized" is provided, then `AnnData.X` | STRONGLY RECOMMENDED | `AnnData.X` |
| scRNA-seq (non-UMI, e.g. SS2) | REQUIRED. Values MUST be one of read counts (e.g. FeatureCounts) or  estimated fragments (e.g. output of RSEM). Each cell MUST contain at least one non-zero value. All non-zero values MUST be positive integers stored as `numpy.float32`. | `AnnData.raw.X` unless no "normalized" is provided, then `AnnData.X` | STRONGLY RECOMMENDED | `AnnData.X` |
| Accessibility (e.g. ATAC-seq, mC-seq) | NOT REQUIRED | | REQUIRED | `AnnData.X` | STRONGLY RECOMMENDED |
|||||

## Integration Metadata

CELLxGENE requires ontology terms to enable search, comparison, and integration of data.
Ontology terms for cell metadata MUST use [OBO-format identifiers](http://www.obofoundry.org/id-policy.html), meaning a CURIE (prefixed identifier) of the form **Ontology:Identifier**.
For example, [EFO:0000001](https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0000001) is a term in the Experimental Factor Ontology (EFO).

The most accurate ontology term MUST always be used. If an exact or approximate ontology term is not available, a new term may be requested:

- For the [Cell Ontology], data submitters may [suggest a new term](https://github.com/obophenotype/cell-ontology/issues/new?assignees=bvarner-ebi&labels=new+term+request%2C+cellxgene&template=a_adding_term_cellxgene.md&title=%5BNTR-cxg%5D) and [notify the curation team](mailto:cellxgene@chanzuckerberg.com) of the pending term request, so that the datasets can be updated once the term is available.

  To meet CELLxGENE schema requirements, the most accurate available CL term MUST be used until the new term is available. For example if `cell_type_ontology_term_id` describes a relay interneuron, but the most accurate available term in the CL ontology is [CL:0000099](https://www.ebi.ac.uk/ols4/ontologies/cl/classes?obo_id=CL%3A0000099) for *interneuron*, then the interneuron term can be used to fulfill this requirement and ensures that users searching for "neuron" are able to find these data.  If no appropriate term can be found (e.g. the cell type is unknown), then `"unknown"` MUST be used. Users will still be able to access more specific cell type annotations that have been submitted with the dataset (but aren't required by the schema).

   
- For all other ontologies, data submitters may submit a [request to the curation team](mailto:cellxgene@chanzuckerberg.com) during the submission process.

Terms documented as obsolete in an ontology MUST NOT be used. For example, [EFO:0009310](https://www.ebi.ac.uk/ols4/ontologies/efo/classes/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0009310) for *obsolete_10x v2* was marked as obsolete in EFO version 3.31.0 and replaced by [EFO:0009899](https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0009899) for *10x 3' v2*.

### Required Ontologies

The following ontology dependencies are *pinned* for this version of the schema.

| Ontology | OBO Prefix | Release | Download |
|:--|:--|:--|:--|
| [Cell Ontology] | CL |  [2024-01-04] | [cl.owl]|
| [Experimental Factor Ontology] | EFO | [2024-01-15 EFO 3.62.0] | [efo.owl]
| [Human Ancestry Ontology] | HANCESTRO | [3.0] | [hancestro-base.owl] |
| [Human Developmental Stages] |  HsapDv | 2020-03-10 | [hsapdv.owl] |
| [Mondo Disease Ontology] | MONDO | [2024-01-03] | [mondo.owl] |
| [Mouse Developmental Stages]| MmusDv | 2020-03-10 | [mmusdv.owl] |
| [NCBI organismal classification] |  NCBITaxon | [2023-06-20] | [ncbitaxon.owl] |
| [Phenotype And Trait Ontology] | PATO | [2023-05-18] | [pato.owl]  |
| [Uberon multi-species anatomy ontology] |  UBERON | [2024-01-18] | [uberon.owl] |
| | | | |

[Cell Ontology]: http://obofoundry.org/ontology/cl.html
[2024-01-04]: https://github.com/obophenotype/cell-ontology/releases/tag/v2024-01-04
[cl.owl]: https://github.com/obophenotype/cell-ontology/releases/download/v2024-01-04/cl.owl

[Experimental Factor Ontology]: http://www.ebi.ac.uk/efo
[2024-01-15 EFO 3.62.0]: https://github.com/EBISPOT/efo/releases/tag/v3.62.0
[efo.owl]: https://github.com/EBISPOT/efo/releases/download/v3.62.0/efo.owl

[Human Ancestry Ontology]: http://www.obofoundry.org/ontology/hancestro.html
[3.0]: https://github.com/EBISPOT/hancestro/releases/tag/3.0
[hancestro-base.owl]: https://github.com/EBISPOT/hancestro/blob/3.0/hancestro-base.owl

[Human Developmental Stages]: http://obofoundry.org/ontology/hsapdv.html
[hsapdv.owl]: http://purl.obolibrary.org/obo/hsapdv.owl

[Mondo Disease Ontology]: http://obofoundry.org/ontology/mondo.html
[2024-01-03]: https://github.com/monarch-initiative/mondo/releases/tag/v2024-01-03
[mondo.owl]: https://github.com/monarch-initiative/mondo/releases/download/v2024-01-03/mondo.owl

[Mouse Developmental Stages]: http://obofoundry.org/ontology/mmusdv.html
[mmusdv.owl]: http://purl.obolibrary.org/obo/mmusdv.owl

[NCBI organismal classification]: http://obofoundry.org/ontology/ncbitaxon.html
[2023-06-20]: https://github.com/obophenotype/ncbitaxon/releases/tag/v2023-06-20
[ncbitaxon.owl]: https://github.com/obophenotype/ncbitaxon/releases/download/v2023-06-20/ncbitaxon.owl.gz

[Phenotype And Trait Ontology]: http://www.obofoundry.org/ontology/pato.html
[2023-05-18]: https://github.com/pato-ontology/pato/releases/tag/v2023-05-18
[pato.owl]: https://github.com/pato-ontology/pato/blob/v2023-05-18/pato.owl

[Uberon multi-species anatomy ontology]: http://www.obofoundry.org/ontology/uberon.html
[2024-01-18]: https://github.com/obophenotype/uberon/releases/tag/v2024-01-18
[uberon.owl]: https://github.com/obophenotype/uberon/releases/download/v2024-01-18/uberon.owl

### Required Gene Annotations

ENSEMBL identifiers are required for genes and [External RNA Controls Consortium (ERCC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4978944/) identifiers for [RNA Spike-In Control Mixes] to ensure that all datasets measure the same features and can therefore be integrated.

The following gene annotation dependencies are *pinned* for this version of the schema. For multi-organism experiments, cells from any Metazoan organism are allowed as long as orthologs from the following organism annotations are used.

| Source | Required version | Download |
|:--|:--|:--|
| [GENCODE (Human)] | Human reference GRCh38.p14 (GENCODE v44/Ensembl 110) | [gencode.v44.primary_assembly.annotation.gtf] |
| [GENCODE (Mouse)] | Mouse reference GRCm39 (GENCODE vM33/Ensembl 110) | [gencode.vM33.primary_assembly.annotation.gtf] |
| [ENSEMBL (COVID-19)] | SARS-CoV-2 reference (ENSEMBL assembly: ASM985889v3) | [Sars\_cov\_2.ASM985889v3.101.gtf] |
| [ThermoFisher ERCC Spike-Ins] | ThermoFisher ERCC RNA Spike-In Control Mixes (Cat # 4456740, 4456739) | [cms_095047.txt] |

[RNA Spike-In Control Mixes]: https://www.thermofisher.com/document-connect/document-connect.html?url=https%3A%2F%2Fassets.thermofisher.com%2FTFS-Assets%2FLSG%2Fmanuals%2Fcms_086340.pdf&title=VXNlciBHdWlkZTogRVJDQyBSTkEgU3Bpa2UtSW4gQ29udHJvbCBNaXhlcyAoRW5nbGlzaCAp

[GENCODE (Human)]: https://www.gencodegenes.org/human/
[gencode.v44.primary_assembly.annotation.gtf]: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz

[GENCODE (Mouse)]: https://www.gencodegenes.org/mouse/
[gencode.vM33.primary_assembly.annotation.gtf]: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.primary_assembly.annotation.gtf.gz

[ENSEMBL (COVID-19)]: https://covid-19.ensembl.org/index.html
[Sars\_cov\_2.ASM985889v3.101.gtf]: https://ftp.ensemblgenomes.org/pub/viruses/gtf/sars_cov_2/Sars_cov_2.ASM985889v3.101.gtf.gz

[ThermoFisher ERCC Spike-Ins]: https://www.thermofisher.com/order/catalog/product/4456740#/4456740
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
            the most accurate child of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0002772"><code>"EFO:0002772"</code></a> for <i>assay by molecule</i>
          </li>
          <li>
            the most accurate child of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010183"><code>"EFO:0010183"</code></a>  for <i>single cell library construction</i>
          </li></ul>
        An assay based on 10X Genomics products SHOULD either be <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008995"><code>"EFO:0008995"</code></a> for <i>10x technology</i> or <b>preferably</b> its most accurate child. An assay based on <i>SMART (Switching Mechanism at the 5' end of the RNA Template) or SMARTer technology</i> SHOULD either be <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010184"><code>"EFO:0010184"</code></a> for <i>Smart-like</i> or preferably its most accurate child.<br><br>
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
              <td><i>10x 5' v1</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0011025"><code>"EFO:0011025"</code></a></td>
            </tr>
            <tr>
              <td><i>10x 5' v2</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0009900"><code>"EFO:0009900"</code></a></td>
            </tr>
            <tr>
              <td><i>Smart-seq2</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008931"><code>"EFO:0008931"</code></a></td>
            </tr>
            <tr>
              <td><i>Visium Spatial Gene Expression</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010961"><code>"EFO:0010961"</code></a></td>
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
        <td>categorical with <code>str</code> categories. This MUST be a CL term or <code>"unknown"</code> if no appropriate term can be found (e.g. the cell type is unknown). The following terms MUST NOT be used:
        <ul><li>
          <a href="https://www.ebi.ac.uk/ols4/ontologies/cl/terms?obo_id=CL:0000255"><code>"CL:0000255"</code></a> for <i>eukaryotic cell</i>
        </li>
        <li>
          <a href="https://www.ebi.ac.uk/ols4/ontologies/cl/terms?obo_id=CL:0000257"><code>"CL:0000257"</code></a> for <i>Eumycetozoan cell</i>
        </li>
        <li>
            <a href="https://www.ebi.ac.uk/ols4/ontologies/cl/terms?obo_id=CL:0000548"><code>"CL:0000548"</code></a> for <i>animal cell</i>
         </li></ul>
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
        <td>categorical with <code>str</code> categories. If unavailable, this MUST be <code>"unknown".</code> <br><br>If <code>organism_ontolology_term_id</code> is <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9606"><code>"NCBITaxon:9606"</code></a> for <i>Homo sapiens</i>, this MUST be the most accurate HsapDv term with the following STRONGLY RECOMMENDED:
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
              <td>Embryonic stage</td>
              <td>A term from the set of <a href="http://www.ontobee.org/search?ontology=HSAPDV&keywords=carnegie&submit=Search+terms">Carnegie stages 1-23</a><br>(up to 8 weeks after conception; e.g. <a href="https://www.ebi.ac.uk/ols4/ontologies/hsapdv/classes?obo_id=HsapDv%3A0000003">HsapDv:0000003</a>)</td>
            </tr>
            <tr>
              <td>Fetal development</td>
              <td>A term from the set of <a href="http://www.ontobee.org/search?ontology=HSAPDV&keywords=post-fertilization&submit=Search+terms">9 to 38 week post-fertilization human stages</a><br>(9 weeks after conception and before birth; e.g. <a href="https://www.ebi.ac.uk/ols4/ontologies/hsapdv/classes?obo_id=HsapDv%3A0000046">HsapDv:0000046</a>)</td>
            </tr>
            <tr>
              <td>After birth for the<br>first 12 months</td>
              <td>A term from the set of <a href="http://www.ontobee.org/search?ontology=HSAPDV&keywords=month-old&submit=Search+terms">1 to 12 month-old human stages</a><br>(e.g. <a href="https://www.ebi.ac.uk/ols4/ontologies/hsapdv/classes?obo_id=HsapDv%3A">HsapDv:0000174)</a></td>
            </tr>
            <tr>
              <td>After the first 12<br>months post-birth</td>
              <td>A term from the set of <a href="http://www.ontobee.org/search?ontology=HSAPDV&keywords=year-old&submit=Search+terms">year-old human stages</a><br>(e.g. <a href="https://www.ebi.ac.uk/ols4/ontologies/hsapdv/classes?obo_id=HsapDv%3A0000246">HsapDv:0000246)</a></td>
            </tr>
          </tbody></table>
          <br>If <code>organism_ontolology_term_id</code> is 
          <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A10090"><code>"NCBITaxon:10090"</code></a> for <i>Mus musculus</i>, this MUST be the most accurate MmusDv term with the following STRONGLY RECOMMENDED:
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
              <td>From the time of conception<br>to 1 month after birth</td>
              <td>A term from the set of <a href="http://www.ontobee.org/search?ontology=MMUSDV&keywords=theiler+stage&submit=Search+terms">Theiler stages</a><br>(e.g. <a href="https://www.ebi.ac.uk/ols4/ontologies/mmusdv/classes?obo_id=MmusDv%3A0000003">MmusDv:0000003</a>)</td>
            </tr>
            <tr>
              <td>From 2 months after birth</td>
              <td>A term from the set of <a href="http://www.ontobee.org/search?ontology=MMUSDV&keywords=month-old&submit=Search+terms"> month-old stages</a><br>(e.g. <a href="https://www.ebi.ac.uk/ols4/ontologies/mmusdv/classes?obo_id=MmusDv%3A0000062">MmusDv:0000062)</a></td>
            </tr>
          </tbody></table>
          <br> Otherwise, for all other organisms this MUST be the most accurate child of <a href="https://www.ebi.ac.uk/ols4/ontologies/uberon/classes?obo_id=UBERON%3A0000105"<code>UBERON:0000105</code></a> for <i>life cycle stage</i>, excluding <a href="https://www.ebi.ac.uk/ols4/ontologies/uberon/classes?obo_id=UBERON%3A0000071"<code>UBERON:0000071</code></a> for <i>death stage</i>.
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
          <li>the most accurate child of <a href="https://www.ebi.ac.uk/ols4/ontologies/mondo/classes?obo_id=MONDO%3A0000001"><code>"MONDO:0000001"</code></a> for <i>disease</i></li>
          <li><a href="https://www.ebi.ac.uk/ols4/ontologies/mondo/classes?obo_id=MONDO%3A0021178"><code>"MONDO:0021178"</code></a> for <i>injury</i> or <b>preferably</b> its most accurate child</li>
       </ul>
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
        <td>categorical with <code>str</code> categories. This MUST be free-text that identifies a unique individual that data were derived from. It is STRONGLY RECOMMENDED that this identifier be designed so that it is unique to:<br><br>
          <ul><li>a given individual within the collection of datasets that includes this dataset</li>
          <li>a given individual across all collections in CELLxGENE Discover</li></ul><br>
          It is STRONGLY RECOMMENDED that <code>"pooled"</code> be used  for observations from a sample of multiple individuals that were not confidently assigned to a single individual through demultiplexing.<br><br>It is STRONGLY RECOMMENDED that <code>"unknown"</code> ONLY be used for observations in a dataset when it is not known which observations are from the same individual.<br><br>
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
        <td><code>bool</code>. This MUST be <code>True</code> if this is the canonical instance of this cellular observation and <code>False</code> if not. This is commonly <code>False</code> for meta-analyses reusing data or for secondary views of data.
        </td>
    </tr>
</tbody></table>
<br>

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
        <td>categorical with <code>str</code> categories. This MUST be a child of <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A33208"<code>NCBITaxon:33208</code></a> for <i>Metazoa</i>.
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
        categorical with <code>str</code> categories. If
        <code>organism_ontolology_term_id</code> is
        <code>"NCBITaxon:9606"</code> for <i>Homo sapiens</i>,
        the value MUST be formatted as one or more comma-separated (with no leading or trailing spaces) HANCESTRO
        terms in ascending lexical order or <code>"unknown"</code> if unavailable.<br><br>For example, if the terms are <code>"HANCESTRO:0014</code> and <code>HANCESTRO:0005"</code> then the value of <code>self_reported_ethnicity_ontology_term_id</code> MUST be <code>"HANCESTRO:0005,HANCESTRO:0014"</code>.<br><br>The following terms MUST NOT be used:<br /><br />
        <ul>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0002?lang=en"
              ><code>"HANCESTRO:0002"</code></a
            >
            for <i>regions</i> and its children
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0003?lang=en"
              ><code>"HANCESTRO:0003"</code></a
            >
            for <i>country</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0004?lang=en"
              ><code>"HANCESTRO:0004"</code></a
            >
            for <i>ancestry category</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0018?lang=en"
              ><code>"HANCESTRO:0018"</code></a
            >
            for <i>uncategorised population</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0290?lang=en"
              ><code>"HANCESTRO:0290"</code></a
            >
            for <i>genetically isolated population</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0304?lang=en"
              ><code>"HANCESTRO:0304"</code></a
            >
            for <i>ancestry status</i> and its children
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0323?lang=en"
              ><code>"HANCESTRO:0323"</code></a
            >
            for <i>Finnish founder</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0324?lang=en"
              ><code>"HANCESTRO:0324"</code></a
            >
            for <i>Dutch founder</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0551?lang=en"
              ><code>"HANCESTRO:0551"</code></a
            >
            for <i>genetically homogenous Irish</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0554?lang=en"
              ><code>"HANCESTRO:0554"</code></a
            >
            for <i>Silk Road founder</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0555?lang=en"
              ><code>"HANCESTRO:0555"</code></a
            >
            for <i>Arab Israeli founder</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0557?lang=en"
              ><code>"HANCESTRO:0557"</code></a
            >
            for <i>Costa Rican founder</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0558?lang=en"
              ><code>"HANCESTRO:0558"</code></a
            >
            for <i>French Canadian founder</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0559?lang=en"
              ><code>"HANCESTRO:0559"</code></a
            >
            for <i>Italian founder</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0560?lang=en"
              ><code>"HANCESTRO:0560"</code></a
            >
            for <i>Northern Finnish founder</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0561?lang=en"
              ><code>"HANCESTRO:0561"</code></a
            >
            for <i>Romanian founder</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0564?lang=en"
              ><code>"HANCESTRO:0564"</code></a
            >
            for <i>Vis founder</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0565?lang=en"
              ><code>"HANCESTRO:0565"</code></a
            >
            for <i>Split founder</i>
          </li>
          <li>
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0566?lang=en"
              ><code>"HANCESTRO:0566"</code></a
            >
            for <i>undefined ancestry population</i>
          </li>
          <li>
            The imported GEO term
            <a
              href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FGEO_000000374?lang=en"
              ><code>"GEO:000000374"</code></a
            >
            for <i>continent</i> and its children:
            <ul>
              <li>
                <a
                  href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0029?lang=en"
                  ><code>"HANCESTRO:0029"</code></a
                >
                for <i>Africa</i>
              </li>
              <li>
                <a
                  href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0030?lang=en"
                  ><code>"HANCESTRO:0030"</code></a
                >
                for <i>Asia</i>
              </li>
              <li>
                <a
                  href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0031?lang=en"
                  ><code>"HANCESTRO:0031"</code></a
                >
                for <i>Europe</i>
              </li>
              <li>
                <a
                  href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0032?lang=en"
                  ><code>"HANCESTRO:0032"</code></a
                >
                for <i>Oceania</i>
              </li>
              <li>
                <a
                  href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0033?lang=en"
                  ><code>"HANCESTRO:0033"</code></a
                >
                for <i>Latin America and the Caribbean</i>
              </li>
              <li>
                <a
                  href="https://www.ebi.ac.uk/ols4/ontologies/hancestro/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHANCESTRO_0034?lang=en"
                  ><code>"HANCESTRO:0034"</code></a
                >
                for <i>Northern America</i>
              </li>
            </ul>
          </li>
        </ul>
        <br />Otherwise, for all other organisms the
        <code>str</code> value MUST be <code>"na"</code>.
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
        <td>categorical with <code>str</code> categories. This MUST be a child of <a href="https://www.ebi.ac.uk/ols4/ontologies/pato/classes?obo_id=PATO%3A0001894">PATO:0001894</a> for  <i>phenotypic sex</i> or <code>"unknown"</code> if unavailable.
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
              <td><i>10x transcription profiling</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030080"><code>EFO:0030080</code></a>] and its children</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr> 
            <tr>
              <td><i>ATAC-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0007045"><code>EFO:0007045</code></a>] and its children</td>
              <td><code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>BD Rhapsody Whole Transcriptome Analysis</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0700003"><code>EFO:0700003</code></a>]</td>
              <td><code>"cell"</code></td>
           </tr>
            <tr>
              <td><i>BD Rhapsody Targeted mRNA</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0700004"><code>EFO:0700004</code></a>]</td>
              <td><code>"cell"</code></td>
           </tr>
            <tr>
              <td><i>CEL-seq2</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010010"><code>EFO:0010010</code></a>]</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>CITE-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0009294"><code>EFO:0009294</code></a>] and its children</td>
              <td><code>"cell"</code></td>
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
              <td><i>microwell-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030002"><code>EFO:0030002</code></a>]</td>
              <td><code>"cell"</code></td>
           </tr>    
            <tr>
              <td><i>Patch-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008853"><code>EFO:0008853</code></a>]</td>
              <td><code>"cell"</code></td>
           </tr>
            <tr>
              <td><i>sci-Plex</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030026"><code>EFO:0030026</code></a>]</td>
              <td><code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>sci-RNA-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010550"><code>EFO:0010550</code></a>]</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>Seq-Well</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008919"><code>EFO:0008919</code></a>] and its children</td>
              <td><code>"cell"</code></td>
           </tr>
            <tr>
              <td><i>Smart-like</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010184"><code>EFO:0010184</code></a>] and its children</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>smFISH</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0009918"><code>EFO:0009918</code></a>] and its children</td>
              <td><code>"na"</code></td>
           </tr>   
            <tr>
              <td><i>snmC-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008939"><code>EFO:0008939</code></a>]</td>
              <td><code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>snmC-seq2</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0030027"><code>EFO:0030027</code></a>]</td>
              <td><code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>spatial proteomics</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0700000"><code>EFO:0700000</code></a>] and its children</td>
              <td><code>"na"</code></td>
           </tr>
            <tr>
              <td><i>spatial transcriptomics</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008994"><code>EFO:0008994</code></a>] and its children</td>
              <td><code>"na"</code></td>
           </tr> 
            <tr>
              <td><i>SPLiT-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0009919"><code>EFO:0009919</code></a>]</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr> 
            <tr>
              <td><i>STRT-seq</i> [<a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008953"><code>EFO:0008953</code></a>]</td>
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
        <td>categorical with <code>str</code> categories. This MUST be <code>"tissue"</code>, <code>"organoid"</code>, or <code>"cell culture"</code>.
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
        <td>categorical with <code>str</code> categories. If <code>tissue_type</code> is <code>"tissue"</code> or <code>"organoid"</code>, this MUST be the most accurate child of <a href="https://www.ebi.ac.uk/ols4/ontologies/uberon/classes?obo_id=UBERON%3A0001062"><code>UBERON:0001062</code></a> for <i>anatomical entity</i>.<br><br> If <code>tissue_type</code> is <code>"cell culture"</code> this MUST follow the requirements for <code>cell_type_ontology_term_id<code>.</td>
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
        <td>categorical with <code>str</code> categories. This MUST be <code>"unknown"</code> if the value of <code>cell_type_ontology_term_id</code> is <code>"unknown"</code>; otherwise, this MUST be the human-readable name assigned to the value of <code>cell_type_ontology_term_id</code>.
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
        <td>categorical with <code>str</code> categories. This MUST be <code>"unknown"</code> if the value of <code>development_stage_ontology_term_id</code> is <code>"unknown"</code>; otherwise, this MUST be the human-readable name assigned to the value of <code>development_stage_ontology_term_id</code>.
        </td>
    </tr>
</tbody></table>
<br>

### disease

<table><tbody>
    <tr>
      <th>Key</th>
      <td>disease</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be the human-readable name assigned to the value of <code>disease_ontology_term_id</code>.
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
        <td>categorical with <code>str</code> categories. This MUST be <code>"na"</code> if the value of <code>self_reported_ethnicity_ontology_term_id</code> is <code>"na"</code>. This MUST be <code>"unknown"</code> if the value of <code>self_reported_ethnicity_ontology_term_id</code> is <code>"unknown"</code>. Otherwise, this MUST be one or more comma-separated (with no leading or trailing spaces) human-readable names for the terms in <code>self_reported_ethnicity_ontology_term_id</code> in the same order.<br><br> For example, if the value of <code>self_reported_ethnicity_ontology_term_id</code> is <code>"HANCESTRO:0005,HANCESTRO:0014"</code> then the value of <code>self_reported_ethnicity</code> is <code>"European,Hispanic or Latin American"</code>.
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
        <td>categorical with <code>str</code> categories. This MUST be <code>"unknown"</code> if the value of  <code>sex_ontology_term_id</code> is <code>"unknown"</code>; otherwise, this MUST be the human-readable name assigned to the value of <code>sex_ontology_term_id</code>.
        </td>
    </tr>
</tbody></table>
<br>

### tissue

<table><tbody>
    <tr>
      <th>Key</th>
      <td>tissue</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>CELLxGENE Discover MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be the human-readable name assigned to the value of <code>tissue_ontology_term_id</code>.
        </td>
    </tr>
</tbody></table>
<br>

## `obsm` (Embeddings)


The size of the ndarray stored for a key in `obsm` MUST NOT be zero.


To display a dataset in CELLxGENE Explorer, Curators MUST annotate **one or more** embeddings of at least two-dimensions (e.g. tSNE, UMAP, PCA, spatial coordinates) as `numpy.ndarrays` in `obsm`.<br><br>

### X_{suffix}

<table><tbody>
    <tr>
      <th>Key</th>
      <td>X_{suffix} with the following requirements:<br><br>
      <ul>
        <li>{suffix} MUST be at least one character in length.</li>
        <li>The first character of {suffix} MUST be a letter of the alphabet and the remaining characters MUST be alphanumeric characters, <code>'_'</code>, <code>'-'</code>, or <code>'.'</code> (This is equivalent to the regular expression pattern <code>"^[a-zA-Z][a-zA-Z0-9_.-]*$"</code>.)</li>
      </ul><br>
      {suffix} is presented as text to users in the <b>Embedding Choice</b> selector in CELLxGENE Explorer so it is STRONGLY RECOMMENDED that it be descriptive.<br><br>See also <code>default_embedding</code> in <code>uns</code>.</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
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
        <td><code>str</code>. If the feature is a gene then this MUST be an ENSEMBL term. If the feature is a RNA Spike-In Control Mix then this MUST be an ERCC Spike-In identifier (e.g. <code>"ERCC-0003"</code>).<br><br> The index of the <code>pandas.DataFrame</code> MUST contain unique identifiers for features. If present, the index of <code>raw.var</code> MUST be identical to the index of <code>var</code>.<br><br></td>
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
        <td><code>bool</code>. This MUST be <code>True</code> if the feature was filtered out in the normalized matrix (<code>X</code>) but is present in the raw matrix (<code>raw.X</code>). The value for all cells of the given feature in the normalized matrix MUST be <code>0</code>.  <br><br>Otherwise, this MUST be <code>False</code>. </td>
    </tr>
</tbody></table>
<br>

Curators MUST NOT annotate the following columns in the `var` dataframe and if present, the `raw.var` dataframe.

When a dataset is uploaded, CELLxGENE Discover MUST automatically add the matching human-readable name for the corresponding feature biotype, identifier, and the NCBITaxon term for the reference organism to the `var` and `raw.var` dataframes. In addition, it MUST
add the feature length. 

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
        <code>uint</code> number of base-pairs (bps). If the <code>feature_biotype</code> is <code>"gene"</code>, then the value is calculated by creating non-overlapping concatenated exons across all isoforms of the gene, and then adding up their length in base-pairs. This approach is modeled on the "length of merged exons of isoforms of a gene" from <a href="https://doi.org/10.1093/bioinformatics/btac561">GTFtools: a software package for analyzing various features of gene models.</a><br><br> If <code>feature_biotype</code> is NOT <code>"gene"</code> then the value MUST be set to 0.
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
        <td><code>str</code>. If the <code>feature_biotype</code> is <code>"gene"</code> then this MUST be the human-readable ENSEMBL gene name assigned to the feature identifier in <code>var.index</code>. If the <code>feature_biotype</code> is <code>"spike-in"</code> then this MUST be the ERCC Spike-In identifier appended with <code>" (spike-in control)"</code>.
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
              <td><i>Homo sapiens</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9606"><code>"NCBITaxon:9606"</code></a></td>
            </tr>
            <tr>
              <td><i>Mus musculus</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A10090"><code>"NCBITaxon:10090"</code></a></td>
            </tr>
            <tr>
              <td><i>SARS-CoV-2</i></td>
              <td><a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A2697049"><code>"NCBITaxon:2697049"</code></a></td>
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

## `varm`

The size of the ndarray stored for a key in `varm` MUST NOT be zero.
<br>

## `varp`
The size of the ndarray stored for a key in `varp` MUST NOT be zero.
<br>

## `uns` (Dataset Metadata)

`uns` is a ordered dictionary with a `str` key. The data stored as a value for a key in `uns` MUST be `True`, `False`, `None`, or its size MUST NOT be zero.

Curators MUST annotate the following keys and values in `uns`:

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
        <li>organism</li>
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
        <code>numpy.ndarray</code>. This MUST be a 1-D array of shape <code>(, c)</code>, where <code>c</code> is greater than or equal to the<br> number of unique categories in the {column} as calculated by:<br><br>
           <samp>len(anndata.obs.{column}.unique())</samp><br><br>
        The color code at the Nth position in the <code>ndarray</code> corresponds to the Nth category of <samp>anndata.obs.{column}.unique()</samp>.<br><br>For example, if <code>cell_type_ontology_term_id</code> includes two unique categories:<br><br>
        <samp>anndata.obs.cell_type_ontology_term_id.unique()</samp><br><br>
        <samp>['CL:0000057', 'CL:0000115']<br>Categories (2, object): ['CL:0000057', 'CL:0000115']</samp><br><br>then <code>cell-type_ontology_term_id_colors</code> MUST contain two or more colors such as:<br><br>
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
          This MUST be <code>"https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.0.0/schema.md"</code>.
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
          This MUST be <code>"5.0.0"</code>.
        </td>
    </tr>
</tbody></table>

## Appendix A. Changelog

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
  * Updated the requirements for `assay_ontology_term_id` to not allow  the parent terms `EFO:0002772` for _assay by molecule_ and `EFO:0010183` for _single cell library construction_. Their most accurate children are still valid. 
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

* The canonical data format was updated from AnnData 0.7 to 0.8.
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
