
# Schema

Contact: brianraymor@chanzuckerberg.com

Document Status: _Draft_

Version: 3.0.0

The key words "MUST", "MUST NOT", "REQUIRED", "SHALL", "SHALL NOT", "SHOULD", "SHOULD NOT", "RECOMMENDED", "NOT RECOMMENDED" "MAY", and "OPTIONAL" in this document are to be interpreted as described in [BCP 14](https://tools.ietf.org/html/bcp14), [RFC2119](https://www.rfc-editor.org/rfc/rfc2119.txt), and [RFC8174](https://www.rfc-editor.org/rfc/rfc8174.txt) when, and only when, they appear in all capitals, as shown here.

## Background

cellxgene aims to support the publication, sharing, and exploration of single-cell datasets. Building on those published datasets, cellxgene seeks to create references of the phenotypes and composition of cells that make up human tissues.

Creating references from multiple datasets requires some harmonization of metadata and features in the cellxgene Data Portal. But if that harmonization is too onerous, it will burden the goal of rapid data sharing. cellxgene balances publishing and reference creation needs by requiring datasets hosted in the cellxgene Data Portal to include a small set of metadata readily available from data submitters.

This document describes the schema, a type of contract, that cellxgene requires all datasets to adhere to so that it can enable searching, filtering, and integration of datasets it hosts.

Note that the requirements in the schema are just the minimum required information. Datasets often have additional metadata, which is preserved in datasets submitted to the cellxgene Data Portal.

## Overview

This schema supports multiple assay types. Each assay takes the form of one or more two-dimensional matrices whose values are quantitative measures of the phenotypes of cells.

The schema additionally describes how the dataset, genes, and cells are annotated to describe the biological and technical characteristics of the data.

This document is organized by:

* [General requirements](#general-requirements)
* [`X` (Matrix layers)](#x-matrix-layers), which describe the data required for different assays
* [`obs` (Cell metadata)](#obs-cell-metadata), which describe each cell in the dataset
* [`var` and `raw.var` (Gene metadata)](#var-and-rawvar-gene-metadata), which describe each gene in the dataset
* [`obsm` (Embeddings)](#obsm-embeddings), which describe each embedding in the dataset
* [`uns` (Dataset metadata)](#uns-dataset-metadata), which describe the dataset as a whole

## General Requirements

**AnnData.** The canonical data format for the cellxgene Data Portal is HDF5-backed [AnnData](https://anndata.readthedocs.io/en/latest) as written by version 0.8 of the anndata library.  Part of the rationale for selecting this format is to allow cellxgene to access both the data and metadata within a single file. The schema requirements and definitions for the AnnData `X`, `obs`, `var`, `raw.var`, `obsm`, and `uns` attributes are described below.

All data submitted to the cellxgene Data Portal is automatically converted to a Seurat V3 object that can be loaded by the R package Seurat. See the [Seurat encoding](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/2.0.0/seurat_encoding.md) for further information.

**Organisms**. Data MUST be from a Metazoan organism or SARS-COV-2 and defined in the NCBI organismal classification. For data that is neither Human, Mouse, nor SARS-COV-2, features MUST be translated into orthologous genes from the pinned Human and Mouse gene annotations.

**Reserved Names**. The names of metadata keys MUST NOT start with `"__"`. The names of the metadata keys specified by the schema are reserved and MUST be unique. For example, duplicate <code>"feature_biotype"</code> keys in AnnData <code>var</code> are not allowed. 

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

**No Personal Identifiable Information (PII)**.  This is not strictly enforced by validation because it is difficult for software to predict what is and is not PII; however, curators MUST agree to the data submission policies of the cellxgene Data Portal on behalf of data submitters which includes this requirement:

> It is my responsibility to ensure that this data is not identifiable. In particular, I commit that I will remove any [direct personal identifiers](https://docs.google.com/document/d/1sboOmbafvMh3VYjK1-3MAUt0I13UUJfkQseq8ANLPl8/edit) in the metadata portions of the data, and that CZI may further contact me if it believes more work is needed to de-identify it.

This includes names, emails, or other PII for researchers or curators involved in the data generation and submission.

#### *Note on types*
The types below are python3 types. Note that a python3 `str` is a sequence of Unicode code points, which is stored null-terminated and UTF-8-encoded by anndata.

## `X` (Matrix Layers)

The data stored in the `X` data matrix is the data that is viewable in cellxgene Explorer. cellxgene does not impose any additional constraints on the `X` data matrix.

In any layer, if a matrix has 50% or more values that are zeros, it is STRONGLY RECOMMENDED that the matrix be encoded as a [`scipy.sparse.csr_matrix`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html).

cellxgene's matrix layer requirements are tailored to optimize data reuse. Because each assay has different characteristics, the requirements differ by assay type. In general, cellxgene requires submission of "raw" data suitable for computational reuse when a standard raw matrix format exists for an assay. It is STRONGLY RECOMMENDED to also include a "normalized" matrix with processed values ready for data analysis and suitable for visualization in cellxgene Explorer. So that cellxgene's data can be provided in download formats suitable for both R and Python, the schema imposes the following requirements:

*   All matrix layers MUST have the same shape, and have the same cell labels and gene labels.
*   Because it is impractical to retain all barcodes in raw and normalized matrices, any cell filtering MUST be applied to both.
    By contrast, those wishing to reuse datasets require access to raw gene expression values, so genes SHOULD NOT be filtered from either dataset.
    Summarizing, any cell barcodes that are removed from the data MUST be filtered from both raw and normalized matrices and genes SHOULD NOT be filtered from the raw matrix.
*   Any genes that publishers wish to filter from the normalized matrix MAY have their values replaced by zeros and MUST be flagged in the column [`feature_is_filtered`](#feature_is_filtered) of [`var`](#var-and-rawvar-gene-metadata), which will mask them from exploration.
*   Additional layers provided at author discretion MAY be stored using author-selected keys, but MUST have the same cells and genes as other layers. It is STRONGLY RECOMMENDED that these layers have names that accurately summarize what the numbers in the layer represent (e.g. `"counts_per_million"`, `"SCTransform_normalized"`, or `"RNA_velocity_unspliced"`).

The following table describes the matrix data and layers requirements that are **assay-specific**. If an entry in the table is empty, the schema does not have any other requirements on data in those layers beyond the ones listed above.

| Assay | "raw" required? | "raw" location | "normalized" required? | "normalized" location |
|-|-|-|-|-|
| scRNA-seq (UMI, e.g. 10x v3) | REQUIRED. Values MUST be de-duplicated molecule counts. | `AnnData.raw.X` unless no "normalized" is provided, then `AnnData.X` | STRONGLY RECOMMENDED | `AnnData.X` |
| scRNA-seq (non-UMI, e.g. SS2) | REQUIRED. Values MUST be one of read counts (e.g. FeatureCounts) or  estimated fragments (e.g. output of RSEM). | `AnnData.raw.X` unless no "normalized" is provided, then `AnnData.X` | STRONGLY RECOMMENDED | `AnnData.X` |
| Accessibility (e.g. ATAC-seq, mC-seq) | NOT REQUIRED | | REQUIRED | `AnnData.X` | STRONGLY RECOMMENDED |
|||||

## Integration Metadata

cellxgene requires ontology terms to enable search, comparison, and integration of data.
Ontology terms for cell metadata MUST use [OBO-format identifiers](http://www.obofoundry.org/id-policy.html), meaning a CURIE (prefixed identifier) of the form **Ontology:Identifier**.
For example, [EFO:0000001](https://www.ebi.ac.uk/ols/ontologies/efo/terms?short_form=EFO_0000001) is a term in the Experimental Factor Ontology (EFO).

The most accurate ontology term MUST always be used. If an exact or approximate ontology term is not available, a new term may be requested:

- For the [Cell Ontology], data submitters may [suggest a new term](https://github.com/obophenotype/cell-ontology/issues/new?assignees=bvarner-ebi&labels=new+term+request%2C+cellxgene&template=a_adding_term_cellxgene.md&title=%5BNTR-cxg%5D) and [notify the curation team](mailto:cellxgene@chanzuckerberg.com) of the pending term request, so that the datasets can be updated once the term is available.

  To meet cellxgene schema requirements, the most accurate available CL term MUST be used until the new term is available. For example if `cell_type_ontology_term_id` describes a relay interneuron, but the most accurate available term in the CL ontology is [CL:0000099](https://www.ebi.ac.uk/ols/ontologies/cl/terms?obo_id=CL:0000099) for *interneuron*, then the interneuron term can be used to fulfill this requirement and ensures that users searching for "neuron" are able to find these data. If no appropriate high-level term can be found or the cell type is unknown, then the most accurate term is [CL:0000003](https://www.ebi.ac.uk/ols/ontologies/cl/terms?obo_id=CL:0000003) for *native cell*. Users will still be able to access more specific cell type annotations that have been submitted with the dataset (but aren't required by the schema).

   
- For all other ontologies, data submitters may submit a [request to the curation team](mailto:cellxgene@chanzuckerberg.com) during the submission process.

Terms documented as obsolete in an ontology MUST NOT be used. For example, [EFO:0009310](http://www.ebi.ac.uk/efo/EFO_0009310) for *obsolete_10x v2* was marked as obsolete in EFO version 3.31.0 and replaced by [EFO:0009899](http://www.ebi.ac.uk/efo/EFO_0009899) for *10x 3' v2*.

### Required Ontologies

The following ontology dependencies are *pinned* for this version of the schema.

| Ontology | OBO Prefix | Required version |
|:--|:--|:--|
| [Cell Ontology] | CL | [cl.owl] : [2021-08-10]|
| [Experimental Factor Ontology] | EFO | [efo.owl] : [2021-08-16 EFO 3.33.0]
| [Human Ancestry Ontology] | HANCESTRO |[hancestro.owl] : [2021-01-04 (2.5)] |
| [Human Developmental Stages] |  HsapDv | [hsapdv.owl] : 2020-03-10 |
| [Mondo Disease Ontology] | MONDO |[mondo.owl] : [2021-08-11] |
| [Mouse Developmental Stages]| MmusDv |  [mmusdv.owl] : 2020-03-10 |
| [NCBI organismal classification] |  NCBITaxon | [ncbitaxon.owl] : [2021-06-10] |
| [Phenotype And Trait Ontology] | PATO | [pato.owl] : [2021-08-06] |  |
| [Uberon multi-species anatomy ontology] |  UBERON | [uberon.owl] : [2021-07-27] |
| | | |

[Cell Ontology]: http://obofoundry.org/ontology/cl.html
[2021-08-10]: https://github.com/obophenotype/cell-ontology/releases/tag/v2021-08-10
[cl.owl]: https://github.com/obophenotype/cell-ontology/blob/v2021-08-10/cl.owl

[Experimental Factor Ontology]: http://www.ebi.ac.uk/efo
[2021-08-16 EFO 3.33.0]: https://github.com/EBISPOT/efo/releases/tag/v3.33.0
[efo.owl]: https://github.com/EBISPOT/efo/releases/download/v3.33.0/efo.owl

[Human Ancestry Ontology]: http://www.obofoundry.org/ontology/hancestro.html
[2021-01-04 (2.5)]: https://github.com/EBISPOT/ancestro/releases/tag/2.5
[hancestro.owl]: https://github.com/EBISPOT/ancestro/blob/2.5/hancestro.owl

[Human Developmental Stages]: http://obofoundry.org/ontology/hsapdv.html
[hsapdv.owl]: http://purl.obolibrary.org/obo/hsapdv.owl

[Mondo Disease Ontology]: http://obofoundry.org/ontology/mondo.html
[2021-08-11]: https://github.com/monarch-initiative/mondo/releases/tag/v2021-08-11
[mondo.owl]: https://github.com/monarch-initiative/mondo/releases/download/v2021-08-11/mondo.owl

[Mouse Developmental Stages]: http://obofoundry.org/ontology/mmusdv.html
[mmusdv.owl]: http://purl.obolibrary.org/obo/mmusdv.owl

[NCBI organismal classification]: http://obofoundry.org/ontology/ncbitaxon.html
[2021-06-10]: https://github.com/obophenotype/ncbitaxon/releases/tag/v2021-06-10
[ncbitaxon.owl]: https://github.com/obophenotype/ncbitaxon/releases/download/v2021-06-10/ncbitaxon.owl.gz

[Phenotype And Trait Ontology]: http://www.obofoundry.org/ontology/pato.html
[2021-08-06]: https://github.com/pato-ontology/pato/releases/tag/v2021-08-06
[pato.owl]: https://github.com/pato-ontology/pato/blob/v2021-08-06/pato.owl

[Uberon multi-species anatomy ontology]: http://www.obofoundry.org/ontology/uberon.html
[2021-07-27]: https://github.com/obophenotype/uberon/releases/tag/v2021-07-27
[uberon.owl]: https://github.com/obophenotype/uberon/blob/v2021-07-27/uberon.owl

### Required Gene Annotations

cellxgene requires ENSEMBL identifiers for genes and [External RNA Controls Consortium (ERCC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4978944/) identifiers for [RNA Spike-In Control Mixes] to ensure that all datasets it stores measure the same features and can therefore be integrated.

The following gene annotation dependencies are *pinned* for this version of the schema. For multi-organism experiments, cells from any Metazoan organism are allowed as long as orthologs from the following organism annotations are used.

| Source | Required version | Download |
|:--|:--|:--|
| [GENCODE (Human)] | Human reference GRCh38 (GENCODE v38/Ensembl 104) | [gencode.v38.primary_assembly.annotation.gtf] |
| [GENCODE (Mouse)] | Mouse reference GRCm39 (GENCODE vM27/Ensembl 104) | [gencode.vM27.primary_assembly.annotation.gtf] |
| [ENSEMBL (COVID-19)] | SARS-CoV-2 reference (ENSEMBL assembly: ASM985889v3) | [Sars\_cov\_2.ASM985889v3.101.gtf] |
| [ThermoFisher ERCC Spike-Ins] | ThermoFisher ERCC RNA Spike-In Control Mixes (Cat # 4456740, 4456739) | [cms_095047.txt] |

[RNA Spike-In Control Mixes]: https://www.thermofisher.com/document-connect/document-connect.html?url=https%3A%2F%2Fassets.thermofisher.com%2FTFS-Assets%2FLSG%2Fmanuals%2Fcms_086340.pdf&title=VXNlciBHdWlkZTogRVJDQyBSTkEgU3Bpa2UtSW4gQ29udHJvbCBNaXhlcyAoRW5nbGlzaCAp

[GENCODE (Human)]: https://www.gencodegenes.org/human/
[gencode.v38.primary_assembly.annotation.gtf]: http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz

[GENCODE (Mouse)]: https://www.gencodegenes.org/mouse/
[gencode.vM27.primary_assembly.annotation.gtf]: http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.primary_assembly.annotation.gtf.gz

[cellranger 2020-A (July 7, 2020) release]: https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build

[ENSEMBL (COVID-19)]: https://covid-19.ensembl.org/index.html
[Sars\_cov\_2.ASM985889v3.101.gtf]: https://ftp.ensemblgenomes.org/pub/viruses/gtf/sars_cov_2/Sars_cov_2.ASM985889v3.101.gtf.gz

[ThermoFisher ERCC Spike-Ins]: https://www.thermofisher.com/order/catalog/product/4456740#/4456740
[cms_095047.txt]: https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095047.txt


## `obs` (Cell Metadata)

`obs` is a [`pandas.DataFrame`](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html).

Curators MUST annotate the following columns in the `obs` dataframe:

### assay_ontology_term_id

<table><tbody>
    <tr>
      <th>Key</th>
      <td>assay_ontology_term_id</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be an EFO term and either:<br><br>
          <ul><li>
            <a href="http://www.ebi.ac.uk/efo/EFO_0002772"><code>"EFO:0002772"</code></a> for <i>assay by molecule</i> or <b>preferably</b> its most accurate child
          </li>
          <li>
            <a href="http://www.ebi.ac.uk/efo/EFO_0010183"><code>"EFO:0010183"</code></a>  for <i>single cell library construction</i> or <b>preferably</b> its most accurate child
          </li></ul>
        An assay based on 10X Genomics products SHOULD either be <a href="http://www.ebi.ac.uk/efo/EFO_0008995"><code>"EFO:0008995"</code></a> for <i>10x technology</i> or <b>preferably</b> its most accurate child. An assay based on <i>SMART (Switching Mechanism at the 5' end of the RNA Template) or SMARTer technology</i> SHOULD either be <a href="http://www.ebi.ac.uk/efo/EFO_0010184"><code>"EFO:0010184"</code></a> for <i>Smart-like</i> or preferably its most accurate child.<br><br>
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
              <td><a href="http://www.ebi.ac.uk/efo/EFO_0009899"><code>"EFO:0009899"</code></a></td>
            </tr>
            <tr>
              <td><i>10x 3' v3</i></td>
              <td><a href="http://www.ebi.ac.uk/efo/EFO_0009922"><code>"EFO:0009922"</code></a></td>
            </tr>
            <tr>
              <td><i>10x 5' v1</i></td>
              <td><a href="http://www.ebi.ac.uk/efo/EFO_0011025"><code>"EFO:0011025"</code></a></td>
            </tr>
            <tr>
              <td><i>10x 5' v2</i></td>
              <td><a href="http://www.ebi.ac.uk/efo/EFO_0009900"><code>"EFO:0009900"</code></a></td>
            </tr>
            <tr>
              <td><i>Smart-seq</i></td>
              <td><a href="http://www.ebi.ac.uk/efo/EFO_0008930"><code>"EFO:0008930"</code></a></td>
            </tr>
            <tr>
              <td><i>Smart-seq2</i></td>
              <td><a href="http://www.ebi.ac.uk/efo/EFO_0008931"><code>"EFO:0008931"</code></a></td>
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
      <td>Curator</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be a CL term.
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
      <td>Curator</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. If unavailable, this MUST be <code>"unknown".</code> <br><br>If <code>organism_ontolology_term_id</code> is <code>"NCBITaxon:9606"</code> for <i>Homo sapiens</i>, this MUST be the most<br>accurate HsapDv term with the following STRONGLY RECOMMENDED:
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
              <td>A term from the set of <a href="http://www.ontobee.org/search?ontology=HSAPDV&keywords=carnegie&submit=Search+terms">Carnegie stages 1-23</a><br>(up to 8 weeks after conception; e.g. <a href="http://purl.obolibrary.org/obo/HsapDv_0000003">HsapDv:0000003</a>)</td>
            </tr>
            <tr>
              <td>Fetal development</td>
              <td>A term from the set of <a href="http://www.ontobee.org/search?ontology=HSAPDV&keywords=post-fertilization&submit=Search+terms">9 to 38 week post-fertilization human stages</a><br>(9 weeks after conception and before birth; e.g. <a href="http://purl.obolibrary.org/obo/HsapDv_0000046">HsapDv:0000046</a>)</td>
            </tr>
            <tr>
              <td>After birth for the<br>first 12 months</td>
              <td>A term from the set of <a href="http://www.ontobee.org/search?ontology=HSAPDV&keywords=month-old&submit=Search+terms">1 to 12 month-old human stages</a><br>(e.g. <a href="http://purl.obolibrary.org/obo/HsapDv_0000174">HsapDv:0000174)</a></td>
            </tr>
            <tr>
              <td>After the first 12<br>months post-birth</td>
              <td>A term from the set of <a href="http://www.ontobee.org/search?ontology=HSAPDV&keywords=year-old&submit=Search+terms">year-old human stages</a><br>(e.g. <a href="http://purl.obolibrary.org/obo/HsapDv_0000246">HsapDv:0000246)</a></td>
            </tr>
          </tbody></table>
          <br>If <code>organism_ontolology_term_id</code> is <code>"NCBITaxon:10090"</code> for <i>Mus musculus</i></code>, this MUST be the most<br>accurate MmusDv term with the following STRONGLY RECOMMENDED:
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
              <td>A term from the set of <a href="http://www.ontobee.org/search?ontology=MMUSDV&keywords=theiler+stage&submit=Search+terms">Theiler stages</a><br>(e.g. <a href="http://purl.obolibrary.org/obo/MmusDv_0000003">MmusDv:0000003</a>)</td>
            </tr>
            <tr>
              <td>From 2 months after birth</td>
              <td>A term from the set of <a href="http://www.ontobee.org/search?ontology=MMUSDV&keywords=month-old&submit=Search+terms"> month-old stages</a><br>(e.g. <a href="http://purl.obolibrary.org/obo/MmusDv_0000062">MmusDv:0000062)</a></td>
            </tr>
          </tbody></table>
          <br> Otherwise, for all other organisms this MUST be the most accurate child of <a href="http://purl.obolibrary.org/obo/UBERON_0000105"<code>UBERON:0000105</code></a> for <i>life cycle stage</i>, excluding <a href="http://purl.obolibrary.org/obo/UBERON_0000071"<code>UBERON:0000071</code></a> for <i>death stage</i>.
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
      <td>Curator</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be a MONDO term or <a href="http://purl.obolibrary.org/obo/PATO_0000461"><code>"PATO:0000461"</code></a> for <i>normal</i> or <i>healthy</i>.
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
      <td>Curator</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be free-text that identifies a unique individual that data were derived from. It is STRONGLY RECOMMENDED that this identifier be designed so that it is unique to:<br><br>
          <ul><li>a given individual within the collection of datasets that includes this dataset</li>
          <li>a given individual across all collections in the cellxgene Data Portal</li></ul>
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
      <td>Curator</td>
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
      <td>Curator</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be a child of <a href="http://purl.obolibrary.org/obo/NCBITaxon_33208"<code>NCBITaxon:33208</code></a> for <i>Metazoa</i>.
        </td>
    </tr>
</tbody></table>
<br>

### self_reported_ethnicity_ontology_term_id

<table><tbody>
    <tr>
      <th>Key</th>
      <td>self_reported_ethnicity_ontology_term_id</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. If <code>organism_ontolology_term_id</code> is <code>"NCBITaxon:9606"</code> for <i>Homo sapiens</i>, this MUST be either a HANCESTRO term or <code>"unknown"</code> if unavailable. <br><br>Otherwise, for all other organisms this MUST be <code>"na"</code>.
        </td>
    </tr>
</tbody></table>
<br>

### sex_ontology_term_id

<table><tbody>
    <tr>
      <th>Key</th>
      <td>sex_ontology_term_id</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be a child of <a href="http://purl.obolibrary.org/obo/PATO_0001894">PATO:0001894</a> for  <i>phenotypic sex</i> or <code>"unknown"</code> if unavailable.
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
      <td>Curator</td>
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
              <td><i>10x transcription profiling</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0030080"><code>EFO:0030080</code></a>] and its children</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr> 
            <tr>
              <td><i>ATAC-seq</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0007045"><code>EFO:0007045</code></a>] and its children</td>
              <td><code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>CEL-seq2</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0010010"><code>EFO:0010010</code></a>]</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>CITE-seq</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0009294"><code>EFO:0009294</code></a>] and its children</td>
              <td><code>"cell"</code></td>
           </tr>
            <tr>
              <td><i>DroNc-seq</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0008720"><code>EFO:0008720</code></a>]</td>
              <td><code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>Drop-seq</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0008722"><code>EFO:0008722</code></a>]</td>
              <td><code>"cell"</code></td>
           </tr>
            <tr>
              <td><i>microwell-seq</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0030002"><code>EFO:0030002</code></a>]</td>
              <td><code>"cell"</code></td>
           </tr>    
            <tr>
              <td><i>Patch-seq</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0008853"><code>EFO:0008853</code></a>]</td>
              <td><code>"cell"</code></td>
           </tr>
            <tr>
              <td><i>sci-Plex</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0030026"><code>EFO:0030026</code></a>]</td>
              <td><code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>sci-RNA-seq</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0010550"><code>EFO:0010550</code></a>]</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>Seq-Well</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0008919"><code>EFO:0008919</code></a>]</td>
              <td><code>"cell"</code></td>
           </tr>
            <tr>
              <td><i>Smart-like</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0010184"><code>EFO:0010184</code></a>] and its children</td>
              <td><code>"cell"</code> or <code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>smFISH</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0009918"><code>EFO:0009918</code></a>] and its children</td>
              <td><code>"na"</code></td>
           </tr>   
            <tr>
              <td><i>snmC-seq</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0008939"><code>EFO:0008939</code></a>]</td>
              <td><code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>snmC-seq2</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0030027"><code>EFO:0030027</code></a>]</td>
              <td><code>"nucleus"</code></td>
           </tr>
            <tr>
              <td><i>spatial proteomics</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0700000"><code>EFO:0700000</code></a>] and its children</td>
              <td><code>"na"</code></td>
           </tr>
            <tr>
              <td><i>spatial transcriptomics by high-throughput sequencing</i> [<a href="http://www.ebi.ac.uk/efo/EFO_0030005"><code>EFO:0030005</code></a>] and its children</td>
              <td><code>"na"</code></td>
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
      <td>Curator</td>
    </tr>
   <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be the UBERON or CL term that best describes the tissue that this cell was derived from, depending on the type of biological sample:
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
              <td>Tissue</td>
              <td>STRONGLY RECOMMENDED to be an UBERON term<br>(e.g. <a href="http://purl.obolibrary.org/obo/UBERON_0008930"><code>"UBERON:0008930"</code></a> for a <i>sematosensory cortex</i> tissue sample)</td>
            </tr>
            <tr>
              <td>Cell Culture</td>
              <td>MUST be a CL term appended with <code>" (cell culture)"</code><br>(e.g. <code><a href="http://purl.obolibrary.org/obo/CL_0000057">"CL:0000057</a> (cell culture)"</code> for the <i>WTC-11 cell line</i>)</td>
            </tr>
            <tr>
              <td>Organoid</td>
              <td>MUST be an UBERON term appended with <code>" (organoid)"</code><br>(e.g. <code><a href="http://purl.obolibrary.org/obo/UBERON_0000955">"UBERON:0000955</a> (organoid)"</code> for a <i>brain organoid</i>)</td>
            </tr>
            <tr>
              <td>Enriched,<br>Sorted,or<br>Isolated<br>Cells from<br>a Tissue</td>
              <td>MUST be an UBERON or CL term and SHOULD NOT use terms that do not capture<br> the tissue of origin<br>(e.g. In the case of <i>CD3+ kidney cells</i>, use <a href="https://www.ebi.ac.uk/ols/ontologies/uberon/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FUBERON_0002113"><code>"UBERON:0002113"</code></a> for <i>kidney</i><br> instead of <a href="https://www.ebi.ac.uk/ols/ontologies/cl/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FCL_0000084"><code>"CL:000084"</code></a> for <i>T cell</i>. However, in the case of <i>EPCAM+ cervical cells</i>,<br>use <a href="https://www.ebi.ac.uk/ols/ontologies/cl/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FCL_0000066"><code>"CL:000066"</code></a> for <i>epithelial cell</i> of the cervix.)
              </td>
            </tr>
          </tbody></table>
        </td>
    </tr>
</tbody></table>
<br>

When a dataset is uploaded, the cellxgene Data Portal MUST automatically add the matching human-readable name for the corresponding ontology term to the `obs` dataframe. Curators MUST NOT annotate the following columns.

### assay

<table><tbody>
    <tr>
      <th>Key</th>
      <td>assay</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Data Portal</td>
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
      <td>Data Portal</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be the human-readable name assigned to the value of <code>cell_type_ontology_term_id</code>.
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
      <td>Data Portal</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be <code>"unknown"</code> if set in <code>development_stage_ontology_term_id</code>; otherwise, this MUST be the human-readable name assigned to the value of <code>development_stage_ontology_term_id</code>.
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
      <td>Data Portal</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be the human-readable name assigned to the value of <code>disease_ontology_term_id</code>.
        </td>
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
      <td>Data Portal</td>
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
      <td>Data Portal</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be <code>"na"</code> or <code>"unknown"</code> if set in <code>self_reported_ethnicity_ontology_term_id</code>; otherwise, this MUST be the human-readable name assigned to the value of <code>self_reported_ethnicity_ontology_term_id</code>.
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
      <td>Data Portal</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be <code>"unknown"</code> if set in <code>sex_ontology_term_id</code>; otherwise, this MUST be the human-readable name assigned to the value of <code>sex_ontology_term_id</code>.
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
      <td>Data Portal</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>categorical with <code>str</code> categories. This MUST be the human-readable name assigned to the value of <code>tissue_ontology_term_id</code>. <code>" (cell culture)"</code> or <code>" (organoid)"</code> MUST be appended if present in <code>tissue_ontology_term_id</code>.<br><br>
       For example, if the <code>tissue_ontology_term_id</code> was curated as <code>"CL:0000057 (cell culture)"</code>, then the value would be <code>"fibroblast (cell culture)"</code>.
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
      <td>Curator</td>
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
      <td>Curator</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>bool</code>. This MUST be <code>True</code> if the feature was filtered out in the normalized matrix (<code>X</code>) but is present in the raw matrix (<code>raw.X</code>). The value for all cells of the given feature in the normalized matrix MUST be <code>0</code>.  <br><br>Otherwise, this MUST be <code>False</code>. </td>
    </tr>
</tbody></table>
<br>


When a dataset is uploaded, cellxgene Data Portal MUST automatically add the matching human-readable name for the corresponding feature biotype, identifier, and the NCBITaxon term for the reference organism to the `var` and `raw.var` dataframes. Curators MUST NOT annotate the following columns:

### feature_biotype

<table><tbody>
    <tr>
      <th>Key</th>
      <td>feature_biotype</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Data Portal</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>This MUST be <code>"gene"</code> or <code>"spike-in"</code>.  
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
      <td>Data Portal</td>
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
      <td>Data Portal</td>
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
              <td><code>"NCBITaxon:9606"
              </code></td>
            </tr>
            <tr>
              <td><i>Mus musculus</i></td>
              <td><code>"NCBITaxon:10090"</code></td>
            </tr>
            <tr>
              <td><i>SARS-CoV-2</i></td>
              <td><code>"NCBITaxon:2697049"</code></td>
            </tr>
            <tr>
              <td><i>ERCC Spike-Ins</i></td>
              <td><code>"NCBITaxon:32630"</code></td>
            </tr>
          </tbody></table>
        </td>
    </tr>
</tbody></table>
<br>

## `obsm` (Embeddings)

For each `str` key, `obsm` stores a `numpy.ndarray` of shape `(n_obs, m)`, where `n_obs` is the number of rows in `X` and `m >= 1`.

To display a dataset in cellxgene Explorer, Curators MUST annotate **one or more** two-dimensional (`m >= 2`) embeddings (e.g. tSNE, UMAP, PCA, spatial coordinates) as `numpy.ndarrays` in `obsm`. The keys for these embedding MUST be prefixed with `"X_"`. The text that follows this prefix is presented to users in the *Embedding Choice* selector in cellxgene Explorer.

To illustrate, the [Krasnow Lab Human Lung Cell Atlas, 10X dataset](https://cellxgene.cziscience.com/e/krasnow_lab_human_lung_cell_atlas_10x-1-remixed.cxg/) in the [A molecular cell atlas of the human lung from single cell RNA sequencing collection](https://cellxgene.cziscience.com/collections/5d445965-6f1a-4b68-ba3a-b8f765155d3a) defines two embeddings in `obsm`:

* `"X_Compartment_tSNE"`
* `"X_tSNE"`

Users can then choose which embedding is visualized in cellxgene Explorer:

![Embeddings Illustration](images/embeddings.png
)

See also `default_embedding` in `uns`.

## `uns` (Dataset Metadata)

`uns` is a ordered dictionary with a `str` key. Curators MUST annotate the following keys and values in `uns`:

### schema_version

<table><tbody>
    <tr>
      <th>Key</th>
      <td>schema_version</td>
    </tr>
    <tr>
      <th>Annotator</th>
      <td>Curator</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          This MUST be <code>"3.0.0"</code>.
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
      <td>Curator</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          <code>str</code>. This text describes and differentiates the dataset from other datasets in the same collection. It is displayed on a page in the cellxgene Data Portal that also has the collection name. To illustrate, the first dataset name in the <a href="https://cellxgene.cziscience.com/collections/b52eb423-5d0d-4645-b217-e1c6d38b2e72">Cells of the adult human heart collection</a> is "All — Cells of the adult human heart".<br><br>It is STRONGLY RECOMMENDED that each dataset <code>title</code> in a collection is unique and does not depend on other metadata such as a different  <code>assay</code> to disambiguate it from other datasets in the collection.
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
      <td>Curator</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          <code>list[str]</code>. <code>str</code> values MUST refer to cell metadata keys in <code>obs</code>. Together, these keys define the <i>batches</i> that a normalization or integration algorithm should be aware of. For example if <code>"patient"</code> and <code>"seqBatch"</code> are keys of vectors of cell metadata, either <code>["patient"]</code>, <code>["seqBatch"]</code>, or <code>["patient", "seqBatch"]</code> are valid values.
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
      <td>Curator</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          <code>str</code>. The value MUST match a key to an embedding in <code>obsm</code> for the embedding to display by default in cellxgene Explorer.
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
      <td>Curator</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          <code>str</code>. cellxgene runs a heuristic to detect the approximate distribution of the data in X so that it can accurately calculate statistical properties of the data. This field enables the curator to override this heuristic and specify the data distribution explicitly. The value MUST be <code>"count"</code> (for data whose distributions are best approximated by counting distributions like Poisson, Binomial, or Negative Binomial) or <code>"normal"</code> (for data whose distributions are best approximated by the Gaussian distribution.)
        </td>
    </tr>
</tbody></table>
<br>

## Appendix A. Changelog

schema v3.0.0

* The canonical data format was updated from AnnData 0.7 to 0.8.
* All references to the "final" matrix has been replaced with "normalized" for clarity.
* General Requirements
  * Reserved Names from previous schema versions that have since been deprecated MUST NOT be present.
* obs (Cell metadata)
  * Removed guidance in `assay_ontology_term_id` that allowed clarifying text enclosed in parentheses if there was not an exact match for an assay.
  * Added `donor_id`
  * Renamed `ethnicity_ontology_term_id` to `self_reported_ethnicity_ontology_term_id`
  * Renamed `ethnicity` to `self_reported_ethnicity`
  * Added `suspension_type`
* var and raw.var (Gene metadata)
  * `feature_biotype` must be annotated by the Data Portal and not the Curator.
* uns (Dataset metadata)
  * Updated `schema_version`
  * Deprecated `X_normalization`

schema v2.0.0 substantially *remodeled* schema v1.1.0:

* "must", "should", and select other words have a defined, standard meaning.

* Curators are responsible for annotating ontology and gene identifiers. The cellxgene Data Portal adds the assigned human-readable names for all identifiers.

* Documented and *pinned* the required versions of ontologies and gene annotations used in schema validation.

* General Requirements
  * AnnData is now the canonical data format. The schema outline and descriptions are AnnData-centric.

  * Metazoan multi-organism data is accepted by the cellxgene Data Portal. For data that is neither Human, Mouse, nor SARS-COV-2, features MUST be translated into orthologous genes from the Human and Mouse gene annotations. 

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
