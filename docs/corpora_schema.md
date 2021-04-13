# cellxgene Data Integration Schema

Contact: acarr@chanzuckerberg.com

Document Status: _Approved_

Version: 1.1.0

Date Last Modified: 2020-12-11


## Background

cellxgene aims to support the publication, sharing, and exploration of single-cell datasets.
Building on those published datasets, cellxgene seeks to create references of the phenotypes and composition of cells that make up human tissues.
Creating references from multiple datasets requires some harmonization of metadata and features in the cellxgene Data Portal. But if that harmonization is too onerous, it
will burden the goal of rapid data sharing.

We balance publishing and reference creation needs by requiring datasets hosted in the cellxgene Data Portal to follow a small schema with only a few required fields. These fields
are expected to be very useful for data integration and also simply and readily available from data submitters.

Note that the requirements in the schema are just the minimum required information. Datasets often have additional metadata, and this
is preserved in datasets submitted to the Data Portal.

## Curation and Validation

When a submitter is preparing a dataset for the cellxgene Data Portal, they can use the `cellxgene-schema apply` command line tool to apply changes to their dataset so that it
follows the schema. When that tool successfully completes, it will give the submitter a message that the dataset is ready to be submitted to the Data Portal,
and it will write the schema version number into the dataset's metadata.

Then, when the submitter uploads the dataset to the Data Portal, the Data Portal will verify that the dataset does indeed follow an accepted version of the
schema. If it does not, it will reject the dataset with an appropriate error message.

## Schema

### Basic Requirements

*   **Unique feature identifiers**. Every gene feature MUST be assigned a unique identifier. This is occasionally not present because
    of one-to-many mappings between gene symbols and other gene ids. In cases where there are duplicated feature identifiers, they will need to be appropriately
    combined before submission. For example, raw counts will be summed and logged counts will be exponentiated, summed, and logged.
    This is needed for the explorer to function. When a user requests gene expression information, the explorer needs to be able to unambiguously return a
    single value.
*   **No PII**. Curators agree to this requirement as part of the data submission policy. However, it is not strictly enforced in our validation tooling because it
*   is difficult for software to predict what is and is not PII. It is up to the submitter to ensure that no metadata can be personally identifiable: no names,
*   dates of birth, specific locations, etc. There's a [list](https://docs.google.com/document/d/1nlUuRiQ8Awo_QsiVFi8OWdhGONjeaXE32z5GblVCfcI/edit?usp=sharing).

### Matrix Layers

cellxgene's data requirements are tailored to optimize data reuse. Because each assay has different characteristics, our requirements differ by assay type. In general,
cellxgene requires submission of "raw" data suitable for computational reuse, when a standard form exists, and strongly suggests that a "final" matrix suitable for
visualization in the explorer be included. cellxgene uses `AnnData` for data ingestion, and uses `AnnData`'s `layers` functionality to accept multiple matrix layers,
such as "raw" and "final". This imposes some requirements on data of all assay types:

*   Anndata provides a [raw](https://anndata.readthedocs.io/en/latest/anndata.AnnData.raw.html#anndata.AnnData.raw) attribute. If a raw matrix is required for an assay
    type, it SHOULD be stored in `AnnData.raw`, but it MAY be instead stored in `AnnData.layers["raw"]`.
*   AnnData requires that matrix layers MUST have the same dimension, so raw count matrices MUST include the same cells and genes as the final. Because it is impractical to
    retain all barcodes in raw and final matrices, we require that cells be filtered from both. By contrast, those wishing to reuse datasets require access to raw gene expression values,
    so we require that genes are present in both datasets. Summarizing, any cell barcodes that are removed from the data MUST be filtered from both raw and final matrices.
    By contrast, genes MUST NOT be filtered from the raw matrix. Any genes that publishers wish to filter from the final matrix MAY have their values replaced by np.nan,
    which will mask them from exploration.

In addition to these general requirements, the following table describes the matrix data and layers requirements that are assay-specific. If cellxgene does not support
an assay you would like to publish, please post an issue on this repository to start a conversation about extending the schema. If an entry in the table is empty,
the cellxgene schema does not have any other requirements on data in those layers beyond the ones listed above. This is usually the case when there are many ways to
produce the matrix layer in question.

| Assay                                 | Raw requirements                                                                                                                              | Raw Required? | Final Requirements                                                                                 | Final Required?      | Other Layers |
|---------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------|---------------|----------------------------------------------------------------------------------------------------|----------------------|--------------|
| scRNA-seq (UMI, e.g. 10x v3)          | Values MUST be de-duplicated molecule counts.                                                                                                 | REQUIRED      |                                                                                                    | STRONGLY RECOMMENDED | OPTIONAL     |
| scRNA-seq (non-UMI, e.g. SS2)         | Values MUST be one of read counts (e.g. FeatureCounts) or  estimated fragments (e.g. output of RSEM).                                         | REQUIRED      |                                                                                                    | STRONGLY RECOMMENDED | OPTIONAL     |
| Accessibility (e.g. ATAC-seq, mC-seq) |                                                                                                                                               | NOT REQUIRED  | Values MUST correspond to gene features (e.g. HGNC symbols if human)                               | REQUIRED             | OPTIONAL     |

### Schema Version

Datasets in the Data Portal MUST store the version of the schema they follow (that is, the version of this document) as
well as the version of the particular encoding used. The encoding is documented
[elsewhere](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema_h5ad_implementation.md) and describes techincal details
of how the schema should be serialized in a particular file format.

**Field name**|**Constraints**
:--|:--
corpora_schema_version|Follows [SemVer](https://semver.org/) conventions
corpora_encoding_version|Follows [SemVer](https://semver.org/) conventions

### Integration Metadata

To support data integration, each cell MUST have the following metadata:

**Field name**|**Constraints**
:--|:--
tissue|string. This field SHOULD be appended with " (cell culture)" or " (organoid)" if appropriate.
assay|string
disease|string
cell\_type|string.
sex|"male", "female", "mixed", "unknown", or "other"
ethnicity|string, "na" if non-human, "unknown" if not available
development\_stage|string, "unknown" if not available

In addition to these free text fields (except sex), each cell MUST also have ontology values:

**Field name**|**Constraints**
:--|:--
tissue\_ontology\_term\_id|UBERON term
assay\_ontology\_term\_id|EFO term
disease\_ontology\_term\_id|MONDO term or [PATO:0000461](http://bioportal.bioontology.org/ontologies/PATO?p=classes&conceptid=PATO%3A0000461)
cell\_type\_ontology\_term\_id|CL term
ethnicity\_ontology\_term\_id|HANCESTRO term, "na" if non-human
development\_stage\_ontology\_term\_id|HsapDv term if human, child of EFO:0000399 otherwise

The `tissue_ontology_term_id` field must be appended with " (cell culture)" or " (organoid)" if appropriate.

cellxgene requires ontology terms to enable search, comparison, and integration of data. When no appropriate ontology value is available, then the most
precise accurate term MUST be used. For example if the `cell_type` field describes a relay interneuron, but the most specific available term in the CL
ontology is CL:0000099 ("Interneuron"), then the interneuron term can be used to fulfill this requirement, and ensures that users searching for "neuron"
are able to find these data. Users will still be able to access any additional more specific cell type annotations that have been submitted with the
data (but aren't required by the schema).

Ontology terms MUST use [OBO-format ID](http://www.obofoundry.org/id-policy.html), meaning they are a CURIE
where the prefix identifies the ontology. For example `EFO:0000001` is a term in the `EFO` ontology.

If the features of the dataset are human genes, then the feature ids MUST be [HGNC](https://www.genenames.org/about/guidelines/#!/#tocAnchor-1-7) approved
gene symbols.

Similarly if the features are mouse genes, then the feature ids SHOULD be [MGI](http://www.informatics.jax.org/mgihome/nomen/gene.shtml) gene symbols.

Finally, the whole dataset MUST be annotated with fields that indicate the organism and describe the meanings of the [matrix layers](#Matrix-Layers):

**Field name**|**Constraints**
:--|:--
organism|String
organism\_ontology\_term\_id|NCBITaxon term
layer\_descriptions|Mapping from {layer\_name: layer\_description, ...} Each description is free text, though one layer must be described as "raw".


### Presentation Metadata

There are also fields that are required so that the cellxgene Data Portal and Explorer can present datasets appropriately.

* Each dataset MUST have at least one **embedding**, a mapping from each cell to a tuple of floats of length at least 2. These are usually generated by algorithms
  like umap or tsne, but can also coordinates of cells in spatial assays. They are used to display the dataset in the Explorer.
* Each dataset MUST have a title. This is a string that describes and differentiates the dataset from others in the collection. It will be displayed on a page that
  also has the collection name. For example, in the collection
  [Cells of the human heart](https://cellxgene.cziscience.com/collections/b52eb423-5d0d-4645-b217-e1c6d38b2e72), the first dataset name is "All — Cells of the
  adult human heart".

#### Presentation Hints

The metadata fields below are optional. They aren't needed for integration, and cellxgene can display the data fine without them, but if they are
included cellxgene will do something with them. This allows submitters to fine-tune how their datasets are presented, which is a common request.


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

### **Scanpy/AnnData**

The Data Portal requires submitted count matrices and associated metadata to be in [Anndata](https://anndata.readthedocs.io/en/stable/) format, the hdf5-based
file format used by scanpy. There is a python library for interacting with it. The count data is stored in an attribute `X` of shape (# of cells, # of genes).
`X` can either be a numpy.ndarray or a scipy.sparse.spmatrix.

Information about cells and genes are stored in `obs` and `var` dataframes, respectively. Each of those has an index which can serve as cell and gene ids
if all the values are unique.

The `obs` dataframe must store all the cell-level metadata except for embeddings. So, there must be columns named "tissue", "assay", "disease", and
"cell_type".

Embeddings are stored in `obsm` with a key name of "X_{description}" where description  provides some information about how the embedding was generated,
X_umap for example.

Finally, the dataset-level metadata is stored in `uns`, which is just key-value pairs.

Note that anndata supports "layers" and "raw" values for counts. Those are permitted, but cellxgene will treat `X` as the "final" matrix for further
visualization. Once anndata [unifies its treatment of layers](https://github.com/theislab/anndata/issues/244), cellxgene will use the "default" as the
final matrix, however that ends up being specified. In anndata files, the layer_descriptions dictionary should have a key "X" and optionally "raw.X" to
describe those layers.
