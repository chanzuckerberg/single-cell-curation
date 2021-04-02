
## Overview

*— it made sense until I reached a period. Then what I thought I had comprehended vanished.*

**Ben Lerner – The Reflections of a Reading**

---

Curators may add gene sets to their data collections hosted on the [cellxgene data portal](https://cellxgene.cziscience.com), by uploading a file that is *properly formatted* based on the requirements in this document.

When the portal successfully validates a gene set file, individual gene sets are extracted and stored in the database. If validation fails for an individual gene set in the file, then the entire operation fails. No gene sets from the file are stored in the database.

The portal **does not** store the uploaded file. 

The order of the gene symbols in each gene set is maintained in the database. (The order of gene symbols in the download for a specific gene set name matches the order of the gene symbols in the upload for a specific gene set name.)

Multiple gene set files may be uploaded to the portal by the curator. 

The key words "MUST", "MUST NOT", "REQUIRED", "SHALL", "SHALL NOT", "SHOULD", "SHOULD NOT", "RECOMMENDED", "NOT RECOMMENDED" "MAY", and "OPTIONAL" in this document are to be interpreted as described in BCP 14, [RFC2119](https://www.rfc-editor.org/rfc/rfc2119.txt), and [RFC8174](https://www.rfc-editor.org/rfc/rfc8174.txt) when, and only when, they appear in all capitals, as shown here.

## cellxgene gene set data format
  
The cellxgene gene set data format is a [*Tidy* CSV](./gene_sets_example.csv) (comma-separated values) file using ASCII encoding. Multiple gene sets MAY be included in the file similar to the [Gene Matrix Transposed](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) or [Gene Matrix](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMX:_Gene_MatriX_file_format_.28.2A.gmx.29) formats.



Example:

| gene_set_name | gene_set_description | gene_symbol | gene_description | provenance1      | provenance1_description |
|---------------|----------------------|-------------|------------------|------------------|-------------------------|
| club.cell     | description          | CCKAR       | description      | Pubmed ID XYZ123 | Primary Pubmed ID       |
| club.cell     |                      | SCGB3A2     | description      | Pubmed ID ABC456 | Primary Pubmed ID       |
| club.cell     |                      | CYP2F2      | description      | Pubmed ID DCF678 | Primary Pubmed ID       |
| macrophage    | description          | CD68        | description      |                  |                         |
| macrophage    |                      | CD163       | description      |                  |                         |


## Mandatory Header

The first row MUST contain header columns using reserved names in the following order:

* `gene_set_name`
* `gene_set_description`
* `gene_symbol`
* `gene_description`

Publishers MAY include additional header columns. It is RECOMMENDED that these custom columns observe the same self-documenting style - `column_name` and `column_description` for easier comprehension by data consumers.

## Rows

Each subsequent row describes a gene in a gene set.
<br><br>

### Values

The values for `gene_set_name`, `gene_set_description`, and `gene_symbol` MUST NOT contain illegal ASCII characters or sequences. If the following cases are detected, validation MUST display an error message and fail the upload:

* control characters (decimal 0-31)
* DEL (decimal 127)
* leading spaces (decimal 32) in a field - "     This is an example"
* trailing spaces (decimal 32) in a field - "This is an example     " 
* multiple spaces (decimal 32) "internal" to a field - "This     is an example"

### `gene_set_name`

---

The `gene_set_name` column MUST contain a value.

If the `gene_set_name` is missing, validation MUST display an error message and fail the upload. This is illustrated by **~~?~~** in the example:

| gene_set_name | gene_set_description | gene_symbol | gene_description |
|---------------|----------------------|-------------|------------------|
| club.cell     | description          | CCKAR       | description      |
|     **~~?~~** |                      | SCGB3A2     | description      |
| club.cell     |                      | CYP2F2      | description      |


`gene_symbol(s)` for a `gene_set_name` MAY exist on noncontiguous rows. An out-of-order `gene_symbol` MUST be added to an existing `gene_set_name`. In the example below, **~~CD163~~** is added to the **club.cell** gene set:

---
| gene_set_name | gene_set_description | gene_symbol | gene_description |
|---------------|----------------------|-------------|------------------|
| club.cell     | description          | CCKAR       | description      |
| club.cell     |                      | SCGB3A2     | description      |
| club.cell     |                      | CYP2F2      | description      |
| macrophage    | description          | CD68        | description      |
| club.cell     |                      | **~~CD163~~**       | description      |
<br>

When new gene sets are being added to a data collection on the portal, validation MUST detect `gene_set_name` collisions with current gene sets in the collection, display an error message, and fail the upload.  <br><br>

### `gene_set_description`

---

The first instance of a `gene_set_description` column for a specific `gene_set_name` MUST contain a value. All other instances are ignored in subsequent rows for the same `gene_set_name`.

If the first instance of the `gene_set_description` is missing, validation MUST display an error message and fail the upload. This is illustrated by **~~?~~** in the example:

---

| gene_set_name | gene_set_description | gene_symbol | gene_description |
|---------------|----------------------|-------------|------------------|
| club.cell     |        **~~?~~**              | CCKAR       | description      |
| club.cell     |                      | SCGB3A2     | description      |
| club.cell     |                      | CYP2F2      | description      |


### `gene_symbol`

---

The `gene_symbol` column MUST contain a value that is unique for the `gene_set_name`. Validation MUST display an error message for a duplicate `gene_symbol` and fail the upload. This is illustrated by **~~CCKAR~~** in the example:

| gene_set_name | gene_set_description | gene_symbol | gene_description |
|---------------|----------------------|-------------|------------------|
| club.cell     | description          | CCKAR       | description      |
| club.cell     |                      | SCGB3A2     | description      |
| club.cell     |                      | **~~CCKAR~~**   | description      |


The value for a `gene_symbol` SHOULD follow [cellxgene schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/docs/corpora_schema.md) guidance for gene symbols. <br><br>

### `gene_description`

---

The `gene_description` column MAY contain a value. <br><br>

## Presentation in the cellxgene UX

`gene_set_name` and `gene_set_description` are presented to data consumers  viewing data collections in the portal. They also may download gene sets.

Users can [**explore**](https://cellxgene.cziscience.com/e/6acb6637-ac08-4a65-b2d1-581e51dc7ccf.cxg/) a dataset and its related gene sets in the data collection. In the visualization, `gene_set_name`, `gene_set_description`,  `gene_symbol`, and `gene_description` are presented to the user. The user may color by the mean expression of the gene set, select cells from a histogram showing the distribution of mean expression, and plot gene sets on the scatter plot.
