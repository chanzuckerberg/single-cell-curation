# cellxgene Data Integration Schema AnnData Implementation
​
Authors: acarr@chanzuckerberg.com
​
Document Status: _Approved_
​
Version: 1.1.1
​
Date Last Modified: 2021-04-26
​
cellxgene Data Portal datasets are stored in the HDF5-backed AnnData format as written by version 0.7 of the anndata library.  
One goal of using this format is to have all metadata stored in the same file as the data, so cellxgene will have access to all metadata within a single file.
​
The schema requirements and definitions for the AnnData `X`, `uns`, `obs`, and `obsm` attributes are below.
​
An h5ad that follows the schema will have a `version` attribute in `uns` which is a `dict`, and that dict will
have a key "corpora_schema_version".
​
### Note on types
The types below are python3 types.  
Note that a python3 `str` is a sequence of Unicode code points, which is stored null-terminated and UTF-8-encoded by anndata.
​
## `X`
​
cellxgene does not impose any additional constraints on the `X` data matrix.  
So it may be sparse or dense and any numeric `numpy.dtype`.
​
The cellxgene schema has the same layer requirements for different scRNA-seq assay types:

| Assay                                 | "raw" required? | "raw" location                                                                     | "final" required?     | "final" location | Other layers |
|---------------------------------------|-----------------|------------------------------------------------------------------------------------|-----------------------|------------------|--------------|
| scRNA-seq (UMI, e.g. 10x v3)          | REQUIRED        | AnnData.raw or AnnData.layers["raw"] unless no "final" is provided, then AnnData.X | STRONGLY RECOMMENDED  | AnnData.X        | OPTIONAL     |
| scRNA-seq (non-UMI, e.g. SS2)         | REQUIRED        | AnnData.raw or AnnData.layers["raw"] unless no "final" is provided, then AnnData.X | STRONGLY RECOMMENDED  | AnnData.X        | OPTIONAL     |

## `uns`
​
Recall that `uns` is a mapping with `str`s as keys. The cellxgene data integration schema requires the following keys and values in `uns`:
​
**Key**|**Value Type**|**Notes**
-----|-----|-----
version|`dict`|Contains keys "corpora_schema_version" and "corpora_encoding_version" indicating the version of the schema and schema encoding used to create the h5ad.
title|`str`|
layer\_descriptions|`dict` with string keys and values|One key must be "X" which describes the transformations (if any) performed to produce the X matrix cellxgene displays.
organism|`str`|
organism\_ontology\_term\_id|`str`|
batch\_condition|`str` OR `list[str]`. `str` values must refer to keys in `obs`.

​
The cellxgene schema defines the following optional keys and values in `uns`.  
If the key is present, then the value must not be empty.
​
**Key**|**Value Type**|**Notes**
-----|-----|-----
default\_embedding|`str`|Must match a key to an embedding in `obsm`.
default\_field|`str`|Must match a column name in `obs`
tags|`list` of `str`s|
<obs\_column>\_colors|`list` of color specifications (see anndata/cellxgene documentation)|<obs\_column> must be a column name in `obs`. There may be multiple keys like this.
​
## `obsm`
​
`obsm` is a mapping from `str`s to numpy arrays of shape `(n_obs, m)` where `n_obs` is the number of rows in `X` and `m >= 1`.
cellxgene requires one value in `obsm` to be an at least two-dimensional embedding, meaning `m >= 2`, and the key for that value must begin with `X_`.
​
## `obs`
​
The cellxgene schema requires a number of columns be present in the `obs` dataframe:
​
**Key**|**Value Type**|**Notes**
-----|-----|-----
tissue|`str` or categorical with `str` categories|
assay|`str` or categorical with `str` categories|
disease|`str` or categorical with `str` categories|
cell\_type|`str` or categorical with `str` categories|
sex|`str` or categorical with `str` categories|
ethnicity|`str` or categorical with `str` categories|
development\_stage|`str` or categorical with `str` categories|
tissue\_ontology\_term\_id|`str` or categorical with `str` categories|Not necessarily human-readable
assay\_ontology\_term\_id|`str` or categorical with `str` categories|Not necessarily human-readable
disease\_ontology\_term\_id|`str` or categorical with `str` categories|Not necessarily human-readable
cell\_type\_ontology\_term\_id|`str` or categorical with `str` categories|Not necessarily human-readable
ethnicity\_ontology\_term\_id|`str` or categorical with `str` categories|Not necessarily human-readable
development\_stage\_ontology\_term\_id|`str` or categorical with `str` categories|Not necessarily human-readable
​
## `var`
​
cellxgene requires unique feature identifiers, so the index of `var` must not contain any duplicate values.
Moreover, `var.index` must contain the human-readable display names for features, for example HGNC symbols.

## Appendix A. Changelog

* [#43](https://github.com/chanzuckerberg/single-cell-curation/issues/43) obsm MUST be encoded as numpy arrays
* AnnData MUST contain a raw layer.
