# Changelog
All notable changes to the python package `cellxgene-schema` are documented in this file.

The format of this changelog is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.1.2] - 2022-04-13

### Changed
- Feature names and reference organims id's are added to `raw.var`
- Duplicate genes with suffix "PAR_Y" are filtered
- Correct eroneous gene names introduced in GENCODE v38

## [2.1.1] - 2022-01-27

### Changed
- Gene reference files now have an extra column for gene length.
- Gene reference files no longer have transcripts.
- Improved CLI responsiveness in certain situations

## [2.1.0] - 2021-11-12    
### Added 
- A check for convertibility to Seurat format. Prints warnings if conversion is not possible.      

## [2.0.4] - 2021-11-01    
    
### Added    
    
- A check during raw data validation in `validate.Validator`. When both `adata.X` and `adata.raw.X` are present, but `adata.raw.X` contains non-integer values a warning is added.         

## [2.0.3] - 2021-10-08
### Changed
- Replace print with logging in `Validator` and `AnnDataLabelAppender`.
- Return the list of errors from `validate` so they can be used by callers.
- Add verbose option to CLI, it prints a progress log.

## [2.0.2] - 2021-10-07

### Changed

- Adds the assay ontology term EFO:0008939 (snmC-seq) to the list of assays for which validation of raw data is bypassed.

## [2.0.1] - 2021-10-06

### Changed

- Improves the performance of `Validator._is_raw()` in `validate.py`, and addresses timeouts when validating datasets with high number of non-zero values in the expression matrix.

## [2.0.0] - 2021-09-15

### Added

- All **MUST** requirements in [schema version 2.0.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/2.0.0/corpora_schema.md) are strictly enforced. For failures, an error message is displayed and validation fails. 
- Some **STRONGLY RECOMMENDED** requirements in the schema version 2.0.0 are checked. Warnings are displayed when recommended best practices are not observed.
- Pinned versions of ontology and feature references (see `./cellxgene_schema/ontology_files/`).
- Downloader and parser for ontology and feature references (see `./scripts/`).
- Option to apply human-readable labels for ontology and feature references: `cellxgene-schema validate --add-labels`.
- Tests that mirror the requirements in schema version 2.0.0.

### Changed

- `cellxgene-schema validate` validates schema version 2.0.0. The implementation is *from scratch*. 
- Ontology validation and label retrieval depend on downloaded references instead of the EBI Ontology Service (see `./cellxgene_schema/ontology.py`).
- Gene/feature validation and label retrieval depend on downloaded references (see `./cellxgene_schema/ontology.py`).

### Removed

- Subcommand `cellxgene-schema apply`.
- Support for schema versions 1.x.x.
