# Changelog
All notable changes to the python package `cellxgene-schema` are documented in this file.

The format of this changelog is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
