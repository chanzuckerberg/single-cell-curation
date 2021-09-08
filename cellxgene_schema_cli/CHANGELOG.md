# Changelog
All notable changes to the python package `cellxgene-schema` will be documented in this file.

The format of this changelog is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.0.0] - 2021-09-09
### Added
- All MUSTs in schema specification are validated, if any is violated a specific error is shown. Refer to the [schema specification](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/2.0.0/corpora_schema.md) for specifc changes.
- Some SHOULDs in the schema specification are validated, if any is violated a specific warning is shown. Refer to the [schema specification](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/2.0.0/corpora_schema.md) for specifc changes.
- Reference ontology files  (`./cellxgene_schema/ontology_files/`).
- Reference feature files (`./cellxgene_schema/ontology_files/`).
- Downloader and parser for ontology and feature files
- `cellxgene-schema validate --add-labels` (`./scripts/`).
- Tests that mirror the schema specification.


### Changed
- All code related to `cellxgene-schema validate` was re-written.
- `cellxgene-schema validte` validates schema version 2.0.0.
- Ontology validation and label retrieval are done using pinned reference files instead of EBI's REST API (`./cellxgene_schema/ontology.py`).
- Gene/feature validation and label retrieval are done through pinned reference files (`./cellxgene_schema/ontology.py`).

### Removed
- Subcommand `cellxgene-schema apply` was deprecated.
- Schema versions 1.x.x are no longer supported.
