# cellxgene internal curation of hosted data

This repository contains documents and code used by cellxgene's curation team. Issues/suggestions pretaining to datasets and how they interact with cellxgene may be created here. 

For information/issues about cellxgene and its portal please refer to:

- [cellxgene](https://github.com/chanzuckerberg/cellxgene)
- [Corpora data portal](https://github.com/chanzuckerberg/corpora-data-portal)

# Getting started

The main set of tools hosted in this repo is for injecting the [cellxgene schema](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema_h5ad_implementation.md) into single-cell data to be hosted at [cellxgene's portal](https://cellxgene.cziscience.com/). 

We provide a CLI with one command to inject the schema (`cellxgene schema apply)` and one for validation (`cellxgene schema validate`).

Please refer to the [manual](docs/schema_guide.md) for detailed usage instructions.

# Repository organization

The folder structure is:

## datasets

Contains step-by-step curation instructions for each dataset we have curated, each dataset has its own independent folder and readme. 

In principle anyone could reproduce our curation process following the dataset's readme, which starts from downloading data (usually publicly available) and finishes at creating one or more `*.h5ad` files that follow [cellxgene schema](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema_h5ad_implementation.md) and are ready to be hosted at [cellxgene's portal](https://cellxgene.cziscience.com/).

## docs

General documents, files, or scripts that have been used or could be used in the future for curation or integration processes.

## Schema CLI Tool

In the `cellxgene_schema_cli` folder, there is a command line tool for applying and validating the cellxgene integration
schema. In can be installed with

```
make install
```

and you can run the tests with

```
make unit-test
```

# Contributing

Please read our contributing [guidelines](CONTRIBUTING.md) and make sure adhere to the Contributor Covenant [code of conduct](https://github.com/chanzuckerberg/.github/blob/master/CODE_OF_CONDUCT.md).

# Reporting Security Issues

Please read our [security reporting policy](SECURITY.md)


