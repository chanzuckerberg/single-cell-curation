# cellxgene curation tools

This repository contains documents and code used by cellxgene's curation team. Issues/suggestions pertaining to datasets and how they interact with cellxgene should be created here. 

For information/issues about cellxgene and its portal please refer to:

- [cellxgene](https://github.com/chanzuckerberg/cellxgene)
- [Corpora data portal](https://github.com/chanzuckerberg/corpora-data-portal)

## Installation

The central tool provided here is a CLI that validaties datasets follow the [cellxgene schema](single-cell-curation/schema/2.0.0/corpora_schema.md) so they can be hosted at [cellxgene's portal](https://cellxgene.cziscience.com/). 

It is available through pip:

```
pip install cellxgene-schema
```

It can also be installed from the source by cloning this repository and running:

```
make install 
```

And you can run the tests with:

```
make unit-test
```


## Quick start

The CLI validates an [AnnData file](https://anndata.readthedocs.io/en/latest/) (\*.h5ad) based on the cellxgene schema specifications.

You can run the validation by doing the following:

```
cellxgene-schema validate input.h5ad
```

If the validation is succesful there will be a zero exit status, otherwise there will be error messages indicating why validation was unsuccesful  along with a non-zero exit status.

A detailed manual for the CLI can be found [here](docs/schema_guide.md).

## Datasets curated by cellxgene’s curation team

Scripts demonstrating how the cellxgene team has curated datasets for hosting on the portal are stored in this repository. They provide worked examples that provide additional demonstrations of how the tool can be used.

The `datasets` folder contains step-by-step curation instructions for each dataset we have curated, each dataset has its own independent folder and readme. 
In principle anyone could reproduce our curation process following the dataset's readme, which starts from downloading data (usually publicly available) and finishes by creating one or more `*.h5ad` files that follow cellxgene schema and are ready to be hosted at cellxgene's portal.

The `docs` folder contains guides, general documentsconf, files, or scripts that have been used or could be used in the future for curation or integration processes.

## Contributing

Please read our contributing [guidelines](CONTRIBUTING.md) and make sure adhere to the Contributor Covenant [code of conduct](https://github.com/chanzuckerberg/.github/blob/master/CODE_OF_CONDUCT.md). 

## Reporting Security Issues                     
                                                
Please read our [security reporting policy](SECURITY.md)
