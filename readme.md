# cellxgene curation tools

This repository contains documents and code used by cellxgene's curation team. Issues/suggestions pertaining to datasets and how they interact with cellxgene should be created here. 

For information/issues about cellxgene and its portal please refer to:

- [cellxgene](https://github.com/chanzuckerberg/cellxgene)
- [Single-cell data portal](https://github.com/chanzuckerberg/single-cell-data-portal)

## Installation

The primary curation tool is the `cellxgene-schema` CLI. It enables curators to perform [schema](single-cell-curation/schema/2.0.0/corpora_schema.md) validation for datasets to be hosted on the [cellxgene Data Portal](https://cellxgene.cziscience.com/).

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

## Usage

The CLI validates an [AnnData file](https://anndata.readthedocs.io/en/latest/) (\*.h5ad) to ensure that it addresses the schema requirements.

Datasets can be validated using the following command line:

```
cellxgene-schema validate input.h5ad
```

If the validation succeeds, the command returns a zero exit code; otherwise, it returns a non-zero exit code and prints validation failure messages.


---

The data portal runs the following in the backend:

```
cellxgene-schema validate --add-labels output.h5ad input.h5ad
```

This execution validates the dataset as above AND adds the human-readable labels for the ontology and gene IDs as defined in the schema. If the validation is successful, a new AnnData file (output.h5ad) is written to disk with the labels appended.

This option SHOULD NOT be used by data contributors.


## Datasets curated by cellxgene’s curation team

Scripts demonstrating how the cellxgene team has curated datasets for hosting on the portal are stored in this repository. They provide worked examples that provide additional demonstrations of how the tool can be used.

The `datasets` folder contains step-by-step curation instructions for each dataset we have curated, each dataset has its own independent folder and readme. 
In principle anyone could reproduce our curation process following the dataset's readme, which starts from downloading data (usually publicly available) and finishes by creating one or more `*.h5ad` files that follow cellxgene schema and are ready to be hosted at cellxgene's portal.

The `docs` folder contains guides, general documents, files, or scripts that have been used or could be used in the future for curation or integration processes.

## Contributing

Please read our contributing [guidelines](CONTRIBUTING.md) and make sure adhere to the Contributor Covenant [code of conduct](https://github.com/chanzuckerberg/.github/blob/master/CODE_OF_CONDUCT.md). 

## Reporting Security Issues                     
                                                
Please read our [security reporting policy](SECURITY.md)
