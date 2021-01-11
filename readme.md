# cellxgene internal curation of hosted data

This repository contains documents and code used by cellxgene's curation team. Issues/suggestions pretaining to datasets and how they interact with cellxgene may be created here. 

The folder structure is:

### datasets

Contains step-by-step curation instructions for each dataset we have curated, each dataset has he's own independent folder and readme. 

In principle anyone could reproduce our curation process following the dataset's readme, which starts from downloading data (usually publicly available) and finishes at creating one or more `*.h5ad` files that follow [cellxgene schema](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema_h5ad_implementation.md) and are ready to be hosted at [cellxgene's portal](https://cellxgene.cziscience.com/).

### docs

General documents, files, or scripts that have been used or could be used in the future for curation or integration processes.
