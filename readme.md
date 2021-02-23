# cellxgene curation tools

This repository contains documents and code used by cellxgene's curation team. Issues/suggestions pertaining to datasets and how they interact with cellxgene should be created here. 

For information/issues about cellxgene and its portal please refer to:

- [cellxgene](https://github.com/chanzuckerberg/cellxgene)
- [Corpora data portal](https://github.com/chanzuckerberg/corpora-data-portal)

## Installation

The central tool provided here is a CLI for augmenting datasets with the [cellxgene schema](docs/corpora_schema.md) so they can be hosted at [cellxgene's portal](https://cellxgene.cziscience.com/). 

It is available through pip:

```
pip install cellxgene-schema
```

It can also e installed from the source by cloning this repository and running:

```
make install 
```

And you can run the test with:

```
make unit-test
```


## Quick start

The CLI augments  an [AnnData file](https://anndata.readthedocs.io/en/latest/) (\*.h5ad) with cellxgene schema required ontology terms using the logic defined in a yaml config file. This yaml file should indicate the values for the schema slots or mappings between the original values and the corresponding schema slots.

An example of a yaml config file looks like this:

```
obs:
  assay_ontology_term_id: EFO:0010550
  ethnicity_ontology_term_id: unknown
  sex: male
  tissue_ontology_term_id: UBERON:0000970
  cell_type_ontology_term_id:
    sub_cluster_name:
      Adrenocortical cells-1: CL:0002097
      Photoreceptor cells-1: CL:0000210
uns:
  version:
    corpora_schema_version: 1.1.0
    corpora_encoding_version: 0.1.0
  organism: Homo sapiens
  organism_ontology_term_id: NCBITaxon:9606
  layer_descriptions:
    X: log1p
   raw.X: raw
  publication_doi: http://dx.doi.org/10.1126/science.aba7721
  title: Survey of human embryonic development
fixup_gene_symbols:
    X: log1p
   raw.X: raw
```

You can use the config file to augment a dataset with  the schema using:

```
cellxgene-schema apply --source-h5ad original.h5ad --remix-config config.yml --output-filename remixed.h5ad
```

And then verify that the schema was properly added with:

```
cellxgene-schema validate remixed.h5ad
```

A detailed manual for the CLI and the config yaml file can be found [here](docs/schema_guide.md).

## Datasets curated by cellxgene’s curation team

Scripts demonstrating how the cellxgene team has curated datasets for hosting on the portal are stored in this repository. They provide worked examples that provide additional demonstrations of how the tool can be used.

The `datasets` folder contains step-by-step curation instructions for each dataset we have curated, each dataset has its own independent folder and readme. 
In principle anyone could reproduce our curation process following the dataset's readme, which starts from downloading data (usually publicly available) and finishes by creating one or more `*.h5ad` files that follow cellxgene schema and are ready to be hosted at cellxgene's portal.

The `docs` folder contains guides, general documentsconf, files, or scripts that have been used or could be used in the future for curation or integration processes.

## Contributing

Please read our contributing [guidelines](CONTRIBUTING.md) and make sure adhere to the Contributor Covenant [code of conduct](https://github.com/chanzuckerberg/.github/blob/master/CODE_OF_CONDUCT.md). 

## Reporting Security Issues                     
                                                
Please read our [security reporting policy](SECURITY.md)
