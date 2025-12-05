# cellxgene curation tools

[![codecov](https://codecov.io/gh/chanzuckerberg/single-cell-curation/branch/main/graph/badge.svg?token=J8OT7OXKHJ)](https://codecov.io/gh/chanzuckerberg/single-cell-curation)

This repository contains documents and code used by cellxgene's curation team. Issues/suggestions pertaining to datasets and how they interact with cellxgene should be created here. 

For information/issues about cellxgene and its portal please refer to:

- [single-cell-explorer](https://github.com/chanzuckerberg/single-cell-explorer)
- [single-cell-data-portal](https://github.com/chanzuckerberg/single-cell-data-portal)

## Installation

The primary curation tool is the `cellxgene-schema` CLI. It enables curators to perform [schema](./schema/3.0.0/schema.md) validation for datasets to be hosted on the [cellxgene Data Portal](https://cellxgene.cziscience.com/).

It requires Python >= 3.8. It is available through pip:

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

This experimental validator also offers the option to annotate required columns `cell_type_ontology_term_id` and `tissue_ontology_term_id` in Zebrafish, Fruit Fly, or C. Elegans anndata BEFORE running validation commands above. 

This relies on your anndata having the appropriate species-specific ontology terms (e.g. ZFA, FbBT, WBbt) labeled in `organism_cell_type_ontology_term_id` and `organism_tissue_ontology_term_id`, respectively. 

```
cellxgene-schema map-species output.h5ad input.h5ad
```

The command will find the closest CL (for cell_type) or UBERON (for tissue) mapping for the given term, offering either an exact match for the given term or a match from the closest possible ancestor term. This is based on the CL and UBERON [SSSOM](https://mapping-commons.github.io/sssom/toolkit/) mappings.

If there are multiple closest ancestors of the same distance with a match, the command will NOT annotate those rows and instead log your closest ancestor match options for your manual curation.

---

This experimental validator also offers the option to annotate genetic perturbations with genomic locations and target genes. This is useful for CRISPR perturbation datasets where you want to add `target_genomic_regions` and `target_features` annotations.

This command requires [guidescan2](https://github.com/pritykinlab/guidescan-cli) to be installed on your system before use.

```
cellxgene-schema annotate-perturbations input.h5ad output.h5ad
```

The command will extract guide sequences from `genetic_perturbations`, run guidescan2 to find genomic matches, identify overlapping genes using bioframe, and update the h5ad file with the annotations.

---

The data portal runs the following in the backend:

```
cellxgene-schema validate --add-labels output.h5ad input.h5ad
```

This execution validates the dataset as above AND adds the human-readable labels for the ontology and gene IDs as defined in the schema. If the validation is successful, a new AnnData file (output.h5ad) is written to disk with the labels appended.

This option SHOULD NOT be used by data contributors.

## Contributing

Please read our contributing [guidelines](CONTRIBUTING.md) and make sure adhere to the Contributor Covenant [code of conduct](https://github.com/chanzuckerberg/.github/blob/master/CODE_OF_CONDUCT.md). 

## Reporting Security Issues                     
                                                
Please read our [security reporting policy](SECURITY.md)

## Code of Conduct

This project adheres to the Contributor Covenant [code of conduct](https://github.com/chanzuckerberg/.github/blob/master/CODE_OF_CONDUCT.md).
By participating, you are expected to uphold this code. 
Please report unacceptable behavior to [opensource@chanzuckerberg.com](mailto:opensource@chanzuckerberg.com).
