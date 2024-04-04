# Updating reference files

## Quick start 
To update the gene and ontology references the following files need to be modified with the updated references:

- `cellxgene_schema_cli/cellxgene_schema/gencode_files/gene_info.yml`
- Bump cellxgene-ontology-guide version to latest release with updated ontologies in `cellxgene_schema_cli/cellxgene_schema/requirements.txt`

- Then the following will command download and process the reference files, and then update the respective package files under `cellxgene_schema_cli/cellxgene_schema/gencode_files/`

```bash
make update-references
```

The reference files are used internally by the class `cellxgene-schema.ontology.GeneChecker`

## Updating ontologies

Ontology information for cellxgene-schema is derived from cellxgene-ontology-guide. To update the ontologies in 
conjunction with a schema release, update [ontology_info.json in cellxgene-ontology-guide](https://github.com/chanzuckerberg/cellxgene-ontology-guide/blob/main/ontology-assets/ontology_info.json),
and then generate a new release with the subsequently generated ontology reference files. Bump the version of cellxgene-ontology-guide in 
`cellxgene_schema_cli/requirements.txt` to this latest release.

## Updating gene references

Gene references are stored under  `cellxgene_schema_cli/cellxgene_schema/gencode_files/`

The following gene references are currently stored (to see specific versions refer to schema definition):

1. Humnan `genes_homo_sapines.csv.gz`
2. Mouse `genes_mus_musculus_sapines.csv.gz`
3. sars‑CoV‑2 `genes_sars_cov_2.csv.gz`
4. ERCC Spike-Ins `genes_ercc.csv.gz`

To update these gene references, a developer has to do run the following:

```bash
make gene_processing
```

This will:

1. Download the reference gene files.
2. For 10X data, uncompress it, get the GTF files and remove the rest.
3. At the end the following files will be in this folder:

    - `mus_musculus.gtf.gz`
    - `homo_sapiens.gtf.gz`
    - `sars_cov_2.gtf.gz`
    - `ercc.txt`

The reference GTF files for mouse and human are obtained from the 10X reference data. The GTF file for sars‑CoV‑2 is obtained from ENSEMBL. The ERCC IDs are obtained from ThermoFisher.
To download all these files the following has to be called:

If reference gene versions need to be updated then `Makefile` has to be changed accordingly under the rules corresponding to each organism.

## Cleanup

If instead of the rule `update-references` individual make rules were executed, then make sure to remove all unused downloads using:
```bash
make clean
```