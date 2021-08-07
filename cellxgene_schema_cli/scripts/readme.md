# Updating reference files

## Quick start 
To update the gene and ontology references the following files need to be modified with the updated references:

- Genes `cellxgene_schema_cli/Makfile` (see rules for each organism)
- Ontologies `cellxgene_schema_cli/cellxgene_schema/ontology_files/owl_info` 

Then the following will command download and process the reference files, and then update the respective package files under `cellxgene_schema_cli/cellxgene_schema/ontology_files/`

```bash
make update-references
```

The reference files are used internally by the classes `cellxgene-schema.ontology.GeneChecker` and 
`cellxgene-schema.ontology.OntologyChecker`, those classes get the references file paths from `cellxgene-schema.env`

For more detailed information please read the rest of this readme.

## Updating ontologies

Ontology information for cellxgene-schema is  stored in `cellxgene_schema_cli/cellxgene_schema/ontology_files/all_ontology.json.gz` 

To update these ontologies, a developer has to:

### 1. Modify owl references

The yaml file `cellxgene_schema_cli/cellxgene_schema/ontology_files/owl_info.yml` contains the current ontologies and versions used.

Foe example, for CL:

```yaml
CL:
  latest: 2021-06-21
  urls:
    2021-06-21: https://github.com/obophenotype/cell-ontology/raw/v2021-06-21/cl.owl
```

`latest` references the version currently used by cellxgene-schema. When updating an ontology, the new version has to 
be added to `urls` under the current date and `latest` has to be to point at the correct date. It would like look this:

```yaml
CL:
  latest: new_date
  urls:
    2021-06-21: https://github.com/obophenotype/cell-ontology/raw/v2021-06-21/cl.owl
    new_date: https://foo.owl
```

**Note:** The processing scripts currently support `.owl` and `.owl.gz` files

### 2. Download and parse ontologies

From this directory (`cellxgene_schema_cli/scripts`) the following script has to be called:

```bash
make download-ontologies
```

This will:

1. Download ontology owl files to `cellxgene_schema_cli/cellxgene_schema/ontology_files/`
2. Parse each owl, extract important information and store all of it into **one** file `cellxgene_schema_cli/cellxgene_schema/ontology_files/all_ontology.json`

## Updating gene references

Gene references are stored under  `cellxgene_schema_cli/cellxgene_schema/ontology_files/`

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

## Updating gene references

If  instead of the rule `update-references` individual make rules were executed, then make sure to remove all unused downloads using:
```bash
make clean
```
