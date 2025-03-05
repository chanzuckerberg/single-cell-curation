# Updating reference files

## Quick start 
To update the gene references, update `cellxgene_schema_cli/cellxgene_schema/gencode_files/gene_info.yml` with new references.

To update ontologies, update the `cellxgene-ontology-guide` version in `cellxgene_schema_cli/requirements.txt`.

- Then the following command will download and process the gene references and update the respective package files under `cellxgene_schema_cli/cellxgene_schema/gencode_files/`

```bash
make update-gene-references
```

The gene reference files are used internally by the class `cellxgene-schema.gencode:GeneChecker`

## Updating ontologies

Ontology information for cellxgene-schema is derived from `cellxgene-ontology-guide`. To update the ontologies in 
conjunction with a schema release, update [ontology_info.json in cellxgene-ontology-guide](https://github.com/chanzuckerberg/cellxgene-ontology-guide/blob/main/ontology-assets/ontology_info.json),
and then generate a new release with the subsequently generated ontology files. Bump the version of `cellxgene-ontology-guide` in 
`cellxgene_schema_cli/requirements.txt` to this latest release.

## Updating gene references
Per-species gene information is managed and stored in the  `cellxgene_schema_cli/cellxgene_schema/gencode_files/` directory. The source of truth for per-species gene versions is in the `gene_info.yml` file, therein. 

To update these gene references, a developer must update the `gene_info.yml` file and run the following:
```bash
make update-gene-references
```

This will managed the downloading, unpacking, and extraction of per-species gene level information and write the results for each species seperately e.g. `genes_ercc.csv.gz` and `genes_homo_sapiens.csv.gz`.
