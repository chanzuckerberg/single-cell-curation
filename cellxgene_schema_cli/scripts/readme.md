# Updating reference files

## Updating ontologies

Ontology information for cellxgene-schema is derived from `cellxgene-ontology-guide`. To update the ontologies in 
conjunction with a schema release, update [ontology_info.json in cellxgene-ontology-guide](https://github.com/chanzuckerberg/cellxgene-ontology-guide/blob/main/ontology-assets/ontology_info.json),
and then generate a new release with the subsequently generated ontology files. Manually bump the version of `cellxgene-ontology-guide` in 
`cellxgene_schema_cli/requirements.txt` to match the version of this latest release.

```
cellxgene-ontology-guide==1.4.2 # if pointing to
```

## Updating gene references
Per-species gene information is managed and stored in the  `cellxgene_schema_cli/cellxgene_schema/gencode_files/` directory. The source of truth for per-species gene versions is in the `gene_info.yml` file. The targeted versions should match those specified in the `schema/<cxg-schema-version>/schema.md` file.

To update these gene references, a developer must update the `gene_info.yml` file and run the following:
```bash
make update-gene-references
```

This will manage the downloading, unpacking, and extraction of per-species gene level information and write the results for each species seperately e.g. `genes_ercc.csv.gz` and `genes_homo_sapiens.csv.gz`.
